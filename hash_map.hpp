#pragma once

#include <upcxx/upcxx.hpp>
#include "kmer_t.hpp"

struct HashMap {
  upcxx::global_ptr<upcxx::global_ptr<kmer_pair>> data;
  upcxx::global_ptr<uint64_t> owning_rank;
  upcxx::global_ptr<uint64_t> locks;
  upcxx::global_ptr<upcxx::global_ptr<int>> used;
  // std::vector <kmer_pair> data;

  size_t my_size;
  size_t per_node_size;

  size_t size() const noexcept;
  HashMap();
  HashMap(size_t size, int n_proc);

  // Most important functions: insert and retrieve
  // k-mers from the hash table.
  bool insert(const kmer_pair &kmer);
  bool find(const pkmer_t &key_kmer, kmer_pair &val_kmer);

  // Helper functions
  uint64_t compute_data_offset(uint64_t i);
  uint64_t compute_segment_index(uint64_t i);
  void add_data_segment(upcxx::global_ptr<kmer_pair> segment, upcxx::global_ptr<int> used_seg, int rank);
  void lock_acquire(uint64_t slot);
  void lock_release(uint64_t slot);
  // Write and read to a logical data slot in the table.
  void write_slot(uint64_t slot, const kmer_pair &kmer);
  kmer_pair read_slot(uint64_t slot);

  // Request a slot or check if it's already used.
  bool request_slot(uint64_t slot);
  bool slot_used(uint64_t slot);
};

HashMap::HashMap() {

}

HashMap::HashMap(size_t size, int n_proc) {
  my_size = size;
  per_node_size = size / n_proc;
  data = upcxx::new_array<upcxx::global_ptr<kmer_pair>>(n_proc);
  used = upcxx::new_array<upcxx::global_ptr<int>>(n_proc);
  locks = upcxx::new_array<uint64_t>(n_proc);
  owning_rank = upcxx::new_<uint64_t>(0.0);
}

void HashMap::add_data_segment(upcxx::global_ptr<kmer_pair> segment, upcxx::global_ptr<int> used_seg, int rank) {
  uint64_t seg_idx = upcxx::rget(owning_rank).wait();
  upcxx::rput(segment, data + seg_idx).wait();
  upcxx::rput(used_seg, used + seg_idx).wait();
  uint64_t new_rank = rank + 1.0;
  upcxx::rput(new_rank, owning_rank).wait();
}

uint64_t HashMap::compute_data_offset(uint64_t i) {
  return i / per_node_size;
}

uint64_t HashMap::compute_segment_index(uint64_t i) {
  if (i > size()) {
    return -1;
  }
  // std::cout << "For node " << upcxx::rank_me() << ", Compute segment index:" << i << ", data offset: " << compute_data_offset(i) << ", index into local segment: " << i - (compute_data_offset(i) * per_node_size) << ", per_node_size: " << per_node_size << "\n" << std::flush;
  return i - (compute_data_offset(i) * per_node_size);
}

void HashMap::lock_acquire(uint64_t slot) {
  bool l_acquired = false;
  while (!l_acquired) {
    if (upcxx::atomic_get(locks + compute_data_offset(slot), std::memory_order_seq_cst).wait() == 0) {
      // acquire lock on the segment
      uint64_t acquire_val = 1.0;
      upcxx::atomic_put(locks + compute_data_offset(slot), acquire_val, std::memory_order_seq_cst).wait();
      return;
    }
  }
}

// should only ever be called by the owner of the lock.
void HashMap::lock_release(uint64_t slot) {
  // release the lock on your segment.
  uint64_t release_val = 0.0;
  upcxx::rput(release_val, locks + compute_data_offset(slot)).wait();
}

// ADD DOWNCASTING
bool HashMap::insert(const kmer_pair &kmer) {
  uint64_t hash = kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    success = request_slot(slot);
    if (success) {
      write_slot(slot, kmer);
    }
  } while (!success && probe < size());
  return success;
}

bool HashMap::find(const pkmer_t &key_kmer, kmer_pair &val_kmer) {
  uint64_t hash = key_kmer.hash();
  uint64_t probe = 0;
  bool success = false;
  do {
    uint64_t slot = (hash + probe++) % size();
    if (slot_used(slot)) {
      val_kmer = read_slot(slot);
      // std::cout << val_kmer.kmer.hash() << " : " << hash << "\n" << std::flush;
      if (val_kmer.kmer == key_kmer) {
        success = true;
      }
    }
  } while (!success && probe < size());
  // std::cout << "For node " << upcxx::rank_me() << ", insert succeeded. " << "\n" << std::flush;
  return success;
}

bool HashMap::slot_used(uint64_t slot) {
  upcxx::global_ptr<int> segment_used_array = upcxx::rget(used + compute_data_offset(slot)).wait();
  // this should be an atomic get?
  int used_val = upcxx::rget(segment_used_array + compute_segment_index(slot)).wait();
  if (used_val != 0) {
    return true;
  } else {
    return false;
  }
}

void HashMap::write_slot(uint64_t slot, const kmer_pair &kmer) {
  upcxx::global_ptr<kmer_pair> segment_array = upcxx::rget(data + compute_data_offset(slot)).wait();
  upcxx::rput(kmer, segment_array + compute_segment_index(slot)).wait();
  lock_release(slot);
  // std::cout << "For node " << upcxx::rank_me() << ", write succeeded. " << "\n" << std::flush;
}

kmer_pair HashMap::read_slot(uint64_t slot) {
  upcxx::global_ptr<kmer_pair> segment_array = upcxx::rget(data + compute_data_offset(slot)).wait();
  kmer_pair k = upcxx::rget(segment_array + compute_segment_index(slot)).wait();
  return k;
}

bool HashMap::request_slot(uint64_t slot) {
  // std::cout << "For node " << upcxx::rank_me() << " Requesting Slot:" << slot << "\n" << std::flush;
  upcxx::global_ptr<int> segment_used_array = upcxx::rget(used + compute_data_offset(slot)).wait();
  lock_acquire(slot);
  // std::cout << "For node " << upcxx::rank_me() << " Lock acquired. Slot:" << slot << "\n" << std::flush;
  // int* used_array = segment_used_array.local();
  // std::cout << "For node " << upcxx::rank_me() << " Segment used obtained. Slot:" << slot << "\n" << std::flush;
  // std::cout << "For node " << upcxx::rank_me() << " Segment Index: " << compute_segment_index(slot) << ", Slot: " << slot << "\n" << std::flush;
  int used_val = upcxx::rget(segment_used_array + compute_segment_index(slot)).wait();
  // std::cout << "For node " << upcxx::rank_me() << " Used Val: " << used_val << ", Slot:" << slot << "\n" << std::flush;
  if (used_val != 0) {
    lock_release(slot);
    return false;
  } else { 
    upcxx::rput(1, segment_used_array + compute_segment_index(slot)).wait();
    // std::cout << "For node " << upcxx::rank_me() << " Slot granted: " << slot << "\n" << std::flush;
    return true;
  }
}

size_t HashMap::size() const noexcept {
  return my_size;
}
