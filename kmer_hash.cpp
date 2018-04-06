#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <set>
#include <numeric>
#include <cstddef>
#include <upcxx/upcxx.hpp>

#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "hash_map.hpp"

#include "butil.hpp"

int main(int argc, char **argv) {
  upcxx::init();

  if (argc < 2) {
    BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test]\n");
    upcxx::finalize();
    exit(1);
  }

  std::string kmer_fname = std::string(argv[1]);
  std::string run_type = "";

  if (argc >= 3) {
    run_type = std::string(argv[2]);
  }

  int ks = kmer_size(kmer_fname);

  if (ks != KMER_LEN) {
    throw std::runtime_error("Error: " + kmer_fname + " contains " +
      std::to_string(ks) + "-mers, while this binary is compiled for " +
      std::to_string(KMER_LEN) + "-mers.  Modify packing.hpp and recompile.");
  }

  size_t n_kmers = line_count(kmer_fname);

  // Load factor of 0.5
  size_t hash_table_size = n_kmers * (1.0 / 0.5);
  size_t n_per_node = hash_table_size / upcxx::rank_n();
  // make the hash table evenly divisible by the number of nodes
  hash_table_size = n_per_node * upcxx::rank_n();
  upcxx::global_ptr<HashMap> hashmap_ptr = nullptr;
  if (upcxx::rank_me() == 0) {
    HashMap* h = new HashMap(hash_table_size, upcxx::rank_n());
    hashmap_ptr = upcxx::allocate<HashMap>(1);
    upcxx::rput(*h, hashmap_ptr).wait();
  }
  hashmap_ptr = upcxx::broadcast(hashmap_ptr, 0).wait();
  upcxx::global_ptr<kmer_pair> local_region = upcxx::new_array<kmer_pair>(n_per_node);
  upcxx::global_ptr<int> used_local_region = upcxx::new_array<int>(n_per_node);
  HashMap tmpHashmap = upcxx::rget(hashmap_ptr).wait();
  while (upcxx::rget(tmpHashmap.owning_rank).wait() != upcxx::rank_n()) {
    if (upcxx::rget(tmpHashmap.owning_rank).wait() == upcxx::rank_me()) {
      tmpHashmap.add_data_segment(local_region, used_local_region, upcxx::rank_me());
      upcxx::rput(tmpHashmap, hashmap_ptr).wait();
    }
  }
  upcxx::barrier();
  HashMap hashmap = upcxx::rget(hashmap_ptr).wait();
  upcxx::barrier();
  std::vector <kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

  if (run_type == "verbose") {
    BUtil::print("Finished reading kmers.\n");
  }

  auto start = std::chrono::high_resolution_clock::now();

  std::vector <kmer_pair> start_nodes;
  int j = kmers.size();
  for (auto &kmer : kmers) {
    bool success = hashmap.insert(kmer);
    // std::cout << "Node " << upcxx::rank_me() << ", number left: " << j << "\n" << std::flush;
    if (!success) {
      throw std::runtime_error("Error: HashMap is full!");
    }

    if (kmer.backwardExt() == 'F') {
      start_nodes.push_back(kmer);
    }
    j -= 1;
  }
  auto end_insert = std::chrono::high_resolution_clock::now();
  upcxx::barrier();

  double insert_time = std::chrono::duration <double> (end_insert - start).count();
  if (run_type != "test") {
    BUtil::print("Finished inserting in %lf\n", insert_time);
  }
  upcxx::barrier();

  auto start_read = std::chrono::high_resolution_clock::now();

  std::list <std::list <kmer_pair>> contigs;
  for (const auto &start_kmer : start_nodes) {
    std::list <kmer_pair> contig;
    contig.push_back(start_kmer);
    while (contig.back().forwardExt() != 'F') {
      kmer_pair kmer;
      bool success = hashmap.find(contig.back().next_kmer(), kmer);
      if (!success) {
        throw std::runtime_error("Error: k-mer not found in hashmap.");
      }
      contig.push_back(kmer);
    }
    contigs.push_back(contig);
  }

  auto end_read = std::chrono::high_resolution_clock::now();
  upcxx::barrier();
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration <double> read = end_read - start_read;
  std::chrono::duration <double> insert = end_insert - start;
  std::chrono::duration <double> total = end - start;

  int numKmers = std::accumulate(contigs.begin(), contigs.end(), 0,
    [] (int sum, const std::list <kmer_pair> &contig) {
      return sum + contig.size();
    });

  if (run_type != "test") {
    BUtil::print("Assembled in %lf total\n", total.count());
  }

  if (run_type == "verbose") {
    printf("Rank %d reconstructed %d contigs with %d nodes from %d start nodes."
      " (%lf read, %lf insert, %lf total)\n", upcxx::rank_me(), contigs.size(),
      numKmers, start_nodes.size(), read.count(), insert.count(), total.count());
  }

  if (run_type == "test") {
    std::ofstream fout("test_" + std::to_string(upcxx::rank_me()) + ".dat");
    for (const auto &contig : contigs) {
      fout << extract_contig(contig) << std::endl;
    }
    fout.close();
  }

  upcxx::finalize();
  return 0;
}
