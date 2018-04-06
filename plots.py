test_dataset_1 = [(16.768505, 8.581043), (13.030488, 4.229887), (13.893945, 2.106382), (37.553939, 1.177379), (13.424993, 0.610240), (13.000536, 0.320153), (15.935454, 0.259916)]

test_dataset_2 = [(13.000536, 0.320153), (22.851817, 1.614905), (22.267222, 1.504), (31.382042, 1.559455)]
test_dataset_4 = []
test_dataset_8 = []

human_1 = [(340.755999, 180.518667), (266.658896, 93.697722)]
salloc -N 8 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 8 -n 32 ./kmer_hash $SCRATCH/my_datasets/test.txt verbose
salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C haswell srun -N 1 -n 2 ./kmer_hash $SCRATCH/my_datasets/human-chr14-synthetic.txt verbose