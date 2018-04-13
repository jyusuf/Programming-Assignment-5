seq: sort.c
	mpicc -o sort sort.c

parA: sort_partA.c
	mpicc -o parAsort sort_partA.c
        
run-seq: seq
	srun -n 1 ./sort 16777216
	srun -n 1 ./sort 33554432
	srun -n 1 ./sort 67108864
	srun -n 1 ./sort 134217728

run-parA: parA
	srun -n 1 ./parAsort 134217728
	srun -n 2 ./parAsort 134217728
	srun -n 4 ./parAsort 134217728
	srun -n 8 ./parAsort 134217728
