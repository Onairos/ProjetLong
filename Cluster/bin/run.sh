#!/bin/bash

/mnt/hmf/mpich2/bin/mpirun -np $1 clusters $2 
#| tee output
ls cluster.final.* | wc -l > nbclusters

