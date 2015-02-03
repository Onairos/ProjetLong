#!/bin/bash

mpirun -np $1 --hostfile ../../hosts/hosts clusters $2 | tee output
ls cluster.final.* | wc -l > nbclusters

