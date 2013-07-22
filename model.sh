#!/bin/bash

for i in `seq 1 2`; do
echo "Rscript model.R $i" | qsub -d `pwd`
done;
