#!/bin/bash

NUMTH=12
IN=PLACOZOA_mt-DNA_rRNA_FL_concat_MUSCLE.fasta
MODEL=GTR+I+G4

# IQTREE
iqtree2 -s $IN -T $NUMTH -m $MODEL -b 1000 -o C0-H0__MH682141
