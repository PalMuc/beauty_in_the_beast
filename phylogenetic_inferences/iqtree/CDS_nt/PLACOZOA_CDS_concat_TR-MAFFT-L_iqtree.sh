#!/bin/bash

source activate iqtree

NUMTH=12
IN=PLACOZOA_CDS_concat_TR-MAFFT-L.fasta


# IQTREE
iqtree2 -s $IN -T $NUMTH -m GTR+I+G4 -b 1000 -o C0-H0__MH682141


conda deactivate