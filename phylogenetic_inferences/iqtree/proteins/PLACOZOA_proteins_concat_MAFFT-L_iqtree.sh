#!/bin/bash

source activate iqtree

NUMTH=12
IN=PLACOZOA_protein_concat_MAFFT-L.fasta


# IQTREE
iqtree2 -s $IN -T $NUMTH -m JTT+I+G4+F -b 1000 -o C0-H0__MH682141


conda deactivate
