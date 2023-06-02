#!/bin/bash


source activate genomics


# pre-sets:
NUMTH=20
MEM=300
SPECIES=Helminthope_sp
RUNMT=SRR12246747
RUNWGS=SRR12246798


L1R1=${RUNMT}_1.fastq
L1R2=${RUNMT}_2.fastq
L2R1=${RUNWGS}_1.fastq
L2R2=${RUNWGS}_2.fastq



## Lib-1
trimmomatic PE -threads $NUMTH -phred33 \
$L1R1 $L1R2 \
${L1R1%.fastq}.cl.fastq ${L1R1%.fastq}.cl.s.fastq \
${L1R2%.fastq}.cl.fastq ${L1R2%.fastq}.cl.s.fastq \
ILLUMINACLIP:/home/ubuntu/tools/anaconda3/envs/genomics/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:3:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:6:20 \
MINLEN:55 \
2> ${SPECIES}_Lib1-mt_trimmomatic.log


## Lib-2
trimmomatic PE -threads $NUMTH -phred33 \
$L2R1 $L2R2 \
${L2R1%.fastq}.cl.fastq ${L2R1%.fastq}.cl.s.fastq \
${L2R2%.fastq}.cl.fastq ${L2R2%.fastq}.cl.s.fastq \
ILLUMINACLIP:/home/ubuntu/tools/anaconda3/envs/genomics/share/trimmomatic-0.39-1/adapters/TruSeq3-PE-2.fa:3:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:6:20 \
MINLEN:55 \
2> ${SPECIES}_Lib2-WGS_trimmomatic.log



# KARECT error correction:
karect -correct  \
-threads=$NUMTH \
-celltype=diploid \
-matchtype=hamming \
-inputfile=${L1R1%.fastq}.cl.fastq \
-inputfile=${L1R2%.fastq}.cl.fastq \
-inputfile=${L2R1%.fastq}.cl.fastq \
-inputfile=${L2R2%.fastq}.cl.fastq \
-inputfile=${L1R1%.fastq}.cl.s.fastq \
-inputfile=${L1R2%.fastq}.cl.s.fastq \
-inputfile=${L2R1%.fastq}.cl.s.fastq \
-inputfile=${L2R2%.fastq}.cl.s.fastq

gzip karect*fastq

rm *.cl.fastq
rm *.cl.s.fastq



# SPAdes assembly:
cat karect*_1.cl.fastq.gz > ${SPECIES}_gDNA_R1.cl.karect.fastq.gz
cat karect*_2.cl.fastq.gz > ${SPECIES}_gDNA_R2.cl.karect.fastq.gz
cat karect*_1.cl.s.fastq.gz *_2.cl.s.fastq.gz > ${SPECIES}_gDNA.cl.karect.s.fastq.gz

# mt-enriched reads:
spades.py --meta -k 33,55,77 -m $MEM -t $NUMTH -o ${SPECIES}_metaSPAdes_mito --only-assembler --disable-rr \
-1 karect_${L1R1%.fastq}.cl.fastq.gz \
-2 karect_${L1R2%.fastq}.cl.fastq.gz \
-s karect_${L1R1%.fastq}.cl.s.fastq.gz \
-s karect_${L1R2%.fastq}.cl.s.fastq.gz


spades.py -k 33,55,77 -m $MEM -t $NUMTH -o ${SPECIES}_SPAdes_mito --only-assembler --disable-rr \
-1 karect_${L1R1%.fastq}.cl.fastq.gz \
-2 karect_${L1R2%.fastq}.cl.fastq.gz \  
-s karect_${L1R1%.fastq}.cl.s.fastq.gz \
-s karect_${L1R2%.fastq}.cl.s.fastq.gz



# all reads:
READS1=${SPECIES}_gDNA_R1.cl.karect.fastq.gz
READS2=${SPECIES}_gDNA_R2.cl.karect.fastq.gz
READSS=${SPECIES}_gDNA.cl.karect.s.fastq.gz


spades.py --meta -k 33,55,77 -m $MEM -t $NUMTH -o ${SPECIES}_metaSPAdes_all --only-assembler --disable-rr \
-1 $READS1 -2 $READS2 -s $READSS

spades.py -k 33,55,77 -m $MEM -t $NUMTH -o ${SPECIES}_SPAdes_all --only-assembler --disable-rr \
-1 $READS1 -2 $READS2 -s $READSS

conda deactivate