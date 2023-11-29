#!/bin/bash

#PBS -N trim
#PBS -A open
#PBS -l nodes=1:ppn=4
#PBS -l pmem=4gb
#PBS -j oe
#PBS -l walltime=2:00:00
#PBS -m e
#PBS -M sjc6663@psu.edu

## run trimmomatic for length & quality
conda activate bioinfo
set -uex

# working directory
WORKDIR=/gpfs/group/evk5387/default/082022_MO7914_Core/usda

# input directory
INDIR=/gpfs/group/evk5387/default/082022_MO7914_Core/usda/fastq

# output directory
OUTDIR=/gpfs/group/evk5387/default/082022_MO7914_Core/usda/trimmed
mkdir -p $OUTDIR

# run trimmomatic
cat $WORKDIR/ids.txt | parallel ~/miniconda3/envs/bioinfo/bin/trimmomatic PE $INDIR/{}.trimmed_1.fastq $INDIR/{}.trimmed_2.fastq -baseout $OUTDIR/{}.trim.fastq SLIDINGWINDOW:4:20 MINLEN:100

# print finished message
echo "trimming complete"
