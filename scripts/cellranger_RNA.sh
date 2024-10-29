#!/bin/bash


module load cellranger-7.1.0

sample_id=$1
reference=$2
fastqs=$3
samples=$4

cellranger count --id $sample_id \
--transcriptome $reference \
--fastqs $fastqs \
--sample $samples \
--localcores 50 --localmem 90 --chemistry ARC-v1