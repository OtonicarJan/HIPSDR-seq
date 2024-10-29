#!/bin/bash


module load cellranger-atac-2.1.0

sample_id=$1
reference=$2
fastqs=$3

cellranger-atac count --id $sample_id \
--reference $reference \
--fastqs $fastqs --localmem 90 --localcores 60