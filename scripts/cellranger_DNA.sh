#!/bin/bash


module load cellranger-arc/2.0.2

sample_id=$1
reference=$2
libraries=$3

cellranger-arc count --id $sample_id \
--reference $reference \
--libraries $libraries --localmem 60 --localcores 50