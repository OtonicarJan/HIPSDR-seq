#!/bin/bash

module load gcc/11.1.0
module load bowtie/1.3.0 

bin_kb=$1
REF=$2

bin=${bin_kb}000

mkdir -p refdata
cd refdata

## Mappability
if [ -f "hg38.150len.fa.map.bw" ]; then
    echo "hg38.150len.fa.map.bw  exists."
else 
    echo "hg38.150len.fa.map.bw  does not exist. Will be created"
    generateMap.pl --window 150 -o hg38.150len.fa.map.bw  -b $REF
    generateMap.pl --window 150 -o hg38.150len.fa.map.bw  $REF
fi

mkdir -p bin_${bin_kb}kb
cd bin_${bin_kb}kb

mapCounter -w ${bin}  -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY ../refdata/hg38.150len.fa.map.bw > hg38.${bin_kb}kb.map.seg

gcCounter -w ${bin}  -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $REF > hg38.${bin_kb}kb.gc.seg