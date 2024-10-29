#!/bin/bash
set -e

cell_numb=$1
cell_id=$2
working_dir=$3
bin_kb=$4
type=$5

bin=${bin_kb}000

mkdir -p $working_dir
cd $working_dir

# Subset to single-barcode BAM file
if [ -f "bam_sliced" ]; then
    echo "bam_sliced exists."
else
    if [ "$type" == "HIPSD" ]; then
        bam_path="../outs/atac_possorted_bam.bam"
    elif [ "$type" == "sciHIPSD" ]; then
        bam_path="../outs/possorted_bam.bam"
    else
        echo "Wrong type."
    fi
    echo "bam_sliced does not exist. Will be created"

    mkdir -p 'bam_sliced'
    cd 'bam_sliced'

    echo $cell_id > cell_id_$cell_numb.csv
    echo "library_id,barcodes_csv" > path_cell_id_$cell_numb.csv
    echo ${cell_id},cell_id_$cell_numb.csv >> path_cell_id_$cell_numb.csv
    
    cellranger-dna bamslice --id=${cell_id} \
                       --csv=path_cell_id_$cell_numb.csv \
                       --bam=$bam_path \
                       --localcores=6 \
                       --localmem=20  

    echo  "Done with cellranger bamslice on cell " $cell_numb " with barcode " cell_$cell_id
    rm cell_id_$cell_numb.csv
    rm path_cell_id_$cell_numb.csv
    cd ..
fi

# Postprocessing of the single-barcode BAM file
CLEAN_BAM_FOLDER='clean_bam'
FLAGSTAT_FOLDER='flagstat'
mkdir -p $CLEAN_BAM_FOLDER
mkdir -p $FLAGSTAT_FOLDER


CELL_BAM='bam_sliced/'$cell_id'/outs/subsets/'$cell_id'.bam'
CLEAN_BAM='clean_bam/'$cell_id'.filtered.bam'
CLEAN_SAM='clean_bam/'$cell_id'.filtered.sam'
CLEAN_SORT_BAM='clean_bam/'$cell_id'.filtered.sorted.bam'
CLEAN_SORT_RH_BAM='clean_bam/'$cell_id'.filtered.sorted.rh.bam'
HEADER='clean_bam/'$cell_id'.header.sam'

if [ -f "$CLEAN_SORT_BAM" ]; then
    echo "$CLEAN_SORT_BAM exists."
else 
    echo "$CLEAN_SORT_BAM does not exist. Will be created"
    samtools view -H $CELL_BAM > $CLEAN_SAM
    samtools view -F 3844 -q 30 $CELL_BAM >> $CLEAN_SAM
    samtools view -bS  $CLEAN_SAM > $CLEAN_BAM
    samtools sort $CLEAN_BAM -o $CLEAN_SORT_BAM
    samtools view -H $CLEAN_SORT_BAM > $HEADER
    sed -i "1s/.*/@HD\tVN:1.3\tSO:coordinate/" $HEADER
    samtools reheader $HEADER $CLEAN_SORT_BAM > $CLEAN_SORT_RH_BAM
    samtools index $CLEAN_SORT_RH_BAM
    samtools flagstat $CLEAN_SORT_BAM > 'flagstat/flagstat.'$cell_id'.filtered.txt'
    samtools flagstat $CELL_BAM > 'flagstat/flagstat.'$cell_id'.original.txt'
    rm $CLEAN_SAM 
    rm $CLEAN_BAM
    rm $HEADER
fi

# Prepare single-barcode wig file with counts per fixed-length genomic bin
READCOUNTS='readCount_filtered_bam'
mkdir -p $READCOUNTS
cd $READCOUNTS

SEGMENT_FILE=readcounts.${bin_kb}kb.cell_bc_${cell_id}.seg
if [ -f "$SEGMENT_FILE" ]; then
    echo "$SEGMENT_FILE exists."
else 
    echo "$SEGMENT_FILE does not exist. Will be created"
    ## Set up the path to readCounter
    readCounter -q 30 -w ${bin} -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY "../$CLEAN_SORT_RH_BAM" > $SEGMENT_FILE
fi
echo  "Done with readcCounter "${bin_kb}"kb window on cell "${cell_id}".bam file " 

