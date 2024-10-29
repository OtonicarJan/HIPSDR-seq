#!/bin/bash

# This script was created specifically for LSF HPC system!
# Edit the code for other HPC type.

module load cellranger-dna/1.1.0
module load samtools
module load gcc/11.1.0

barcodes=$1
workdir=$2
binsize=$3
type="HIPSD"

num_cells=$(cat $barcodes | wc -l)

echo $num_cells' cells'

for cell_numb in $(seq 1 $num_cells)
do
    cell_id=$(sed $cell_numb'q;d' $barcodes)
    segment_file=$workdir'/readCount_filtered_bam/readcounts.'$binsize'kb.cell_bc_'$cell_id'.seg'
    if [ ! -f "$segment_file" ] || [ ! -s "$segment_file" ]; then
        num_jobs=$(bjobs -w | wc -l)
        echo $num_jobs
        if [[ $num_jobs -gt 999 ]]; then
            now=$(date +"%T")
            echo 'too many jobs submitted, going to sleep, Current time:'$now
            sleep 5m
        fi
        echo 'running cell '$cell_id
        bsub -R "rusage[mem=20GB]" -q long -n 6 -e job_err_out.err -o job_err_out.out process_bam.sh $cell_numb $cell_id $workdir $binsize $type
    else
        echo 'Segment file for cell '$cell_id' exists and is not empty.'
    fi
done
