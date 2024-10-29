#!/bin/bash

module load cellranger-dna/1.1.0
module load samtools
module load gcc/11.1.0

path_samples=$1

samples=$(ls $path_samples)
binsize=1000
type="sciHIPSD"
# Iterate over each combinatorial index
for sample in $samples
do
  # At this point, each combinatorial index should be processed separately
  if [[ $sample == LFS* ]]; then
    cellsummary=$sample'/outs/extracted_barcodes.txt'

    num_cells=$(cat $cellsummary | wc -l)

    echo $sample' with '$num_cells' cells'
    

    for cell_numb in $(seq 1 $num_cells)
    do
        cell_id=$(sed $cell_numb'q;d' $cellsummary)
        segment_file=$sample'/readCount_filtered_bam/readcounts.'$binsize'kb.cell_bc_'$cell_id'.seg'

        if [ ! -f "$segment_file" ] || [ ! -s "$segment_file" ]; then
            num_jobs=$(bjobs -w | wc -l)
            echo $num_jobs
            if [[ $num_jobs -gt 999 ]]; then
                now=$(date +"%T")
                echo 'too many jobs submitted, going to sleep, Current time:'$now
                sleep 5m
            fi
            echo 'running cell '$cell_id
            bsub -R "rusage[mem=20GB]" -q long -n 6 -e job_err_out.err -o job_err_out.out process_bam.sh $cell_numb $cell_id $sample $binsize $type
        else
            echo 'Segment file for cell '$cell_id' exists and is not empty.'
        fi
    done
  fi
done
