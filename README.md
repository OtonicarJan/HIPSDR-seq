# Scripts and notebooks for 'HIPSD&R-seq enables scalable genomic copy number and transcriptome profiling'

The analysis for the manuscript was performed in the following steps:

## 1. **Data preprocessing with cellranger pipelines**

   - For RNA component of HIPSD&R-seq, `scripts/cellranger_RNA.sh` was used.
   - For DNA component of HIPSD&R-seq, `scripts/cellranger_DNA.sh` was used.
   - For HIPSD-seq, sciHIPSD-seq and scATAC-seq, `scripts/cellranger_atac.sh` was used. Each combinatorial index was run separately in sciHIPSD-seq.

## 2. **Cell definition**

   - **RNA part**:\
     The cells were defined by cellranger pipeline.
   - **DNA part**:\
     The standard cellranger ATAC pipeline can not define cells from HIPSD-seq output as it bases the cell selection based on peak enrichment. Therefore, selection based on total fragments is used. For the paper we used cellranger's (v7.2.0) algorithm to call for the cells above the highest gradient of barcode rank plot. However, due to licensing limitations, we cannot share the code here. Thus, we provide our own cell calling algorithm. See `cell_calling/cell_definition.ipynb` for an example. With our implemetation we see a minor difference in the number of cells called compared to cellranger's version (1 cell difference for DEFND-seq and 11 cell difference for HIPSD&R-seq).

   All steps from herein are applied only to DNA part.

## 3. **BAM file postprocessing**

   - After you extract barcode associated with real cells, you can extract BAM files corresponding to single barcodes. Each cell BAM file should also should be filtered, sorted and indexed. We provide an example of postprocessing script in     `scripts/process_bam.sh`. To run postprocessing for all the desired barcodes, you can run `scripts/batch_submission.sh` for HIPSD or `scripts/batch_submission_sciHIPSD.sh` (but keep in mind that it's optimised for LSF batch scheduler). 

   - We provide a snakemake pipeline for BAM postprocessing: [HIPSDR-pipeline](https://github.com/OtonicarJan/HIPSDR-pipeline).

   - For step 5 you will also need genome mappability and GC files. These files can also be generated using [hmmcopy_utils](https://github.com/shahcompbio/hmmcopy_utils).
    We provide an example in `scripts/mappability_gc.sh`. If you plan to use single cells for the analysis, we suggest that you use 1000kb as a bin size. If you plan to use metacelling, use 100kb.


## 4. **Count processing (and optional metacells creation)**

  - After you obtain individual single-cell count seg files, we can perform some QC, i.e. number of counts per cell, number of non-empty bins per cell, TSS score enrichment per cell and fragment length distribution. Once you perform cell filtering, we can call CNAs. Optionally, you can perform metacelling at this stage:
    
  ### Metacelling (optional step that enables better resolution)

  - Metacells are produced based on a greedy algorithm. Please see the full notebook `02_HIPSDR_DNA_metacelling_DEMO.ipynb` for details. To execute the notebook, one would need to have the following files:

  - The result of the metacelling is a directory with `.txt` files. In each file there is a path to readcount files that are assigned to the metacell.

## 5. **Copy number calling**

  - Copy number colling is performed using [hmmcopy](https://bioconductor.org/packages/release/bioc/html/HMMcopy.html) R package. 

  - We developed a wrapper script that can be used to reproduce the analysis precisly. The script `scripts/run_hmmcopy.R` and it can be used as follows:

  - First create a file with paths to seg files `find path/to/counts/*.seg -type f > file_with_seg_counts`. Note that the script also uses bedtools.

  - Next, run:
  ```bash
  Rscript run_hmmcopy.R -f file_with_seg_counts -o hmmcopy_results/ -r path/to/reference/gc_and_mappability -l -w 1000
  # or for metacells
  for file in metacell_files/*.txt; do  Rscript run_hmmcopy.R -f $file -o hmmcopy_results/ -r path/to/reference/gc_and_mappability -l -b; done
  ```

  - To view all available options, please run `Rscript run_hmmcopy.R -h`. The script produces `.bed` files for each cell at the desired bin size.
