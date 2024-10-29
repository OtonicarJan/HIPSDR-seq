library(epiAneufinder)

epiAneufinder(input="fragments_sorted.tsv", #Enter path to your fragments.tsv file or the folder containing bam files
              outdir="epianeufinder_1e6", #Path to the directory where results should be written 
              blacklist="hg38-blacklist.v2.bed", #Path to bed file that contains the blacklisted regions of your genome
              windowSize=1e6, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram of sample data", 
              ncores=12,
              minFrags=20000)
