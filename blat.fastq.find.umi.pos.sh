### extract fastq for target panel ###
## 1.1. prepare a chr:start-stop per line file for extraction
## 1.2. make fasta file for blat reference 
  samtools faidx /rsrch3/scratch/mol_cell_onc/jyao/work/gatk_bundle_hg38/original.fasta.copy \
  -r merged.3.panels.names.txt >panel.3.for.blat.fasta
  
### blat 1,000,000 fastq sequences ###
## 2.1. prepare fasta from 10,000 fastq records
##      do this for biggest file and smallest file
##      for high and low duplication ratios
  fh = gzfile('raw.R1.fastq.gz','r'); x = readLines(fh, 40000)
  close(fh); y = x[seq(2, 40000, 4)]; head(y)
  fh = file('10K.record.fasta', 'w')
  for(i in 1:10000){ 
    writeLines(paste0('>seq_',i),fh)
    writeLines(y[i], fh)}; close(fh)
    
### figure out where UMI seqs are ###
