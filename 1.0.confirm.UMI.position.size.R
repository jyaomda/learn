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
blat panel.3.for.blat.fasta ../YY-ctDNA-MCL_PM1304_189/Sample_1224172/10K.record.fasta blat.10K.1224172.psl
table(a$size)
  81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96   97   98   99  100
  25   24   17   28   17   33   18   23   24   19   36   41   41   27 2332 2993 1018  281  100  789
## this means although the UMI tag is 3-bp at 5'-end, some has additional matches (98-100),
## majority are actually having additional 1-2bp unmatching sequences before target sequence
## nevertheless, the right way is to extract raw.R1/2.fastq then combine the 3-bp UMI tag
