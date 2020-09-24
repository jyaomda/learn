#### aim: extract all UMI and find correct UMI sequences ####
##
##  procedure
##  for each fastq record pair
##    record sequence pair under a UMI
##  for each UMI
##    find the consensus sequences for all molecules
##
##  results
##  182Mb/3M- records took 25 min, so larger 3,800Mb file may take 9 hours
##  182Mb/3M records v2 4 min, so larger 3.8Gb/56M records took 8 hours (OK)
##  for v2 w/ combined fasta, it takes 9.7 hours to write 5Gb fasta,
##    whereas blat of this 5Gb fasta to reference panel takes only 5 min
##    therefore blat is not the bottle neck step here, extremely fast!
##
##  NANA.fasta produced, with 350K (up to 10%), what are these?
##  NANA removed in v2, likely from last reading block
##
##  purity v1small v2small  v2large
##       0  384988   33522   675032 <- 1.2% constant 
##       1 2615012 2615012 54865344 
##  from here we know ambiguous UMI 1.2% can be removed
## 
##  even though v4 sounds cleaner, computational wise it is slower
##  the winning version is v3.1. in this case (case closed) ##

### version 3 ###
add_UMI_record = function(s1, s2){
  ## input
  ## s1/s2: PE reads
  ## procedure
  ##   only valid UMI are recorded
  ## note
  ## use size = 100 for this project
  ## assume all fastq reads are 100 nt in length
  
  ### if no read, bad UMI, pass ###
  if(is.na(s1)){ return() }
  umi = paste0(substr(s1,1,3), substr(s2,1,3))
  pos = grep('N', umi)
  if(length(pos)>0){ return() }
  
  ### use PE seq as list index, record freq ###
  pe = paste0(substr(s1,4,100), substr(s2,4,100))
  if(umi %in% names(result)){
    if(pe %in% names(result[[umi]])){
      result[[umi]][[pe]] <<- result[[umi]][[pe]] + 1
    }else{ result[[umi]][[pe]] <<- 1 }
  }else{
    result[[umi]] <<- list()
    result[[umi]][[pe]] <<- 1
  }
}

result = list()
read_block = 2000000; n = 0
fh1 = gzfile('raw.R1.fastq.gz', 'r')
fh2 = gzfile('raw.R2.fastq.gz', 'r')

while(1){
  x1 = readLines(fh1, read_block)
  x2 = readLines(fh2, read_block)
  if(length(x1) == 0){ break }
  seq1 = x1[seq(2, read_block, 4)]
  seq2 = x2[seq(2, read_block, 4)]
  for(i in 1:(read_block/4)){
    add_UMI_record(seq1[i], seq2[i])
  }
  n = n + 1
  print(paste0(n*(read_block/4000000),'M reads, ', length(result),' UMIs'))
}
close(fh1); close(fh2)

### version 3, just make fasta ###
fa_str = NULL
for(u in names(result)){
  z = result[[u]]
  for(i in 1:length(z)){
    fa_str = c(fa_str, paste0('>UMI_', u,'_seq', i,'_count=', z[[i]]))
    fa_str = c(fa_str, names(z)[i])
  }
}
writeLines(fa_str, 'valid.UMI.seq.v3.fasta')

########## evaluate results ###########
### with in umi folder 
fnames = dir('.', 'fasta$')
reads = NULL
for(f in fnames){
  x = readLines(f)
  y = x[seq(1, length(x), 2)]
  z = sapply(strsplit(y,'count='),'[',2)
  reads = c(reads, sum(as.numeric(z)))
}
names(reads) = fnames
a = as.data.frame(reads)
a$good = 1; pos = grep('N', row.names(a))
a[pos, 'good'] = 0
aggregate(a[,1], by=list(a$good), sum)

        
