#### aim: extract all UMI and find correct UMI sequences ####
##
##  procedure
##  for each fastq record pair
##    record sequence pair under a UMI
##  for each UMI
##    find the consensus sequences for all molecules

add_UMI_record = function(s1, s2){
  ## input:
  ## s1/2: PE sequence reads
  ## procedure:
  ## extract UMI
  ## if UMI in result
  ##   append s1/2 sequences (UMI trimmed)
  ## else
  ##   initiate UMI list in result
  ##   initiate s1/2 record
  size = nchar(s1)
  umi = paste0(substr(s1,1,3), substr(s2,1,3))
  if(umi %in% names(result)){
    result[[umi]]$s1 <<- c(result[[umi]]$s1, substr(s1, 4, size))
    result[[umi]]$s2 <<- c(result[[umi]]$s2, substr(s2, 4, size))
  } else {
    result[[umi]]$s1 <<- substr(s1, 4, size)
    result[[umi]]$s2 <<- substr(s2, 4, size)
  }
}

make_UMI_fasta = function(u){
  ## input:
  ## u: UMI record list of s1/s2 PE sequences w/ UMI trimmed
  ## procedure:
  ## concatenate s1/s2, make table
  ## for each unique s1/s2
  ##   write fasta >line with freq
  ##   write fasta sequence s1/s2
  ## return fasta string
  
  s12 = paste0(u$s1, u$s2)
  a = table(s12); f = NULL
  for(i in 1:length(a)){
    f = c(f, paste0('>seq_', i, '_count=', a[i]))
    f = c(f, names(a[i]))
  }
  return(f)
}

result = list()
read_block = 4000000; n = 0
fh1 = gzfile('raw.R1.fastq.gz', 'r')
fh2 = gzfile('raw.R2.fastq.gz', 'r')

while(1){
  x1 = readLines(fh1, read_block)
  if(length(x1) == 0){ break }
  x2 = readLines(fh2, read_block)
  seq1 = x1[seq(2, read_block, 4)]
  seq2 = x2[seq(2, read_block, 4)]
  for(i in 1:(read_block/4)){
    add_UMI_record(seq1[i], seq2[i])
  }
  n = n + 1
  print(paste0(n,'M PE reads processed ... ', length(result),' UMIs found'))
}
close(fh1); close(fh2)

dir.create('umi')
for(u in names(result)){
  savename = paste0('umi/', u, '.fasta')
  writeLines(make_UMI_fasta(result[[u]]), savename)
}

        
        
