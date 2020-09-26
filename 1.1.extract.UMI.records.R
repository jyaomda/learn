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
##  final v3.2 largeted 56M records took 11 hours to extract and blat
##    this is ultimately chosen over combined UMI processing below
## 
##  for v2 w/ combined fasta, it takes 9.7 hours to write 5Gb fasta,
##    whereas blat of this 5Gb fasta to reference panel takes only 5 min
##    therefore blat is not the bottle neck step here, extremely fast!
##    this branch is terminated because 4,096 files provide indexing
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

### version 3.3 as of 9/26/2020 ###

#################### functions for extracting valid UMIs ####################

add_UMI_record = function(s1, s2){
  ## input s1/s2: PE reads
  ## use size 100 for this project
  
  ### exit for NA str or bad UMI ###
  if(is.na(s1)){ return() }
  umi = paste0(substr(s1,1,3), substr(s2,1,3))
  pos = grep('N', umi)
  if(length(pos)>0){ return() }

  ### append PE to existing list ###
  pe = paste0(substr(s1,4,100), substr(s2,4,100))
  result[[umi]] <<- c(result[[umi]], pe)
}

make_UMI_fasta = function(u){
  ## input u: content of PE strings
  ## produce a fasta string:
  ##   >$seqNo_$readCount header, followed by
  ##   PE read1 + read2 (w/o UMI overhangs)
  
  a = table(u); f = NULL
  for(i in 1:length(a)){
    f = c(f, paste0('>s', i, '_read=', a[i]))
    f = c(f, names(a[i]))
  }
  return(f)
}

extract_UMI = function(){

  ### prepare file handles ###
  result = list()
  read_block = 4000000; n = 0
  fh1 = gzfile('raw.R1.fastq.gz', 'r')
  fh2 = gzfile('raw.R2.fastq.gz', 'r')

  ### read fastq by block ###
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
    print(paste0(n*(read_block/4000000),'M PE reads processed'))
  }
  close(fh1); close(fh2)

  ### write fasta/psl file by UMI ###
  dir.create('umi_files')
  for(u in names(result)){
    savename = paste0('umi_files/', u, '.fasta')
    writeLines(make_UMI_fasta(result[[u]]), savename)

    ## blat each fasta file, save psl ##
    cmd_line = 'blat ../../2020_MWang_MCL_PM1304/panel.3.for.blat.fasta'
    cmd_line = paste(cmd_line, savename)
    savename = gsub('.fasta','.psl',savename)
    cmd_line = paste(cmd_line, savename)
    system(cmd_line)
  }
}

#################### main procedure starts here ####################

if(file.exists('umi_files')){
  print('UMI already extracted, proceed to find representative UMIs ...')
} else {
  print('UMI have not been extracted, start extracting UMIs ...')
  extract_UMI()
}
