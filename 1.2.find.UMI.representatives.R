#### processing psl and find UMI representatives ####
### input: psl files w/in umi_files folder
##  prodedure
##  for each umi collect records
##    limit mapping size {80,120} & gap size {0,10}
##    remove singleton dogs
##    for each target region tagged by this umi
##      find significant starts (>33% of population) on two strands
##      collect seq reads having starts on both strands
##      allow two umi tagging only, allow 1-bp start shift
##  merge result in a big table
##  adjust alignment start position !
##  calculate median extended count per start per target region
##    this is used to estimate coverage (amplifcation of each UMI)
##    find cutoff exteanded read to remove non-amplified UMIs
##  filter representative UMI table
##    remove UMIs having reads below extended read cutoff
##  the result UMI table is ready for pileup

folder = 'umi_files'

fname = 'blat.valid.UMIs.psl'
psl = read.delim(fname, skip=5, header=F)
names(psl)=c('match','misMatch','repMatch','Ns','qGapCnt','qGapBases',
           'TgapCnt','TgapBases','strand','Qname','Qsize','Qstart',
           'Qend','Tname','Tsize','Tstart','Tend','blockCnt',
           'blockSizes','qStartMatch','TstartMatch')
psl$umi = sapply(strsplit(x$Qname,'_'),'[',2)
result = data.frame()
for(u in umi){
  
  x = psl[psl$umi == u, ]
  
  ### filter off unpaired reads ###
  x$mapSize = x$Tend - x$Tstart +1
  x = x[x$mapSize>80 & x$mapSize<120,]
  x = x[x$qGapBases<10,]
  x = rm_singleton(x, 'Qname')
  x = rm_singleton(x, 'Tname')
  
  ### process each target ###
  targets = unique(x$Tname)
  res = data.frame()
  for(target in targets){
    y = x[x$Tname == target,]
    y = rm_singleton(y, 'Qname')
    if(nrow(y) == 0){ next }
    if(length(unique(y$strand)) == 1){ next }
    res = rbind(res, parse_target(y))
  }
  result = rbind(result, res)
  print(u)
}
result$extCount = round(result$count/result$pct)
##### testing area #####
target = 'chr11:118500936-118507669'
y = x[x$Tname == target,]; dim(y)
parse_target(y)
a = y[y$strand=='-',]
##### end of testing area #####

rm_singleton = function(x, colName){
  a = table(x[[colName]])
  return(x[x[[colName]] %in% names(a[a>=2]),])
}

parse_target = function(y){
  ## input
  ## y: psl record for one target region
  ## procedure
  ## find significant starts (>33% of population) on two strands
  ## collect seq reads having starts on both strands
  
  ### work on each strand 
  a = parse_target_start(y[y$strand == '+',])
  b = parse_target_start(y[y$strand == '-',])
  
  ### match result from each strand
  common = intersect(a$Qname, b$Qname)
  if(length(common)==0){return(data.frame())}
  z = data.frame()
  for(p in common){
    z = rbind(z, a[a$Qname == p,])
    z = rbind(z, b[b$Qname == p,])
  }
  return(z)
}

parse_target_start = function(a){
  ## input
  ## a: psl record for one target on one strand
  ## find significant starts (>33% of population) on this strand
  
  conditioned_two = function(d){
    ## input: sorted table() result
    ## procedure
    ## pick one
    ## remove, if any, near the first 1
    ## if remaining, pick another one
    pick = names(d)[1]
    neighbor = as.character(c(as.numeric(pick)+1, 
                              as.numeric(pick)-1))
    e = setdiff(d, c(pick, neighbor))
    if(length(e>0)){ pick = c(pick, names(e)[1]) }
    return(pick)
  }
  a$count = as.numeric(sapply(strsplit(a$Qname,'read='),'[',2))
  b = NULL
  for(i in 1:nrow(a)){
    b = c(b, rep(a$Tstart[i], a$count[i]))
  }
  d = sort(table(b),decreasing = T); d
  p = as.numeric(conditioned_two(d))
  res = data.frame()
  for(d in p){
    b = a[abs(a$Tstart - d) <= 1, ]
    ### this collects reads having 1-bp shift from start
    ## important to determine total reads from this start
    b = b[order(b$count, decreasing = T),]
    b$pct = round(b$count[1]/sum(b$count),4)
    res = rbind(res, b[1,])
  }
  res = res[res$pct>0.33, ]
  return(res)
}
