### looking at psl after blat all UMI.fasta ###

fname = 'blat.valid.UMIs.psl'
x = read.delim(fname, skip=5, header=F)
names(x)=c('match','misMatch','repMatch','Ns','qGapCnt','qGapBases',
           'TgapCnt','TgapBases','strand','Qname','Qsize','Qstart',
           'Qend','Tname','Tsize','Tstart','Tend','blockCnt',
           'blockSizes','qStartMatch','TstartMatch')
x$umi = sapply(strsplit(x$Qname),'_'),'[',2)
umi = unique(x$umi)
result = 
for(u in umi){
  y = x[x$umi == u, ]
  ######################################### paused here 
}

parse_target_start = function(a){
  ## input
  ## a: psl record for one target on one strand
  a$count = as.numeric(sapply(strsplit(a$Qname,'count='),'[',2))
  b = NULL
  for(i in 1:nrow(a)){
    b = c(b, rep(a$Tstart[i], a$count[i]))
  }
  d = sort(table(b),decreasing = T); d
  p = as.numeric(names(d[1:2])); p = p[!is.na(p)]; 
  res = data.frame()
  for(d in p){
    b = a[a$Tstart == d, ]
    b = b[order(b$count, decreasing = T),]
    b$pct = round(b$count[1]/sum(b$count),4)
    res = rbind(res, b[1,])
  }
  res = res[res$pct>0.33, ]
  return(res)
}

parse_target = function(y){
  ## input
  ## y: psl record for one target region
  a = table(y$Qname)
  b = names(a[a>=2])
  if(length(b) == 0){return(data.frame())}
  if(length(b) == 1){
    y$count = as.numeric(unlist(strsplit(y$Qname,'count='))[2])
    y$pct = 1
    return(y)
  }
  y = y[y$Qname %in% b,]
  a = parse_target_start(y[y$strand == '+',]); a
  b = parse_target_start(y[y$strand == '-',]); b
  if(a$Qname != b$Qname){
    print('conflict call on PE, inspect please')
    return(data.frame())
  }
  return(rbind(a,b))
}

parse_UMI_psl = function(umi){
  ### load blat psl result ###

 
  ### filter off unpaired reads ###
  x$mapSize = x$Tend - x$Tstart +1
  x = x[x$mapSize>80 & x$mapSize<120,]
  x = x[x$qGapBases<10,]
  a = table(x$Qname)
  x = x[x$Qname %in% names(a[a>=2]),]

  ### process each target ###
  targets = unique(x$Tname)
  result = data.frame()
  for(target in targets){
    y = x[x$Tname == target,]; y
    result = rbind(result, parse_target(y))
    print(target)
  }
  
}

