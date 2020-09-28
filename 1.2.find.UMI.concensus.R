#### processing psl and find UMI representatives ####
### input: psl file within umi_files folder
##  prodedure
##  for each umi collect records
##    limit mapping size {80,120} & gap size {0,10}
##    remove singleton dogs
##    for each target region tagged by this umi
##      find significant starts (>33% of population) on two strands
##      collect seq reads having starts on both strands

#####  1. check amplication fold  #####
#### based on amplication fold set up cutoff for low UMIs
### which is 3 PCR cycle away from the median ampFold

get_density_peak = function(d, low_cut=log2(3)){
  e = data.frame(x = d$x, y = d$y)
  e = e[e$x>=low_cut, ]
  return(e$x[which.max(e$y)])
}

fnames = dir('umi_files', 'psl$', full.names=T); length(fnames)
medCnt = NULL; 
mapStat = data.frame(pos = 0:9, freq=0)
for(fname in fnames){
  x = read.delim(fname,skip=5,header=F)
  names(x)=c('match','misMatch','repMatch','Ns','qGapCnt','qGapBases',
             'TgapCnt','TgapBases','strand','Qname','Qsize','Qstart',
             'Qend','Tname','Tsize','Tstart','Tend','blockCnt',
             'blockSizes','qStartMatch','TstartMatch')
  y = unlist(strsplit(x$Qname, '_'))
  x$sid = y[seq(1,length(y),3)]
  x$read = as.numeric(gsub('read=','',y[seq(2,length(y),3)]))
  pos = match(unique(x$sid), x$sid)
  y = x[pos, c('sid','read')]
  d = density(log2(y$read))
  medCnt = c(medCnt, get_density_peak(d))
  a = table(x$Qstart)[1:10]
  b = as.data.frame(a)
  z = b$Freq[1:10]
  z[is.na(z)] = 0
  mapStat$freq = mapStat$freq + z
  print(gsub('.psl','',basename(fname)))
}
plot(density(medCnt))
x = get_density_peak(density(medCnt), 3); x
cutoff.uread = round(2^(x-3)); cutoff.uread
cutoff.xread = round(2^(x+3)); cutoff.xread
print(paste('median amplified UMI read =',round(2^x)))
print(paste('cutoff for low UMI read =', cutoff.uread))
print(paste('cutoff for high UMI read =', cutoff.xread))
barplot(mapStat$freq, main = 'mapping start positions', col='dark blue',
        names=mapStat$pos, ylab='frequency')
### this means for any UMI with same starting gene position
## UMI with 7+ reads will be kept and with 7 or less will be removed
# UMI with 463+ reads will be doubled in variant calling 

#####  2. find UMI targets  #####
rm_singleton = function(x, colName){
  a = table(x[[colName]])
  return(x[x[[colName]] %in% names(a[a>=2]),])
}
parse_target_start = function(a, cutoff){
  ## input a: psl record for one target on one strand
  ## find family of same start, merge by position +/-1
  ## return families with counts above cutoff
  
  conditioned_pick = function(k){
    ## input: sorted table() result
    ## procedure
    ## pick one
    ##   remove from list, if any, items near the first 1
    ## if remaining, pick another one
    ## do this recursively
    pick = names(k)[1]
    neighbor = as.character(c(as.numeric(pick)+1, 
                              as.numeric(pick)-1))
    e = setdiff(names(k), c(pick, neighbor))
    if(length(e>0)){
      pick = c(pick, conditioned_pick(k[e]))
    } else {
      return(pick) 
    }
    return(pick)
  }
  b = NULL
  for(i in 1:nrow(a)){
    b = c(b, rep(a$Tstart[i], a$read[i]))
  }
  d = sort(table(b),decreasing = T); d
  p = as.numeric(conditioned_pick(d))
  res = data.frame()
  for(d in p){
    b = a[abs(a$Tstart - d) <= 1, ]
    ### this collects reads having 1-bp shift from start
    ## consensus must be >33% of all species from this start
    # the UMI is kept when total reads > cutoff
    b = b[order(b$read, decreasing = T),]
    b$pct = round(b$read[1]/sum(b$read),4)
    res = rbind(res, b[1,])
  }
  res$xread = round(res$read/res$pct)
  res = res[res$pct > 0.33 & res$xread > cutoff,]
  return(res)
}
parse_target = function(y, cutoff){
  ## input y: psl record for one target region
  ## procedure
  ## find significant starts (>33% of population) on two strands
  ## collect seq reads having starts on both strands
  
  ### work on each strand, merge same start UMI family
  a = parse_target_start(y[y$strand == '+',],cutoff)
  b = parse_target_start(y[y$strand == '-',],cutoff)
  
  ### match result from each strand
  common = intersect(a$sid, b$sid)
  if(length(common)==0){return(data.frame())}
  z = data.frame()
  for(p in common){
    z = rbind(z, a[a$sid == p,])
    z = rbind(z, b[b$sid == p,])
  }
  return(z)
}

result = data.frame()
for(fname in fnames){
  u = gsub('.psl','',basename(fname))
  print(u)
  x = read.delim(fname,skip=5,header=F)
  names(x)=c('match','misMatch','repMatch','Ns','qGapCnt','qGapBases',
             'TgapCnt','TgapBases','strand','Qname','Qsize','Qstart',
             'Qend','Tname','Tsize','Tstart','Tend','blockCnt',
             'blockSizes','qStartMatch','TstartMatch')
  y = unlist(strsplit(x$Qname, '_'))
  x$sid = y[seq(1,length(y),3)]
  x$read = as.numeric(gsub('read=','',y[seq(2,length(y),3)]))
  x$mapSize = x$Tend - x$Tstart + 1
  
  ### initial filter ###
  x = x[x$mapSize>80 & x$mapSize<120 & x$qGapBases<10 & x$blockCnt<3,]

  ### process each target ###
  targets = unique(x$Tname)
  res = data.frame()
  for(target in targets){
    y = x[x$Tname == target,]
    y = rm_singleton(y, 'sid')
    if(nrow(y) == 0){ next }
    if(length(unique(y$strand)) == 1){ next }
    res = rbind(res, parse_target(y, cutoff.uread))
  }
  if(nrow(res) == 0){next}
  res$umi = u
  res$xRead = round(res$read/res$pct)
  res = res[res$xRead > cutoff.uread, ]
  result = rbind(result, res)
}
result$uid = paste(result$umi, result$sid, sep='_'); head(result,3)
length(unique(result$uid))  ## 24,317 for smalles lib

p = get_density_peak(density(log2(result$xRead),4))
cutoff.xread = round(2^(p+3)); cutoff.xread
range(result$xRead)
result[result$xRead> cutoff.xread,]

##### find insert size #####
a = sort(table(result$uid), decreasing=T)
result[result$uid == names(a)[1], ]
uid = names(a[a==2]); length(uid)
saveRDS(result, file='result.UMI.consensus.rds')
inserts = NULL
for(u in uid){
  b = result[result$uid == u,]
  inserts = c(inserts, max(b$Tend)-min(b$Tstart))
}
plot(density(log2(inserts)))
hist(log2(inserts), breaks=10000, xlim=c(6,log2(1000)))
hist(inserts, breaks=20000,xlim=c(0,600), main='Sample_1030163')
abline(v=c(1:3)*166,col='grey')
plot(density(log2(inserts)),xlim=c(6,10))
get_density_peak(density(log2(inserts)))
2^7.385    ## 167 basepair !
get_density_peak(density(log2(inserts)),8)
2^8.35847  ## 328 basepair !
get_density_peak(density(log2(inserts)),8.732)
2^8.887899 ## 474 basepair !


x =readRDS('result.UMI.consensus.rds'); dim(x)
a = table(x$uid); length(a[a>2]) ## 0.17% remove
torm = names(a[a>2]); torm
x = x[!x$uid %in% torm,]; dim(x)
table(x$TgapCnt)
table(x$blockCnt)
x = x[x$blockCnt == 1, ]
targets = unique(x$Tname); length(targets) ### 652 target regions


complement = function(sq){
  basepair = function(bp){
    if(bp == 'A'){return('T')}
    if(bp == 'C'){return('G')}
    if(bp == 'G'){return('C')}
    if(bp == 'T'){return('A')}
    if(bp == 'N'){return('N')}
  }
  a = unlist(strsplit(sq,'')); b = NULL
  for(n in rev(a)){ b=c(b, basepair(n))}
  paste(b,collapse='')
}
parse_fasta = function(fa){
  pos = grep('>', fa); length(pos)
  fstart = pos+1
  fend = c(pos[-1]-1,length(fa))
  z = list()
  for(i in 1:length(pos)){
    z[[gsub('>','',fa[pos[i]])]] = paste(fa[fstart[i]:fend[i]],collapse='')
  }
  return(z)
}
initialize_target_table = function(target, fa){
  ## input
  ##   target: target string
  ##   fa: fasta lines (targets w/ seq)
  
  ch = unlist(strsplit(target,':'))
  pos = as.numeric(unlist(strsplit(ch[2],'-'))); pos
  z = data.frame(chr=ch[1], pos=pos[1]:pos[2],
                 ref=unlist(strsplit(fa[[target]],'')),
                 A=0, C=0, G=0, T=0)
  return(z)
}
pileup = function(ttab, z){
  ## input
  ## ttab: target table to fill
  ## z: record from result having psl mapping
  fill_ttab = function(ttab, tstart, tend, sq){
    if(!(nchar(sq) == tend - tstart + 1)){
      print('seq length error'); return()
    }
    a = unlist(strsplit(sq,''))
    for(i in 1:length(a)){
      val = ttab[tstart+i-1, a[i]]
      ttab[tstart+i-1,a[i]] = val +1
    }
    return(ttab)
  }
  ### make an empty ttab copy ###
  d = ttab; d$A=0; d$C=0; d$G=0; d$T=0
  ### add plus strand to ttab copy ###
  a = z[z$strand=='+',]
  fau = parse_fasta(readLines(paste0('umi_files/', a$umi, '.fasta')))
  sq = fau[[a$Qname]]
  if(a$blockCnt==1){
    sq.yes = substr(sq, a$Qstart+1, a$Qend)
    d = fill_ttab(d, a$Tstart+1, a$Tend, sq.yes)
  }else{
    print('2 blocks')
  }
  ### add minus strand to ttab copy ###
  a = z[z$strand=='-',]
  # fau = parse_fasta(readLines(paste0('umi_files/', a$umi, '.fasta')))
  sq = fau[[a$Qname]];
  if(a$blockCnt==1){
    sq.rev = complement(substr(sq,a$Qstart+1,a$Qend))
    d = fill_ttab(d, a$Tstart+1, a$Tend, sq.rev)
  }else{
    print('2 blocks')
  }
  ### merge overlapping reads from short insert ###
  ## then add reads to ttab ##
  d[d==2] = 1
  ttab$A = ttab$A + d$A
  ttab$C = ttab$C + d$C
  ttab$G = ttab$G + d$G
  ttab$T = ttab$T + d$T
  return(ttab)
}

fa = readLines('panel.3.for.blat.fasta')
fa = parse_fasta(fa)
pile = data.frame()
for(target in targets){
  print(target)
  y = x[x$Tname == target,]
  ttab = initialize_target_table(target, fa)
  y = rm_singleton(y,'uid'); dim(y)
  for(u in unique(y$uid)){
    z = y[y$uid == u,]
    ttab = pileup(ttab, z)
  }
  # a = apply(ttab[,4:7],1,sum)
  # plot(a,type='h',main=target)
  # locator(1)
  ttab$depth = apply(ttab[,4:7],1,sum)
  ttab$target = target
  ttab = ttab[ttab$depth>2, ]
  pile = rbind(pile, ttab)
  print(paste0(round(match(target,targets)/length(targets)*100,1),'% done'))
}
saveRDS(pile, file='pileup.rds')
############### call mutation ##############
head(pile)
x = pile[pile$ref=='A' & (pile$C>0 | pile$G>0 | pile$T>0), ]; dim(x)
y = pile[pile$ref=='C' & (pile$A>0 | pile$G>0 | pile$T>0), ]; dim(y)
z = pile[pile$ref=='G' & (pile$C>0 | pile$A>0 | pile$T>0), ]; dim(z)
a = pile[pile$ref=='T' & (pile$C>0 | pile$G>0 | pile$A>0), ]; dim(a)

mut=rbind(x,y,z,a)
mut=mut[order(mut$pos),]
mut=mut[order(mut$chr),]
head(mut)
dim(mut)
saveRDS(mut, file='mut.rds')
hist(log2(pile$depth), breaks=100)
plot(density(log2(pile$depth)), col='blue')
lines(density(log2(mut$depth)), col='red')
abline(v=log2(cutoff.uread),col='grey')
mut= mut[mut$depth>cutoff.uread,]; dim(mut)
tail(mut)




