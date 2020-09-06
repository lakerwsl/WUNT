binning_search<-function(v,start,end){
  m<-length(v)
  b<-numeric(m)
  nbin=length(start)
  for (kk in 1:nbin){
    b[which(v>start[kk] & v<=end[kk])]=kk
  }
  return(b)
}
