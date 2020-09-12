binning<-function(v,nbin,range=c(0,1)){
  m<-length(v)
  size<-floor(m/nbin)
  o<-rank(v,ties.method = 'random')
  bsize<-rep(size,nbin)
  index<-sample(1:nbin,m-size*nbin)
  bsize[index]=bsize[index]+1
  b<-numeric(m)
  loc<-cumsum(bsize)
  locs<-c(0,loc)
  for (kk in 1:nbin){
    b[which(o>locs[kk] & o<=locs[kk+1])]=kk
  }
  cutpt1<-sort(v)[loc[1:(nbin-1)]]
  cutpt2<-sort(v)[loc[1:(nbin-1)]+1]
  cutpt<-cutpt1+(cutpt2-cutpt1)*stats::runif(nbin-1)
  s=c(range[1],cutpt)
  e=c(cutpt,range[2])
  list(b=b,s=s,e=e)
}
