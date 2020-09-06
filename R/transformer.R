transformer<-function(x,newx,marginal=F,type='constant'){
  if (is.vector(x)){
    n<-length(x)
    d<-1
    x<-matrix(x,ncol=1)
    newx<-matrix(newx,ncol=1)
  }else{
    n<-dim(x)[1]
    d<-dim(x)[2]
  }

  if (marginal | d<=1){
    N0<-n

    #transform to [0,1]
    xx<-x
    newxx<-newx
    eps<-runif(1)/n
    for (i in 1:d){
      alli<-c(x[,i],newx[,i])
      maxi<-max(alli)
      mini<-min(alli)
      xx[,i]<-(x[,i]-mini+eps)/(maxi-mini+2*eps)
      newxx[,i]<-(newx[,i]-mini+eps)/(maxi-mini+2*eps)
    }

    #determine intervals
    startm=matrix(0,nrow=n,ncol=d)
    endm=matrix(0,nrow=n,ncol=d)
    binm=matrix(0,nrow=n,ncol=d)

    for (i in 1:d){
      xxi<-xx[,i]
      o<-rank(xxi,ties.method = 'random')
      sxx<-sort(xxi)
      cutpt1<-sxx[1:(n-1)]
      cutpt2<-sxx[2:n]
      cutpt<-cutpt1+(cutpt2-cutpt1)*runif(n-1)
      s=c(0,cutpt)
      e=c(cutpt,1)
      binm[,i]=o
      startm[,i]=s
      endm[,i]=e
    }

    I=endm-startm
    mm=(startm+endm)/2

    #transformation
    #x
    phi=matrix(0,nrow=n,ncol=d)
    for (i in 1:n){
      for (j in 1:d){
        xxij=xx[i,j]
        js=binm[i,j]
        li=js
        st=startm[li,j]
        if (I[li,j]==0) xmij=1
        else xmij=(xxij-mm[li,j])/I[li,j]
        HH=H(xmij,type=type)
        phi[i,j]=(js-1+HH)/n
      }
    }

    #newx
    if (is.vector(newxx)){
      newn<-length(newxx)
      newxx<-matrix(newxx,ncol=1)
    }else{
      newn<-dim(newxx)[1]
    }

    new_phi=matrix(0,nrow=newn,ncol=d)
    for (i in 1:newn){
      for (j in 1:d){
        xxij=newxx[i,j]
        li=which((startm[,j]<xxij) & (endm[,j]>=xxij))
        js=li
        st=startm[li,j]
        if (I[li,j]==0) xmij=1
        else xmij=(xxij-mm[li,j])/I[li,j]
        HH=H(xmij,type=type)
        new_phi[i,j]=(js-1+HH)/n
      }
    }

    return(list(xt=phi,newxt=new_phi))
  }else{
    N0<-floor(n^(1/d))
    size<-rep(N0,d)
    sizeold<-size
    ix<-1
    while(prod(size)<n){
      sizeold<-size
      size[ix]<-N0+1
      ix=ix+1
    }

    if (d>1){
      size<-sizeold
      size<-sample(size,d)
    }

    cumsize<-c(1,cumprod(size))

    #transform to [0,1]^d
    xx<-x
    newxx<-newx
    eps<-runif(1)/n
    for (i in 1:d){
      alli<-c(x[,i],newx[,i])
      maxi<-max(alli)
      mini<-min(alli)
      xx[,i]<-(x[,i]-mini+eps)/(maxi-mini+2*eps)
      newxx[,i]<-(newx[,i]-mini+eps)/(maxi-mini+2*eps)
    }

    #determine intervals
    defm=matrix(0,nrow=cumsize[d+1],ncol=d)
    for (kk in 1:d){
      defm[,kk]<-rep(rep(1:size[kk],each=prod(size[(kk+1):d])),cumsize[kk])
    }
    startm=matrix(0,nrow=cumsize[d+1],ncol=d)
    endm=matrix(0,nrow=cumsize[d+1],ncol=d)
    binm=matrix(0,nrow=n,ncol=d)
    bins=rep(1,n)
    BinInterval<-list()
    for (i in 1:d){
      start=numeric(cumsize[i+1])
      end=numeric(cumsize[i+1])
      newbins=numeric(n)
      for (j in 1:cumsize[i]){
        who=which(bins==j)
        use=xx[who,i]
        res=binning(use,size[i])
        newbins[who]=res$b
        binm[who,i]=res$b
        start[((j-1)*size[i]+1):(j*size[i])]=res$s
        end[((j-1)*size[i]+1):(j*size[i])]=res$e
      }
      BinInterval[[i]]<-cbind(start,end)
      bins=as.factor(bins):as.factor(newbins)
      levels(bins)=1:cumsize[i+1]
      startm[,i]=rep(start,each=prod(size[(i+1):d]))
      endm[,i]=rep(end,each=prod(size[(i+1):d]))
    }

    I=endm-startm
    mm=(startm+endm)/2
    Q=apply(I,1,cumprod)
    revcumsize=c(rev(cumprod(rev(size))[1:(d-1)]),1)

    #transformation
    #x
    phi=matrix(0,nrow=n,ncol=d)
    for (i in 1:n){
      js=binm[i,]
      if (d>1) li=sum((js-1)*revcumsize)+1
      else li=js
      st=startm[li,]
      xxi=xx[i,]
      for (j in 1:d){
        if (I[li,j]==0) xmij=1
        else xmij=(xxi[j]-mm[li,j])/I[li,j]
        HH=H(xmij,type=type)
        phi[i,j]=as.numeric((js[j]-1+HH)/size[j])
      }
    }

    #newx
    if (is.vector(newxx)){
      newn<-length(newxx)
      newxx<-matrix(newxx,ncol=1)
    }else{
      newn<-dim(newxx)[1]
    }

    new_phi=matrix(0,nrow=newn,ncol=d)
    new_binm=matrix(0,nrow=newn,ncol=d)
    new_bins=rep(1,newn)
    for (i in 1:d){
      new_start=BinInterval[[i]][,1]
      new_end=BinInterval[[i]][,2]
      new_newbins=numeric(newn)
      for (j in 1:cumsize[i]){
        who=which(new_bins==j)
        use=newxx[who,i]
        res=binning_search(use,new_start[((j-1)*size[i]+1):(j*size[i])],new_end[((j-1)*size[i]+1):(j*size[i])])
        new_binm[who,i]=res
        new_newbins[who]=res
      }
      new_bins=as.factor(new_bins):as.factor(new_newbins)
      levels(new_bins)=1:cumsize[i+1]
    }

    for (i in 1:newn){
      xxi=newxx[i,]
      js=new_binm[i,]
      li=sum((js-1)*revcumsize)+1
      st=startm[li,]
      for (j in 1:d){
        if (I[li,j]==0) xmij=1
        else xmij=(xxi[j]-mm[li,j])/I[li,j]
        HH=H(xmij,type=type)
        new_phi[i,j]=as.numeric((js[j]-1+HH)/size[j])
      }
    }

    return(list(xt=phi,newxt=new_phi))
  }
}
