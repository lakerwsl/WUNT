wunt <- function(x,z,y,type='constant',marginal=FALSE,extra=FALSE,xextra=NULL,method='Kernel',h=NULL,kernel='gaussian',degree=2,interact=TRUE){
  if (extra) stopifnot(!is.null(xextra))
  #ATT
  if (is.vector(x)) x<-matrix(x,ncol=1)
  n=dim(x)[1]
  d=dim(x)[2]
  treat=z==1
  control=z==0
  if (extra==F){
    phi=transformer(x[control,],x[treat,],marginal,type)
    phi0=phi$xt
    phi1=phi$newxt
  }else{
    phi=transformer(xextra,x,marginal,type)$newxt
    phi0=phi[control,]
    phi1=phi[treat,]
  }

  if (is.vector(phi0)) phi0<-matrix(phi0,ncol=1)
  if (is.vector(phi1)) phi1<-matrix(phi1,ncol=1)
  y0=y[control]
  y1=y[treat]
  if (method=='Kernel'){
    if (is.null(h)) h=1/n^{1/d}
    tsum1=0
    tsum2=0

    weight=numeric(nrow(phi0))
    for(i in 1:nrow(phi0)){
      phi0i=matrix(phi0[i,],ncol=d,nrow=nrow(phi1),byrow = T)
      K=Kernel(phi0i-phi1,h,kernel)
      ww=sum(K)
      tsum1=tsum1+ww*y0[i]
      tsum2=tsum2+ww
      weight[i]=ww
    }

    if (tsum2==0) mcontrol=0
    else  mcontrol=tsum1/tsum2

    mtreat=mean(y1)

    Re=(mtreat-mcontrol)
    return(list(est=Re,weight=weight/tsum2,control_t=phi0,treated_t=phi1))
  }else if (method=='Projection'){

    tsum1=0
    tsum2=0
    mphi0=cbind(rep(1,sum(1-z)),phi0)
    mphi1=cbind(rep(1,sum(z)),phi1)

    if (!interact){
      it=1
      while(it<degree){
        it=it+1
        for (j in 1:d) {
          mphi0<-cbind(mphi0,phi0[,j]^it)
          mphi1<-cbind(mphi1,phi1[,j]^it)
        }
      }
    }else{
      if (degree>=2){
        for (j in 1:d) {
          mphi0<-cbind(mphi0,phi0[,j]*phi0[,j:d])
          if (sum(treat)==1) mphi1<-cbind(mphi1,matrix(phi1[,j]*phi1[,j:d],nrow=1))
          else mphi1<-cbind(mphi1,phi1[,j]*phi1[,j:d])
        }
        if (degree>=3){
          for (j1 in 1:d) {
            for (j2 in j1:d){
              mphi0<-cbind(mphi0,phi0[,j1]*phi0[,j2]*phi0[,j2:d])
              mphi1<-cbind(mphi1,phi1[,j1]*phi1[,j2]*phi1[,j2:d])
            }
          }
          if (degree>=4){
            for (j1 in 1:d) {
              for (j2 in j1:d){
                for (j3 in j2:d){
                  mphi0<-cbind(mphi0,phi0[,j1]*phi0[,j2]*phi0[,j3]*phi0[,j3:d])
                  mphi1<-cbind(mphi1,phi1[,j1]*phi1[,j2]*phi1[,j3]*phi1[,j3:d])
                }
              }
            }
          }else if (degree>4){
            warning('Five way interaction is not recommented.')
          }
        }
      }


    }

    svdphi<-svd(mphi0)
    newphi0=svdphi$u
    newphi1=mphi1 %*% svdphi$v %*% diag(1/svdphi$d)

    d=dim(newphi0)[2]
    for(i in 1:d){
      sum1=mean(newphi0[,i]*y0)*mean(newphi1[,i])
      sum2=mean(newphi0[,i])*mean(newphi1[,i])
      tsum1=tsum1+sum1
      tsum2=tsum2+sum2
    }

    if (tsum2==0) mcontrol=0
    else  mcontrol=tsum1/tsum2
    mtreat=mean(y1)

    Re=(mtreat-mcontrol)

    weight=apply(newphi0%*%t(newphi1),1,sum)
    return(list(est=Re,weight=weight/sum(weight),control_t=phi0,treated_t=phi1))
  }else{
    stop('method should be Kernel or Projection.')
  }
}
