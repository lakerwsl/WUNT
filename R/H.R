H<-function(x,type='constant'){
  if (type=='constant'){
    if (x<=-0.5) return (0)
    if (x>0.5) return (1)
    if ((x<=0.5) & (x>-0.5)) return (x+0.5)
  }else if (type=='linear'){
    if (x<=-0.5) return (0)
    if (x>0.5) return (1)
    if ((x<=0.5) & (x>0)) return (-2*x*x+2*x+1/2)
    if ((x<=0) & (x>-0.5)) return (2*x*x+2*x+1/2)
  }else if (type=='quadratic'){
    a<-16
    b<-2
    c<-16
    if (x<=-0.5) return (0)
    if (x>0.5) return (1)
    if ((x<=-0.25) & (x>-0.5)) return ((x+0.5)^3*a/3)
    if ((x<=0.25) & (x>-0.25)) return (b*x-x^3*c/3+1/2)
    if ((x<=0.5) & (x>0.25)) return ((x-0.5)^3*a/3+1)
  }else{
    stop('type should be one of constant, linear and quadratic')
  }
}
