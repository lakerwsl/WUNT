Kernel <- function(x,h,kernel='gaussian'){
  y=x/h
  d=ncol(x)

  if (kernel=='gaussian') Re=1/sqrt(2*pi)*exp(-y^2/2)
  else if (kernel=='gaussian4') Re=1/sqrt(2*pi)*exp(-y^2/2)*(3/2-y^2/2)
  else if (kernel=='gaussian6') Re=1/sqrt(2*pi)*exp(-y^2/2)*(15/8-5*y^2/4+y^4/8)
  else if (kernel=='gaussian8') Re=1/sqrt(2*pi)*exp(-y^2/2)*(35/16-35*y^2/16+7*y^4/16-y^6/48)
  else if (kernel=='gaussian10') Re=1/sqrt(2*pi)*exp(-y^2/2)*(315/128-105*y^2/32+63*y^4/64-3*y^6/32+y^8/384)
  else if (kernel=='epanechnikov') Re=3/4*(1-y^2)*(abs(y)<=1)
  else if (kernel=='uniform') Re=1/2*(abs(y)<=1)
  else if (kernel=='triangular') Re=(1-abs(y))*(abs(y)<=1)
  else if (kernel=='triweight') Re=35/32*(1-y^2)^3*(abs(y)<=1)
  else if (kernel=='tricube') Re=70/81*(1-abs(y)^3)^3*(abs(y)<=1)
  else if (kernel=='biweight') Re=15/16*(1-y^2)^2*(abs(y)<=1)
  else if (kernel=='cosine') Re=pi/4*cos(pi/2*y)*(abs(y)<=1)
  else if (kernel=='silverman') Re=1/2*exp(-abs(y)/sqrt(2))*sin(abs(y)/sqrt(2)+pi/4)

  if (d>1) Re=apply(Re,1,prod)

  return(Re)
}
