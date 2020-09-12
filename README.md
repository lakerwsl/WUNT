# WUNT
This package implements the new weighting framework by uniform transformer (WUNT), proposed by Yu and Wang (2020).

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("lakerwsl/WUNT")
```

Example usage:

```R
library(WUNT)

x1=rnorm(1000)
x2=rnorm(1000, mean=0.5)
x=cbind(x1,x2)
ix=sample(1:1000,300)
z=numeric(1000)
z[ix]=1
y=x1*x2
res1=wunt(x,z,y,type='constant',marginal=FALSE,method='Kernel')
res1$est
res2=wunt(x,z,y,type='constant',marginal=FALSE,method='Projection')
res2$est
```

#### References
Ruoqi Yu and Shulei Wang.
<b>Covariate Balancing by Uniform Transformer.</b>
2020.
[<a href="https://arxiv.org/pdf/2008.03738.pdf">arxiv</a>]