\name{multispeciesPP-package}
\alias{multispeciesPP-package}
\alias{multispeciesPP}
\docType{package}
\title{
  Fitting proportional bias models to presence-absence and presence-only data.
}
\description{
This package implements a proportional-bias model for presence-only and
presence-absence data from multiple species.
}
\details{
\tabular{ll}{
Package: \tab multispeciesPP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-28\cr
License: \tab What license is it under?\cr
}
}
\author{
William Fithian

Maintainer: William Fithian <wfithian@stanford.edu>
}
\references{
  "A proportional bias model for multispecies distribution modeling."
  William Fithian, Jane Elith, Trevor Hastie, and David A. Keith.
}


\keyword{ package }
\seealso{
}
\examples{
n.pixel <- 1000
n.other.spec <- 20
spec.names <- letters[1:(n.other.spec+1)]

## Geographic covariates affecting species abundance
x <- matrix(rnorm(2*n.pixel),nrow=n.pixel)

## Geographic covariate causing selection bias (correlated with x1)
z <- scale(x[,1] + rnorm(n.pixel)*sqrt(.95^(-2)-1))

## Intercept and slopes for abundance rate
alpha <- c(-2,.3*rnorm(n.other.spec)-2)
beta <- cbind(c(1,-.5),matrix(rnorm(n.other.spec*2)/2,2))

## Intercept and slope for selection bias
gamma <- -4
delta <- -.3

## PO data is impacted by selection bias
po.count <- matrix(rpois(n.pixel*(n.other.spec+1),lambda=exp(rep(alpha,each=n.pixel) + x \%*\% beta + gamma + c(z) * delta)),
                   n.pixel,dimnames=list(NULL,spec.names))

PO.list <- lapply(spec.names,function(sp) data.frame(x1=x[,1],x2=x[,2],z=z)[rep(1:n.pixel,po.count[,sp]),])
names(PO.list) <- spec.names
BG <- data.frame(x1=x[1,],x2=x[,2],z=z)

## PA data is unbiased
n.sites <- 500
pa.samp <- sample(1:n.pixel,n.sites)
pa.count <- matrix(rpois(n.sites*(n.other.spec+1),lambda=exp(rep(alpha,each=n.sites) + x[pa.samp,] \%*\% beta)),
                   n.sites,dimnames=list(NULL,spec.names))
table(pa.count>0)

## PA data doesn't need biasing covariates

PA <- cbind(data.frame(x1=x[pa.samp,1],x2=x[pa.samp,2]),as.data.frame(1*(pa.count>0)))

full.mod <- multispeciesPP(~x1+x2,~z,PA,PO.list,BG,region.size=n.pixel)

summary(full.mod)

plot(c(rbind(alpha,beta)),c(full.mod$species.coef[1:3,]),
     xlab="true coefficients",ylab="fitted coefficients")
abline(0,1)
}
