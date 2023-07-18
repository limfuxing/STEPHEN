# functions for fitting uncorrelated mvnbinom
logl.nbinom <- function(p,x,wt) {
  #wt[is.na(wt) | is.infinite(wt) | wt==0] <- 1e-04
  nlogl <- -sum(wt*dnbinom(x,mu=exp(p[1]),size=exp(p[2]),log=TRUE),na.rm=T)
  #nlogl <- ifelse(is.infinite(nlogl),1e+20,nlogl)
  nlogl 
}

nbinom.wtest <- function(x,wt) {
 exp(optim(p=c(2,2),x=x,wt=wt,fn=logl.nbinom)$p)
}


mstep.mvnbinom <- function (x, wt) {
    k = ncol(wt)
    p = ncol(x)
    lambda = list(mu1=numeric(),mu2=numeric(),size1=numeric(),size2=numeric())
    for (i in 1:k) {
      out1 <- nbinom.wtest(x=x[,1],wt=wt[,i])
      lambda$mu1[i] = out1[1]
      lambda$size1[i] = out1[2]
      out2 <- nbinom.wtest(x=x[,2],wt=wt[,i])
      lambda$mu2[i] = out2[1]
      lambda$size2[i] = out2[2]
    }
lambda
}

rmvnbinom.hsmm <- function(j, model) {
  c(rnbinom(1,mu = model$parms.emission$mu1[j],size=model$parms.emission$size1[j]),rnbinom(1,mu = model$parms.emission$mu2[j],size=model$parms.emission$size2[j]))
}

dmvnbinom.hsmm <- function(x, j, model) {
  dnbinom2(x,mu1=model$parms.emission$mu1[j],mu2=model$parms.emission$mu2[j],size1=model$parms.emission$size1[j],size2=model$parms.emission$size2[j])
}

dnbinom2 <- function(x,mu1,mu2,size1,size2) {
  dnbinom(x[,1],mu=mu1,size=size1)*dnbinom(x[,2],mu=mu2,size=size2)
}

