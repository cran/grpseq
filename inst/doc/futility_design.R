## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=F-------------------------------------------------------------------
#  obj <- fut(alpha, beta, t, gamma, side = 2, si = 0, scale = "CP")

## ----eval=F-------------------------------------------------------------------
#  summary(obj)

## ----setup--------------------------------------------------------------------
## load the package
library(grpseq)
## two-sided level 0.05 test with 80% power;
## evenly spaced three futility looks with predictive power 20%;
## inflate sample size to recoup power.
obj1 <- fut(alpha=0.05,beta=0.2,t=(1:3)/4,gamma=0.2*rep(1,3),side=2,scale="PP",si=1)
obj1

## -----------------------------------------------------------------------------
## do the same thing without sample size inflation
obj2 <- fut(alpha=0.05,beta=0.2,t=(1:3)/4,gamma=0.2*rep(1,3),side=2,scale="PP",si=0)
obj2
## print the summary results
summary(obj2)

## ---- fig.height = 5, fig.width=7.2-------------------------------------------
oldpar <- par(mfrow = par("mfrow"))
par(mfrow=c(1,2))
## plot the futility boundaries by z-value
plot(obj2,scale='z',lwd=2,main="")
## plot the futility boundaries by B-value
plot(obj2,scale='b',lwd=2,main="")
par(oldpar)

## ---- fig.height = 5, fig.width=7.2-------------------------------------------
## plot the power curve as a function of the (local)
## effect size in units of the hypothesized effect size
## ref=TRUE requests the power curve for the original one-time analysis
powerplot(obj2,lwd=2, ref=TRUE)

