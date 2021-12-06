


#' Print basic information about the futility design
#' @description Print the power loss or sample size inflation factor due to
#' the planned futility analysis.
#' @param x An object of class \code{fut}.
#' @param ... Further arguments passed to or from other methods.
#' @return Print the results of \code{fut} object.
#' @seealso \code{\link{fut}}, \code{\link{summary.fut}}
#' @export
#' @keywords fut
#' @examples
#'# see example for fut
print.fut=function(x,...){
  theta=x$theta
  alpha=x$alpha
  beta=x$beta
  side=x$side
  si=x$si
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  if (si==0){
  loss=x$loss
  cat("Power loss due to futility rules: ")
  cat(sum(loss))
  cat(" (planned power: ")
  cat(1-beta)
  cat(").")
  cat("\n")
  }
  else{

  cat("Sample Size inflation factor: ")
  cat(x$IF)
  cat("\n")
  }
}



#' Detailed summary of the futility design
#'
#' Provide key information about the futility design, including B-/z-values, beta (type II error)
#' spent, and power loss at each futility look as well the the sample size distribution under
#' the null hypothesis.
#'
#' @param object An object returned by \code{\link{fut}}.
#' @param ... further arguments passed to or from other methods.
#' @return An object of class \code{summary.fut} with components:
#' \item{t}{A \eqn{K}-dimensional vector of information times.}
#' \item{b}{A \eqn{K}-dimensional vector of B-values at \code{t}.}
#' \item{z}{A \eqn{K}-dimensional vector of z-values at \code{t}.}
#' \item{type2}{A \eqn{K}-dimensional vector of beta spent at \code{t}.}
#' \item{loss}{A \eqn{K}-dimensional vector of power loss at \code{t}.}
#' \item{ess}{Expected sample size at \eqn{H_0}.}
#' \item{...}{}
#' @seealso \code{\link{fut}}, \code{\link{print.fut}}, \code{\link{print.summary.fut}}.
#' @keywords fut
#' @importFrom stats pchisq
#' @examples
#'# see example for fut
#' @export
summary.fut=function(object,...){
theta=object$theta
alpha=object$alpha
beta=object$beta
t=object$t
gamma=object$gamma
gamma1=object$gamma1
side=object$side
si=object$si
loss=object$loss
scale=object$scale

za=qnorm(1-alpha/side)
b=b.value(alpha*(3-side),t,gamma1,theta)
z=b/sqrt(t)
k=length(t)+1
mu=theta*c(t,1)
S=matrix(0,k,k)
tmpt=c(t,1)
for (i in 1:k){
		for (j in 1:k){
		S[i,j]=min(tmpt[i],tmpt[j])
			}
		}
type2=rep(NA,k)

for (i in 1:(k-1)){
	if (i==1){
		type2[1]=pmvnorm(upper=b[1],mean=mu[1],sigma=S[1,1])
		}else{
		type2[i]=pmvnorm(lower=c(b[1:(i-1)],-Inf),mean=
				mu[1:i],sigma=S[1:i,1:i])-
				pmvnorm(lower=c(b[1:(i-1)],b[i]),mean=
				mu[1:i],sigma=S[1:i,1:i])
		}
}
if (side==2){
type2[k]=pmvnorm(lower=c(b,-za),mean=mu,sigma=S)-
				pmvnorm(lower=c(b,za),mean=mu,sigma=S)
		}
	else{type2[k]=pmvnorm(lower=c(b,-Inf),mean=mu,sigma=S)-
				pmvnorm(lower=c(b,za),mean=mu,sigma=S)}

Ess0=expectss(0,object)


result=list(alpha=alpha,beta=beta,t=c(t,1),gamma=c(gamma,NA),
		theta=theta,side=side,call=object$call,b=c(b,za),z=c(z,za),
		type2=type2,si=si,loss=c(loss,NA),Ess0=Ess0,ess=Ess0$ess,scale=scale)

class(result)="summary.fut"
return(result)
}

#' Print method for summary.fut objects
#'
#' Print the detailed summary of the futility design.
#'
#' @param x An object returned by \code{\link{summary.fut}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @seealso \code{\link{fut}}, \code{\link{summary.fut}}.
#' @export
#' @importFrom stats qnorm
print.summary.fut=function(x,...){
cat("\n")
cat("Call:\n")
print(x$call)
cat("\n")
theta=x$theta
alpha=x$alpha
beta=x$beta
side=x$side
si=x$si
loss=x$loss
scale=x$scale

N=(theta/(qnorm(1-alpha*(3-side)/2)+qnorm(1-beta)))^2
cat("Sample Size inflation factor: ")
cat(N)
cat("\n")

k=length(x$t)

if (si==0){
output1=rbind(x$t,x$gamma,x$b,x$z,x$type2,x$loss)
colnames(output1)=rep(NA,k)
for (i in 1:(k-1)){colnames(output1)[i]=paste("Interim", i)}
colnames(output1)[k]="Final"


rownames(output1)=c("Info Time",scale,"B-Value","Z-value","Beta spent","Power Loss")


cat("Nominal type I error = ",x$alpha, ", Power loss = ", round(sum(loss[-k]),4),
" (Planned power: ",1-beta,"),\n", "Side = ",side, "\n\n",sep="")

print(round(output1,4))
cat("\n\n")
}else{
output1=rbind(x$t,x$gamma,x$b,x$z,x$type2)
colnames(output1)=rep(NA,k)
for (i in 1:(k-1)){colnames(output1)[i]=paste("Interim", i)}
colnames(output1)[k]="Final"
rownames(output1)=c("Info Time",scale,"B-Value","Z-value","Beta spent")

cat("Nominal type I error=: ",x$alpha, ", Overall power=", 1-x$beta,", Side=",
	side, "\n\n")
print(round(output1,4))
cat("\n\n")

}

#expected sample size under H0
output2=rbind((x$Ess0)$ss,(x$Ess0)$stop.prob)
colnames(output2)=colnames(output1)
rownames(output2)=c("Sample Size", "Probability")
cat("Sample size distribution under H_0:\n")
cat("(as a ratio to that of the reference test)\n")
cat("Expected sample size:",
	round((x$Ess0)$ess,4),"\n\n")
print(round(output2,4))
cat("\n")
}






#' Plot the planned futility boundaries
#'
#' Plot the planned futility boundaries in B- or z-values as a function
#' of information time.
#'
#' @param x An object returned by \code{\link{fut}}.
#' @param scale \code{"z"}: plot z-values; \code{"b"}: plot B-values.
#' @param add If TRUE, the curve will be overlaid on an existing plot; otherwise,
#' a separate plot will be constructed.
#' @param main A main title for the plot.
#' @param xlim The x limits of the plot.
#' @param ylim The y limits of the plot.
#' @param xlab A label for the x axis, defaults to a description of x.
#' @param ylab A label for the y axis, defaults to a description of y.
#' @param lty Line type for the segments connecting the z-/B-value points.
#' @param pch Point types for the z-/B-values.
#' @param type Plot type. \code{"l"}: only line segments; \code{"p"}: only z-/B-value points;
#' \code{"b"}: both.
#' @param cex Point size.
#' @param ... Other arguments that can be passed to the underlying \code{plot} method.
#' @return No return value, called for side effects.
#' @seealso \code{\link{fut}}, \code{\link{summary.fut}}, \code{\link{powerplot}}.
#' @keywords fut
#' @importFrom stats qnorm
#' @importFrom graphics lines points abline
#' @export
#' @examples
#'# see example for fut
plot.fut=function(x,scale="z", add=FALSE,lty=8,xlab="Info Time",
ylab="z score", type="b",pch=1,cex=1,
main="Futility Boundary for the Planned Test",xlim=c(0,1.1),
 ylim=NULL,...){
theta=x$theta
alpha=x$alpha
beta=x$beta
t=x$t
gamma=x$gamma
side=x$side

za=qnorm(1-alpha/side)
b=b.value(alpha*(3-side),t,gamma,theta)
z=b/sqrt(t)
k=length(t)+1
if (scale=="b"){
	z=b
	ylab="B score"
		}else{
		if (scale!="z"){stop("scale must be 'z' or 'b'")}
			}
z=c(z,za)
t=c(t,1)
if (is.null(ylim)){ylim=1.2*c(-za,za)}
if (add==FALSE){
	plot(t,z,main=main,xlim=xlim,ylim=ylim,
		lty=lty,pch=pch,xlab=xlab,ylab=ylab,
		type=type,cex=cex,...)
		points(1,za,pch=20,cex=cex)
		abline(h=0)
		}else{
	lines(t,z,main=main,xlim=xlim,ylim=ylim,
		lty=lty,pch=pch,xlab=xlab,ylab=ylab,
		type=type,cex=cex,...)
		points(1,za,pch=20,cex=cex)
		abline(h=0)
	}
}

#' Plot the power function of the planned analysis
#'
#' Plot the power curve of the planned futility analysis as a function of
#' the effect size (in units of the hypothesized effect size).
#'
#' @param x An object returned by \code{\link{fut}}.
#' @param ref If TRUE, power curve of the reference test (one that ignores the futility boundaries)
#' will be overlaid.
#' @param add If TRUE, the curve will be overlaid on an existing plot; otherwise,
#' a separate plot will be constructed.
#' @param lty Line type for the power curve of the futility analysis.
#' @param ref.lty Line type for the power curve of the reference if \code{ref=TRUE}.
#' @param lwd Line width.
#' @param main A main title for the plot.
#' @param xlim The x limits of the plot.
#' @param ylim The y limits of the plot.
#' @param xlab A label for the x axis, defaults to a description of x.
#' @param ylab A label for the y axis, defaults to a description of y.
#' @param ... Other arguments that can be passed to the underlying \code{plot} method.
#' @return No return value, called for side effects.
#' @seealso \code{\link{fut}}, \code{\link{summary.fut}}, \code{\link{plot.fut}}.
#' @keywords fut
#' @importFrom graphics legend
#' @export
#' @examples
#'# see example for fut
powerplot=function(x,ref=FALSE,add=FALSE,lty=1,ref.lty=2,lwd=1,xlab=expression(delta),
ylab="Power", main="Power curve of the planned futility analysis",xlim=c(0,1.5),
 ylim=c(0,1),...){
theta=x$theta
IF=x$IF
alpha=x$alpha
beta=x$beta
t=x$t
gamma=x$gamma1
side=x$side
za=qnorm(1-alpha/side)
b=b.value(alpha*(3-side),t,gamma,theta)
z=b/sqrt(t)
k=length(t)+1

delta=seq(0,1.5,by=0.05)
power=rep(NA,length(delta))
for (i in 1:length(delta)){
	power[i]=power.fun(delta[i]*theta,alpha,t,gamma,theta,side=side)
	}
if (!add){
plot(delta,power,lty=lty,xlab=xlab,
ylab=ylab, type='l',
main=main,xlim=xlim,lwd=lwd,
 ylim=ylim,...)
	}else{
		lines(delta,power,lty=lty,xlab=xlab,
		ylab=ylab,
		main=main,xlim=xlim, lwd=lwd,
 		ylim=ylim,...)
		}
if (ref){
	for (i in 1:length(delta)){
	power[i]=1-pnorm(za-delta[i]*(za+qnorm(1-beta)))+ (side==2)*
		pnorm(-za-delta[i]*(za+qnorm(1-beta)))
		}
	lines(delta,power,lty=ref.lty,lwd=lwd)
	legend(0,1,lty=c(lty,ref.lty),lwd=lwd,c("Designed","Reference"))
	}
}



