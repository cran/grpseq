


# power function (when futility rules of enforced) of theta given a conditional
# power (designed for theta_A) futility design
# Two steps:
# 1. b value <- t, gamma, alpha/side, thetaA
# 2. power <- Z > (b, z_a) | (Z_{K-1} >b and Z_K < z_a), where Z \sim MVN(t*theta, t \min t')
power.fun=function(theta,alpha,t,gamma,thetaA,side=2){

	#error check: side=1 or 2
	if (!any(side==(c(1,2)))){
					stop("Side must be 1 or 2!")
					}

	k=length(t)+1
	za=qnorm(1-alpha/side)
	zg=qnorm(1-gamma)
	b=za-thetaA*(1-t)-zg*sqrt(1-t)

	#mean and variance of the B vector under theta
	mu=theta*c(t,1)
	S=matrix(0,k,k)
	tmpt=c(t,1)
	for (i in 1:k){
		for (j in 1:k){
		S[i,j]=min(tmpt[i],tmpt[j])
			}
		}
	if (side==2){
		Psi=(pmvnorm(lower=c(b,za),mean=mu,sigma=S)+
		pmvnorm(lower=c(b,-Inf),mean=mu,sigma=S)-
		pmvnorm(lower=c(b,-za),mean=mu,sigma=S))[1]
			}
	else{
	  Psi=pmvnorm(lower=c(b,za),mean=mu,sigma=S)[1]
	  }

	return(Psi)
}


# b values based on alpha, t, gamma and thetaA
# b=za-thetaA*(1-t)-zg*sqrt(1-t)
b.value=function(alpha,t,gamma,thetaA){
	# k=length(t)+1
	za=qnorm(1-alpha/2)
	zg=qnorm(1-gamma)
	b=za-thetaA*(1-t)-zg*sqrt(1-t)
	return(b)
}


# Power at theta=thetaA when futility rules are enforced (see power.fun)
# smaller than 1-pnorm(thetaA-z_a)
psi.tilde=function(theta,alpha,t,gamma,side=2){
	return(power.fun(theta,alpha,t,gamma,theta,side=side))
}



#' Design of non-binding futility analysis at multiple points
#'
#' @description Design of non-binding futility looks at multiple information times based
#' on conditional power (CP), predictive power (PP), or condition power
#' under current estimate (CPd) (Gallo, Mao, and Shih, 2014).
#'
#' @param alpha Type I error.
#' @param beta Type II error (1 - power).
#' @param t A numeric vector of information times in \eqn{(0, 1)} for futility looks.
#' @param gamma A numeric vector of probabilities (whose meaning depends on
#' \code{scale}) at information times \code{t}.
#' @param scale Character string specifying the scaled used: \code{"CP"}, conditional power;
#'  \code{"PP"}, predictive power; \code{"CPd"}: condition power
#' under current estimate.
#' @param side \code{1}- or \code{2}-sided test.
#' @param increment Error for the numerical solution of the sample size inflation factor.
#' @param si \code{0}: without sample size inflation;
#' \code{1}: with sample size inflation.
#' @param seed Seed number for the randomized evaluation of multivariate normal distribution.
#' @return An object of class \code{fut} with the following components.
#' \code{gamma1}: conditional power at information times \code{t} converted from
#' the supplied \code{gamma} and \code{scale};
#' \code{theta}: local alternative associated with the actual power when the
#' futility rules of enforced;
#' \code{IF}: sample size inflation factor if \code{si}=1;
#' \code{loss}: power loss if \code{si}=0.
#'
#' @seealso \code{\link{print.fut}}, \code{\link{summary.fut}}, \code{\link{plot.fut}},
#' \code{\link{powerplot}}
#' @import mvtnorm
#' @importFrom stats pnorm
#' @export
#' @keywords fut
#' @references Gallo, P.,  Mao, L., and Shih, V.H. (2014).
#' Alternative views on setting clinical trial futility criteria.
#' Journal of Biopharmaceutical Statistics, 24, 976-993.
#' @examples
#' ## load the package
#' library(grpseq)
#' ## two-sided level 0.05 test with 80% power;
#' ## evenly spaced three futility looks with predictive power 20%;
#' ## inflate sample size to recoup power.
#' obj1 <- fut(alpha=0.05,beta=0.2,t=(1:3)/4,gamma=0.2*rep(1,3),side=2,scale="PP",si=1)
#' obj1
#' ## print the summary results
#' summary(obj1)
#'
#' ## do the same thing without sample size inflation
#' obj2 <- fut(alpha=0.05,beta=0.2,t=(1:3)/4,gamma=0.2*rep(1,3),side=2,scale="PP",si=0)
#' obj2
#' ## print the summary results
#' summary(obj2)
#' oldpar <- par(mfrow = par("mfrow"))
#' par(mfrow=c(1,2))
#' ## plot the futility boundaries by z-value
#' plot(obj2,scale='z',lwd=2,main="")
#' ## plot the futility boundaries by B-value
#' plot(obj2,scale='b',lwd=2,main="")
#' par(oldpar)
#' ## plot the power curve as a function of the (local)
#' ## effect size in units of the hypothesized effect size
#' ## ref=TRUE requests the power curve for the original one-time analysis
#' powerplot(obj2,lwd=2, ref=TRUE)

fut=function(alpha,beta,t,gamma,side=2,increment=1e-4,si=0,scale='CP',seed=12345){

  set.seed(seed)
  za <- qnorm(1-alpha/side)
	zb <- qnorm(1-beta)

	theta <- za+zb
	#add error checks here
		#error check: side=1 or 2
	if (!any(side==(c(1,2)))){
					stop("Side must be 1 or 2!")
					}

	# number of analyses (including final)
	k <- length(t)+1
	#initial value of theta
	#lower bound for theta
	a <- qnorm(1-alpha/side)+qnorm(1-beta)

	# convert gamma in "scale" to gamma1 in CP

	if (scale=='PP'){
		gamma1=pp2cp(alpha=alpha,beta=beta,t=t,pp=gamma,side=side)
			}else{
			if (scale=='CPd'){gamma1=cpd2cp(alpha=alpha,beta=beta,t=t,cpd=gamma,side=side)}
			else{gamma1=gamma}
			}

	loss=0
	IF <- 1
	if (si==0){
		theta=a
		loss=power.loss(alpha=alpha,beta=beta,t=t,gamma=gamma1,side=side)$loss
		}
	else{
	#find upper bound for theta
	c=za+4
	Psi=psi.tilde(theta,alpha,t,gamma1,side=side)

	while ( c-a>increment){
		theta=(a+c)/2
		Psi=psi.tilde(theta,alpha,t,gamma1,side=side)
		if (Psi<=1-beta){
			a=theta
			}else{c=theta}
		#print(theta)
	}
	# sample size inflation factor
	IF <- (theta/(qnorm(1-alpha/side)+qnorm(1-beta)))^2
	}



	obj <- list(alpha=alpha,beta=beta,t=t,gamma=gamma,scale=scale,gamma1=gamma1,
		theta=theta,side=side,si=si,loss=loss, IF=IF,call=match.call())

	class(obj) <- "fut"

	return(obj)
}


expectss=function(delta,object){
theta=object$theta
alpha=object$alpha
beta=object$beta
t=object$t
gamma1=object$gamma1
side=object$side
za=qnorm(1-alpha*(3-side)/2)
N=(theta/(za+qnorm(1-beta)))^2
b=b.value(alpha*(3-side),t,gamma1,theta)
z=b/sqrt(t)
k=length(t)+1
mu=delta*theta*c(t,1)
S=matrix(0,k,k)
tmpt=c(t,1)
for (i in 1:k){
		for (j in 1:k){
		S[i,j]=min(tmpt[i],tmpt[j])
			}
		}
stop.prob=rep(NA,k-1)
for (i in 1:(k-1)){
	if (i==1){
		stop.prob[1]=pmvnorm(upper=b[1],mean=mu[1],sigma=S[1,1])
		}else{
		stop.prob[i]=pmvnorm(lower=c(b[1:(i-1)],-Inf),mean=
				mu[1:i],sigma=S[1:i,1:i])-
				pmvnorm(lower=c(b[1:(i-1)],b[i]),mean=
				mu[1:i],sigma=S[1:i,1:i])
		}
	}
stop.prob[k]=1-sum(stop.prob)
t=c(t,1)
ess=sum(stop.prob*t)
return(list(ss=t*N,stop.prob=stop.prob,ess=ess*N))
}

# Power loss if futility rules are enforced with no sample size inflation
# A refined version of power.fun()
# two steps:
# 1. b value <- t, gamma, alpha/side, thetaA
# 2. power loss at t_k: (Z_1,...,Z_{k-1}, Z_K)>(b_1,...,b_{k-1}, za)
# & z_k<= b_k (ignore Z_K<-za)
power.loss=function(alpha,beta,t,gamma,side=2){

	k=length(t)+1
	za=qnorm(1-alpha/side)
	zb=qnorm(1-beta)
	# zg=rep(NA,k-1)
	# for (i in 1:(k-1)){
	# zg[i]=qnorm(1-gamma[i])
	# }
	zg=qnorm(1-gamma)
	theta=za+zb
	b=c(za-theta*(1-t)-zg*sqrt(1-t),za)
  # b[k]=za
  mu=theta*c(t,1)
  S=matrix(0,k,k)
  tmpt=c(t,1)
  for (i in 1:k){
		for (j in 1:k){
		S[i,j]=min(tmpt[i],tmpt[j])
			}
		}
  loss=rep(NA,k-1)
  for (i in 1:(k-1)){
	  if (i==1){
		  loss[1]=pmvnorm(lower=b[k],mean=mu[k],sigma=S[k,k])[1]-
			pmvnorm(lower=c(b[c(1,k)]),mean=mu[c(1,k)],sigma=S[c(1,k),c(1,k)])[1]
		}else{
		loss[i]=pmvnorm(lower=c(b[1:(i-1)],-Inf,b[k]),mean=
				mu[c(1:i,k)],sigma=S[c(1:i,k),c(1:i,k)])[1]-
				pmvnorm(lower=c(b[1:(i-1)],b[i],b[k]),mean=
				mu[c(1:i,k)],sigma=S[c(1:i,k),c(1:i,k)])[1]
		}
	}
b=b[-k]
return(list(theta=theta,b=b,loss=loss))
}


# Evenly distributed power loss -> conditional power and b (z) values
# step 0: get power loss l=tloss/length(t) at each look
# step 1: b1 <- root_b{power-loss(t1, b1)=l}
# ...
# step k: bk <- root_b{power-loss(t1,...t_k,b1,...bk)=l}
# conditional power <- b
# z <- b/sqrt(t)
evenly.powerloss=function(alpha,beta,tloss,t,side=2){

	za=qnorm(1-alpha*(3-side)/2)
	zb=qnorm(1-beta)
	theta=za+zb
	k=length(t)
	b=rep(NA,k+1)
	l=tloss/k
	b[k+1]=za
	mu=theta*c(t,1)
	S=matrix(0,k+1,k+1)
	tmpt=c(t,1)
	for (i in 1:(k+1)){
		for (j in 1:(k+1)){
		S[i,j]=min(tmpt[i],tmpt[j])
			}
		}



for (i in 1:k){
	if (i==1){
	#define the function to solve for b1
	f=function(x){
		result=pmvnorm(lower=b[k+1],mean=mu[k+1],sigma=S[k+1,k+1])-
			pmvnorm(lower=c(x,b[k+1]),mean=mu[c(1,k+1)],
			sigma=S[c(1,k+1),c(1,k+1)])-l
		return(result)
			}
	b[1]=bisection(f,-4*sqrt(t[1]),4*sqrt(t[1]))
		}


	else{
	f=function(x){

		result=pmvnorm(lower=c(b[1:(i-1)],-Inf,b[k+1]),mean=
				mu[c(1:i,k+1)],sigma=S[c(1:i,k+1),c(1:i,k+1)])-
				pmvnorm(lower=c(b[1:(i-1)],x,b[k+1]),mean=
				mu[c(1:i,k+1)],sigma=S[c(1:i,k+1),c(1:i,k+1)])-l

		return(result)
			}
	b[i]=bisection(f,-4*sqrt(t[i]),4*sqrt(t[i]))
			}
		}
	z=b/sqrt(c(t,1))
	gamma=1-pnorm((za-b[1:k]-theta*(1-t))/sqrt(1-t))
return(list(b=b,z=z,gamma=gamma))
	}


#use bisection method to solve f(x)=0
# with f monotone increasing function
bisection=function(f,a,b,err=0.0001){

while(b-a>err){
	x=mean(c(a,b))
	if (f(x)<0){
	a=x
	}else{b=x}
	}
return(x)

}


###program predictive power to conditional power
pp2cp=function(alpha,beta,t,pp,side=2){

	za=qnorm(1-alpha*(3-side)/2)
	zb=qnorm(1-beta)
	theta=za+zb
	k=length(t)
	zp=qnorm(1-pp)
#compute boundary
	b=t*za-sqrt(t*(1-t))*zp
	gamma=1-pnorm((za-b-theta*(1-t))/sqrt(1-t))

return(gamma)
}


###program predictive power to conditional power
cpd2cp=function(alpha,beta,t,cpd,side=2){

	za=qnorm(1-alpha*(3-side)/2)
	zb=qnorm(1-beta)
	theta=za+zb
	k=length(t)
	zp=qnorm(1-cpd)
#compute boundary
	b=t*(za-sqrt(1-t)*zp)
	gamma=1-pnorm((za-b-theta*(1-t))/sqrt(1-t))

return(gamma)
}



