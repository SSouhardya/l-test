trunc_norm_cdf<-function(x,a,b,mu,sigma){
	#a<b is necessary
	if(a>b){
		stop('Lower truncation point larger than higher')
	}
	normalizer = pnorm(a,mean = mu, sd = sigma) + pnorm(b, mean = mu, sd = sigma, lower.tail = FALSE)
	if(x<=a){
		return(pnorm(x,mean = mu, sd = sigma)/normalizer)
	}
	else if(x>=b){
		return( (pnorm(x,mean = mu, sd = sigma) - 1 + normalizer)/normalizer )
	}
	else{
		return(pnorm(a,mean = mu, sd = sigma)/normalizer)
	}
}


liu_full_interval<-function(y,X,ind,lambda, coverage, interval, sigma = 1){
	#interval is the interval where uniroot will find the solution
	require(CVXR)

	alpha = 1-coverage
	n = nrow(X)
	p = ncol(X)

	lambda = n*lambda
	#the initial lambda was chosen for the LASSO objective with (1/2n) in the denominator of the squared error loss.

	beta = Variable(ncol(X))
	obj = sum((y-X%*%beta)^2)/(2) + lambda*p_norm(beta,1)
	prob = Problem(Minimize(obj))
	beta_coef = as.vector( round(solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9) )

	if(beta_coef[ind]==0){
		return(-1)
	}

	eta = X%*%solve(t(X)%*%X)[,ind]
	v = (diag(n) - eta%*%t(eta)/sum(eta^2))%*% y
	r = X[,-ind]%*%beta_coef[-ind] - v
	a = sum(eta^2)*(sum(X[,ind]*r)-lambda)
	b = sum(eta^2)*(sum(X[,ind]*r)+lambda)
	stat = sum(eta*y)

	if(stat<a | stat>b){
		fl<-function(mu){
			return(trunc_norm_cdf(stat,a,b,mu,sigma*sqrt(sum(eta^2)))-alpha/2)
		}
		fu<-function(mu){
			return(trunc_norm_cdf(stat,a,b,mu,sigma*sqrt(sum(eta^2)))-1+alpha/2)
		}
		L = uniroot(fl, interval = interval)$root
		U = uniroot(fu, interval = interval)$root
		return(range(c(L,U)))
	}
	else{
		warning('Statistic value inside the truncation interval but LASSO estimate is non-zero')
		return(NA) 
	}
}
