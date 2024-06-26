source('l_testing.R')

t_test_pval<-function(y,X,ind,sgn){
	y = std(as.vector(y))
	n = nrow(X)
	p = ncol(X)
	fit = lm(y~X)
	pval_both = summary(fit)$coefficients[1+ind,4]
	pval_one = pt(summary(fit)$coefficients[1+ind,3],df = n-p-1, lower.tail = (sgn<0) )
	pval_left = pt(-abs(summary(fit)$coefficients[1+ind,3]),df = n-p-1, lower.tail = TRUE )
	pval_right = pt(abs(summary(fit)$coefficients[1+ind,3]),df = n-p-1, lower.tail = FALSE )
	return(c(pval_both, pval_one, pval_left, pval_right))
}


hetero1<-function(X,par){
	n = nrow(X)
	vals = apply(X,1,function(x){sum(x^2)})
	cutoff = median(vals)
	variance = rep(1,n)
	variance[which(vals>cutoff)] = par
	error = vector(length = n)
	for(i in 1:n){
		error[i] = rnorm(1,sd = sqrt(variance[i]))
	}
	return(error)
}



#COMPARE ROBUSTNESS
compare_robustness<-function(n,p,s,a,b,rho,alpha,times,dist_type,par,set_seed = TRUE){
	require(MASS)
	require(glmnet)

	power = rep(0,3)
	pval_and_beta = matrix(0,nrow = times, ncol = 3)
	pval_and_beta = as.data.frame(pval_and_beta)
	names(pval_and_beta) = c('LASSO','One_sided','Two_sided')
	lambda_vals = vector(length = times)
	for(t in 1:times){
		if(set_seed){
			set.seed(t)
		}
		X = mvrnorm(n, mu = rep(0,p), Sigma = toeplitz(rho^(0:(p-1))))
		X = apply(X,2,g)
		if(dist_type == 'non_linear'){
			error = rnorm(n)
			Xproxy = X
			Xproxy = apply(Xproxy,2,function(x){return(sign(x)*(abs(x)^par))})	#par is the exponent
			amplitudes = (2*rbinom(s,1,0.5)-1)*c(a,rep(b,s-1))
			random_ind = sample(1:p, size = s, replace = FALSE)
			ind = random_ind[1]
			beta = rep(0,p)
			beta[random_ind] = amplitudes
			y = Xproxy%*%beta+error
		} else{
			if(dist_type == 't'){
				error = rt(n,df = par)
				if(par>2){
					error = error/(sqrt(par/(par-2)))
				}
			} else if(dist_type == 'gamma'){
				error = (rgamma(n, shape = par, rate = 1) - par)/sqrt(par)
			} else if(dist_type == 'hetero1'){
				error = hetero1(X,par)
			} else if(dist_type == 'hetero2'){
				error = hetero2(X,par)
			} else if(dist_type == 'gaussian'){
				error = rnorm(n, sd = sqrt(par))
			}
			amplitudes = (2*rbinom(s,1,0.5)-1)*c(a,rep(b,s-1))
			random_ind = sample(1:p, size = s, replace = FALSE)
			ind = random_ind[1]
			beta = rep(0,p)
			beta[random_ind] = amplitudes
			y = X%*%beta+error
		}
		Z = cbind(rep(1,n), X[,-ind])
 		proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		y.hat = proj%*%y
 		sigma.hat = sqrt(sum(((diag(n)-proj)%*%(y))^2))
 		V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
 		u = rnorm(n-ncol(Z))
   		u = u/sqrt(sum(u^2))
    	y_new = y.hat + sigma.hat*V%*%u
    	lambda = cv.glmnet(X[,-ind],y_new,standardize=FALSE)$lambda.min

		lambda_vals[t] = lambda
		pval_lasso = l.test(y,X,ind,lambda_cv = lambda, adjusted = FALSE, smoothed = TRUE)
		pval_t = t_test_pval(y,X,ind, sign(amplitudes[1]))
		pval_t_tails = pval_t[3:4]
		pval_t = pval_t[1:2]
		pval_and_beta[t,] = c(pval_lasso, pval_t)
		power = power + (c(pval_lasso, pval_t)<=alpha)
	}
	return(list(power/times, cbind(pval_and_beta, lambda_vals) ))
}


#--------------------------

par_number = 8
n = c(6,10,16,20,26,30,36,40)
p = n/2
s = rep(2,par_number)
a_vals = rep(0,par_number)
b_vals = rep(3.3, par_number)
rho = rep(0,par_number)

dist_type_vec = c(rep('t',6),rep('gamma',6), rep('hetero1',6), rep('non_linear',6))
par_vec = c(c(30,20,15,10,5,2),c(1,2,4,6,8,10), c(0.01,0.25,0.5,1,4,8),c(0.3,0.5,1,2,3,4))

times = 1000

run_multiple<-function(times,dist_type,par){
	D = matrix(0, nrow = length(a_vals), ncol = 4)
	L = list(length = length(a_vals))
	D[,1] = a_vals
	for(i in 1:length(a_vals)){
		L[[i]] = compare_robustness(n[i],p[i],s[i],a_vals[i],b_vals[i],rho[i],0.05,times,dist_type,par,set_seed = TRUE)
		D[i,-1] = L[[i]][[1]]
	}
	D = as.data.frame(D)
	names(D) = c('a','LASSO','t_both','t_one')
	saveRDS(list(D,L), paste0('robustness_ltest_size_',dist_type,par,'.list'))
	return(D)
}

for(i in 1:length(dist_type_vec)){
	print(paste('Doing:',dist_type_vec[i], par_vec[i]))
	D = run_multiple(times, dist_type_vec[i], par_vec[i])
}

