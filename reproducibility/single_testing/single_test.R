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

#COMPARE THE UNCONDITIONAL TESTS

compare_unconditional_tests<-function(n,p,s,a,b,rho,lambda,alpha,times, set_seed = TRUE){
	require(MASS)
	require(glmnet)
	lambda_s = lambda
	power = rep(0,3)
	pval_and_beta = matrix(0,nrow = times, ncol = 3)
	names(pval_and_beta) = c('LASSO','One_sided','Two_sided')
	pval_and_beta = as.data.frame(pval_and_beta)
	lambda_vals = vector(length = times)
	for(t in 1:times){
		if(set_seed){
			set.seed(t)
		}
		X = mvrnorm(n, mu = rep(0,p), Sigma = toeplitz(rho^(0:(p-1))))
		X = apply(X,2,g)
		amplitudes = (2*rbinom(s,1,0.5)-1)*c(a,rep(b,s-1))
		random_ind = sample(1:p, size = s, replace = FALSE)
		ind = random_ind[1]
		beta = rep(0,p)
		beta[random_ind] = amplitudes
		y = X%*%beta+rnorm(n)

		if(lambda_s == -1){	
			lambda = cv.glmnet(X[,-ind],y,standardize=FALSE)$lambda.min
		} else if(lambda_s == -2){ #this is the recommended choice of lambda
			 Z = cbind(rep(1,n), X[,-ind])
 		 	 proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		 	 y.hat = proj%*%y
 		 	 sigma.hat = sqrt(sum(((diag(n)-proj)%*%(y))^2))
 		 	 V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
 		 	 u = rnorm(n-ncol(Z))
   			 u = u/sqrt(sum(u^2))
    		 y_new = y.hat + sigma.hat*V%*%u
    		 lambda = cv.glmnet(X[,-ind],y_new,standardize=FALSE)$lambda.min
		} else if(lambda_s == -3){
			 Z = cbind(rep(1,n), X[,-ind])
 		 	 proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		 	 y.hat = proj%*%y
 		 	 sigma.hat = sqrt(sum(((diag(n)-proj)%*%(y))^2))
 		 	 V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
 		 	 u = rnorm(n-ncol(Z))
   			 u = u/sqrt(sum(u^2))
    		 y_new = y.hat + sigma.hat*V%*%u
    		 lambda = cv.glmnet( cbind(proj%*%X[,ind],X[,-ind]),y_new,standardize=FALSE)$lambda.min
		} else if(lambda_s == -4){
			 Z = cbind(rep(1,n), X[,-ind])
 		 	 proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		 	 y.hat = proj%*%y
 		 	 sigma.hat = sqrt(sum(((diag(n)-proj)%*%(y))^2))
 		 	 V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
 		 	 u = rnorm(n-ncol(Z))
   			 u = u/sqrt(sum(u^2))
    		 y_new = y.hat + sigma.hat*V%*%u
    		 lambda = cv.glmnet( X,y_new,standardize=FALSE)$lambda.min
		} else if(lambda_s == -5){
			 Z = cbind(rep(1,n), X[,-ind])
 		 	 proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		 	 y.hat = proj%*%y
    		 lambda = cv.glmnet( X[,-ind],y.hat,standardize=FALSE)$lambda.min
		} else if(lambda_s == -6){
			 Z = cbind(rep(1,n), X[,-ind])
 		 	 proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		 	 y.hat = proj%*%y
    		 lambda = cv.glmnet( cbind(proj%*%X[,ind],X[,-ind]),y.hat,standardize=FALSE)$lambda.min
		} else if(lambda_s == -7){
			 Z = cbind(rep(1,n), X[,-ind])
 		 	 proj = Z%*%solve(t(Z)%*%Z,t(Z))
 		 	 y.hat = proj%*%y
    		 lambda = cv.glmnet( X,y.hat,standardize=FALSE)$lambda.min
		} else if(lambda_s == - 8){
			lambda = cv.glmnet(X,y,standardize=FALSE)$lambda.min
		} else if(lambda_s == - 9){
			lambda = cv.glmnet(X[,-ind],y,standardize=FALSE)$lambda.min
		} else if(lambda_s>= 0){
			lambda = lambda_s
		} else{
			stop('Wrong choice')
		}
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

#----------------------------------

#Select one of the block of these parameter choices

#Varying amplitude with n and d far
#par_number = 7
#n = rep(100, par_number)
#p = rep(50, par_number)
#a_vals = c(0,1.3,2.3,3.3,4.3,5.3,6.3)
#s = rep(5,par_number)
#rho = rep(0,par_number)
#file_name = 'unconditional_test_summary_setting1'
#times = 200

#Varying amplitude with n and d close
#par_number = 7
#n = rep(100, par_number)
#p = rep(90, par_number)
#a_vals = c(0,1.3,2.3,3.3,4.3,5.3,6.3)
#s = rep(5,par_number)
#rho = rep(0,par_number)
#file_name = 'unconditional_test_summary_setting2'
#times = 200

#Varying inter-variable correlations
#par_number = 6
#n = rep(100, par_number)
#p = rep(50, par_number)
#a_vals = rep(4.3, par_number)
#s = rep(5,par_number)
#rho = c(0,0.1,0.3,0.5,0.7,0.9)
#file_name = 'unconditional_test_summary_correlation'
#times = 200

#Varying sparsity
#par_number = 8
#n = rep(100, par_number)
#p = rep(50, par_number)
#s = c(1,5,10,15,20,30,40,50)
#a_vals = rep(4.3, par_number)
#rho = rep(0,par_number)
#file_name = 'unconditional_test_summary_sparsity' #this is irrelevant if running on cluster.
#times = 400

cluster = FALSE #indicates whether running on a cluster

if(cluster){
	slurm_ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
	set.seed(slurm_ind)
}

run_multiple<-function(times,cluster = FALSE){
	D = matrix(0, nrow = length(a_vals), ncol = 4)
	L = list(length = length(a_vals))
	D[,1] = a_vals
	for(i in 1:length(a_vals)){
		print(i)
		L[[i]] = compare_unconditional_tests(n[i],p[i],s[i],a_vals[i],a_vals[i],rho[i],-2,0.05,times,set_seed = cluster)
		D[i,-1] = L[[i]][[1]]
	}
	D = as.data.frame(D)
	names(D) = c('a','LASSO','t_both', 't_one')
	if(cluster){
		saveRDS(list(D,L), paste0('single_testing_',slurm_ind,'.list'))
	} else{
		saveRDS(list(D,L), paste0(file_name,'.list'))
	}
	return(D)
}


D = run_multiple(times)

