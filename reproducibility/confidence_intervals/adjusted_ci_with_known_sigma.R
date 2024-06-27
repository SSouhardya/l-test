source('~/utilities.R')

beta_x<-function(x,y,X,ind,lamb, solver = NA){
	require(CVXR)
	n = nrow(X)
	Z = cbind(rep(1,n),X[,-ind])
	beta = Variable(ncol(X))
	obj = sum((y-x*X[,ind]-Z%*%beta)^2)/(2*n) + lamb*p_norm(beta[-1],1)
	prob0 = Problem(Minimize(obj))
	beta.x = round(solve(prob0,feastol=1e-10,reltol=1e-10,abstol=1e-10, solver = solver )$getValue(beta), digits = 9)
	return(beta.x)
}

beta_full<-function(y,X,ind,lamb,solver = NA){
	require(CVXR)
	n = nrow(X)
	Z = cbind(rep(1,n),X)
	beta = Variable(ncol(X)+1)
	obj = sum((y-Z%*%beta)^2)/(2*n) + lamb*p_norm(beta[-1],1)
	prob0 = Problem(Minimize(obj))
	beta.full = round(solve(prob0,feastol=1e-10,reltol=1e-10,abstol=1e-10, solver = solver )$getValue(beta), digits = 9)
	return(beta.full)
}


l.cdf_adjusted<-function(x,y,X,ind,gamma, lambda_cv , lambda , display = FALSE, solver = NA, tail = 'left', smoothed = FALSE, known_sigma = FALSE){
	require(CVXR)
	n = nrow(X)
	p = ncol(X)

	Z = cbind(rep(1,n),X[,-ind])
	proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
	y.hat = proj %*%y
	sigma.hat = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
	second_denom_term = sqrt(sum(((diag(n) - proj)%*%X[,ind])^2))

    sgn<-function(x){
        if(x == 0){
            return(1)
        }
        return(sign(x))
    }

    if(smoothed & (x==0)){
    	beta.x = beta_x(0,y - gamma*X[,ind],X,ind,lambda_cv, solver = solver)
    	mid = -sum(X[,ind]*(y.hat - gamma*(proj)%*%X[,ind] - Z%*%beta.x))/(sigma.hat * second_denom_term)
    	u1 = sum(X[,ind]*(y - y.hat - gamma*(diag(n) - proj)%*%X[,ind] ))/(sigma.hat*second_denom_term)
    	dist_from_mid = abs(mid - u1)
    	if(tail == 'left'){
    		v.x = mid - dist_from_mid
    	} else{
    		v.x = mid + dist_from_mid
    	}
    } else{
		beta.x = beta_x(x,y - gamma*X[,ind],X,ind,lambda_cv, solver = solver)
		v.x = ( -sum(X[,ind]*(y.hat - gamma*(proj)%*%X[,ind]- x*X[,ind] - Z%*%beta.x)) + n*lambda_cv*sgn(x) )/(sigma.hat*second_denom_term)
	}
	beta.0 = beta_x(0,y,X,ind,lambda, solver = solver)

	v1 = ( -sum(X[,ind]*(y.hat + gamma*(diag(n) - proj)%*%X[,ind] - Z%*%beta.0)) - n*lambda )/(sigma.hat*second_denom_term)
    v2 = ( -sum(X[,ind]*(y.hat + gamma*(diag(n) - proj)%*%X[,ind] - Z%*%beta.0)) + n*lambda )/(sigma.hat*second_denom_term)

    denom = 1 - (qhaar(v2,n-ncol(X), lower.tail = TRUE, known_sigma = known_sigma, sigma.hat = sigma.hat) - qhaar(v1,n-ncol(X), lower.tail = TRUE, known_sigma = known_sigma, sigma.hat = sigma.hat) )	# i am calulating this outside as this give sthe dropping probability
    numer_1 = qhaar(v.x, n-ncol(X), lower.tail = TRUE, known_sigma = known_sigma, sigma.hat = sigma.hat)
    numer_2 = 0
    if(v.x>v1){
    	numer_2 = qhaar(min(c(v.x,v2)), n- ncol(X), lower.tail = TRUE, known_sigma = known_sigma, sigma.hat = sigma.hat) - qhaar(v1, n- ncol(X), lower.tail = TRUE, known_sigma = known_sigma, sigma.hat = sigma.hat)
    }
    cond_prob = (numer_1 - numer_2)/denom

    if(tail == 'right'){
    	cond_prob = 1-cond_prob
    }
    return(cond_prob)

    toc = Sys.time()
    if(display){
    	print(toc - tic)
    }
}

l.ci_adjusted<-function(y,X,ind, gamma_range, lambda,  lambda_cv=-1, coverage = 0.95, display = FALSE, smoothed = FALSE, outer_approx = FALSE, outer_grid.length = 10, known_sigma = FALSE){
	require(CVXR)

	g.length = length(gamma_range)
	lambda_cv_flag = lambda_cv

	cond_pval = vector(length = length(gamma_range))
	right.pval = cond_pval
	left.pval = cond_pval
	lambda_vec = left.pval
	lambda_cv_vec = lambda_vec

	n = nrow(X)
	p = ncol(X)

	Z = cbind(rep(1,n),X[,-ind])
	proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
	y.hat = proj %*%y
	V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
	second_denom_term = sqrt(sum(((diag(n) - proj)%*%X[,ind])^2))

	beta = Variable(ncol(X)+1)
	obj = sum(( y - cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda*p_norm(beta[-1],1)
	prob = Problem(Minimize(obj))
	beta.hat = round(solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10,)$getValue(beta), digits = 9)

	if(beta.hat[1+ind] == 0){
		return(-1)	#no confidence interval if there is no selection
	}

	for(i in 1:length(gamma_range)){
		tic1 = Sys.time()
		gamma = gamma_range[i]
		if(lambda_cv_flag <0){
			sigma.hat.gamma = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
			u = rnorm(n-ncol(Z))
			if(known_sigma){
				u = u/sigma.hat.gamma
			} else{
				u = u/sqrt(sum(u^2))
			}
			y_temp = as.vector(y.hat + gamma*(diag(n)-proj)%*%X[,ind] + sigma.hat.gamma*V%*%u)
			glmnet_object_cv = cv.glmnet(X[,-ind], y_temp-gamma*X[,ind], standardize=FALSE)
			if(lambda_cv_flag == -1){
				lambda_cv = glmnet_object_cv$lambda.min
			} else if(lambda_cv_flag == -2){
				lambda_cv = glmnet_object_cv$lambda.1se
			}
		}
		lambda_vec[i] = lambda
		lambda_cv_vec[i] = lambda_cv

		beta = Variable(ncol(X)+1)
		obj = sum(( (y - gamma*X[,ind]) -cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda_cv*p_norm(beta[-1],1)
		prob = Problem(Minimize(obj))
		beta.hat.stat = round(solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9)



		stat = abs(beta.hat.stat[1+ind])
		#print(stat)
		if(beta.hat[1+ind] == 0){
			cond_pval[i] = -1
		} else{
			right.pval[i] = l.cdf_adjusted(stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'right',known_sigma = known_sigma)
			left.pval[i] = l.cdf_adjusted(-stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'left',known_sigma = known_sigma)
			cond_pval[i] = right.pval[i] + left.pval[i]
		}
	}

	inds = which(cond_pval>1-coverage)
	if(length(inds)==0){
		warning('None of the grid elements selected. Try using a different grid?')
		return(vector(length = 0))
	}
	adjusted = TRUE
	if(adjusted & (lambda<0)){
		stop('If adjusted, the selection lambda must be supplied')
	}
	if(outer_approx){
		if(min(inds)!=1){
			inds1 = min(inds)-1
			left_additional = seq(from = gamma_range[inds1], to = gamma_range[min(inds)], length.out = outer_grid.length)

			pvals.left = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				gamma = left_additional[j]
				if(lambda_cv_flag <0){
					sigma.hat.gamma = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
					u = rnorm(n-ncol(Z))
					if(known_sigma){
						u = u/sigma.hat.gamma
					} else{
						u = u/sqrt(sum(u^2))
					}
					y_temp = as.vector(y.hat + gamma*(diag(n)-proj)%*%X[,ind] + sigma.hat.gamma*V%*%u)
					glmnet_object_cv = cv.glmnet(X[,-ind], y_temp-gamma*X[,ind], standardize=FALSE)
					if(lambda_cv_flag == -1){
						lambda_cv = glmnet_object_cv$lambda.min
					} else if(lambda_cv_flag == -2){
						lambda_cv = glmnet_object_cv$lambda.1se
					}
				}			
				beta = Variable(ncol(X)+1)
				obj = sum(( (y - gamma*X[,ind]) -cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda_cv*p_norm(beta[-1],1)
				prob = Problem(Minimize(obj))
				beta.hat.stat = round(solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9)
				stat = abs(beta.hat.stat[1+ind])

				right_tail = l.cdf_adjusted(stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'right', known_sigma = known_sigma)
				left_tail = l.cdf_adjusted(-stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'left', known_sigma = known_sigma)
				pvals.left[j] = left_tail + right_tail
			}
			inds.left = which(pvals.left > 1-coverage)
			inds.left = unique(c(inds.left, length(left_additional)))
			if(min(inds.left)!=1){
				inds.left = c(min(inds.left)-1, inds.left)
			}
			ci.left = left_additional[inds.left[1]]
		} else{
			ci.left = gamma_range[inds[1]]
			warning('CI hit the lower limit. Try increasing the range?')
		}
		if(max(inds)!=g.length){
			inds2 = max(inds)+1
			right_additional = seq(from = gamma_range[max(inds)], to = gamma_range[inds2], length.out = outer_grid.length)
			pvals.right = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				gamma = right_additional[j]
				if(lambda_cv_flag <0){
					sigma.hat.gamma = sqrt(sum(((diag(n) - proj)%*%(y - X[,ind]*gamma))^2))
					u = rnorm(n-ncol(Z))
					if(known_sigma){
						u = u/sigma.hat.gamma
					} else{
						u = u/sqrt(sum(u^2))
					}
					y_temp = as.vector(y.hat + gamma*(diag(n)-proj)%*%X[,ind] + sigma.hat.gamma*V%*%u)
					glmnet_object_cv = cv.glmnet(X[,-ind], y_temp-gamma*X[,ind], standardize=FALSE)
					if(lambda_cv_flag == -1){
						lambda_cv = glmnet_object_cv$lambda.min
					} else if(lambda_cv_flag == -2){
						lambda_cv = glmnet_object_cv$lambda.1se
					}
				}			
				beta = Variable(ncol(X)+1)
				obj = sum(( (y - gamma*X[,ind]) -cbind(rep(1,n), X)%*%beta)^2)/(2*n) + lambda_cv*p_norm(beta[-1],1)
				prob = Problem(Minimize(obj))
				beta.hat.stat = round(solve(prob,feastol=1e-10,reltol=1e-10,abstol=1e-10)$getValue(beta), digits = 9)
				stat = abs(beta.hat.stat[1+ind])

				right_tail = l.cdf_adjusted(stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'right', known_sigma = known_sigma)
				left_tail = l.cdf_adjusted(-stat, y, X, ind, gamma, lambda_cv, lambda, smoothed = smoothed, tail = 'left', known_sigma = known_sigma)
				pvals.right[j] = left_tail + right_tail
			}
			inds.right = which(pvals.right > 1-coverage)
			inds.right = unique(c(1,inds.right))
			if(max(inds.right)!= outer_grid.length){
				inds.right = c(inds.right, max(inds.right)+1)
			}
			ci.right = right_additional[inds.right[length(inds.right)]]
		} else{
			ci.right = gamma_range[inds[length(inds)]]
			warning('CI hit the upper limit. Try increasing the range?')
		}
	} else{
		temp = range( gamma_range[inds] )
		 ci.left = temp[1]
		 ci.right = temp[2]
	}
	return(c(ci.left, ci.right))
}


