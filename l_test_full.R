library(glmnet)
library(MASS)
library(CVXR)


#---------------- auxillary functions --------------------------------
std<-function(y){
	y = as.vector(y)
	m.y = mean(y)
	sd.y = sqrt(var(y)*(length(y) - 1)/length(y))
	return((y - m.y)/sd.y)
}


g<-function(v){
	return(v/sqrt(sum(v^2))  )
}


qhaar <- function(q,n,lower.tail=TRUE,stoperr=FALSE){
	#cdf for the first (or any) element of a random vector distributed uniformly on S_{n-1}, i.e., the (n-1)-dimensional unit spherical shell in ambient dimension n
	if( (abs(q)>1) & stoperr){
		p = NA
		stop("impossible haar quantile")
	}else if(q>=1){
		p = 1
	}else if(q<=-1){
		p = 0
	}else if(q==0){
		p = 0.5
	}else if(q<0){
		p = pt(sqrt((n-1)/(1/q^2-1)),n-1,lower.tail=FALSE)
	}else{
		p = 1-pt(sqrt((n-1)/(1/q^2-1)),n-1,lower.tail=FALSE)
	}
	if(!lower.tail){p = 1-p}
	p
}

#------------------ calculating beta cdfs ----------------------------

#this is cdf calculation using glmnet package. The null is always gamma = 0
l.cdf_glmnet<-function(x,y,X,ind,lambda,lambda_cv, glmnet_object=NULL, glmnet_object_type = 1, adjusted = FALSE, return_both = FALSE, return_dropprob = FALSE){
    # x is the point at which we want to evaluate the conditional cdf
    
    if(return_both){
        adjusted = TRUE
    }

    if(return_dropprob & (adjusted == FALSE) ){
    	adjusted = TRUE
    	warning('adjutsed set to TRUE since return_dropprob is TRUE')
    }
    
    n = nrow(X)
    p = ncol(X)
    
    y.std = std(y)
    Z = cbind(rep(1,n),X[,-ind])
    proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
    y.std.hat = proj %*%y.std
    
    sigma.hat.std = sqrt(sum((y.std - y.std.hat)^2))
    
    second_denom_term = sqrt(sum(((diag(n) - proj)%*%X[,ind])^2))
    
    if(is.null(glmnet_object)){
        glmnet_object_type = 1
        glmnet_object = list()
        glmnet_object[[1]] = cv.glmnet(X,y.std, standardize=FALSE)
        glmnet_object[[2]] = cv.glmnet(X[,-ind],y.std, standardize=FALSE)
    }
    
    ind_to_look = (glmnet_object_type-1)*(1+ind) + (2-glmnet_object_type)*2

    if(lambda == -1){
        lambda = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda == -2){
        lambda = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda == -3){
        lambda = glmnet_object[[1]]$lambda.1se
    } else if(lambda == -4){
        lambda = glmnet_object[[1]]$lambda.1se
    }
    
    if(lambda_cv == -1){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda_cv == -2){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda_cv == -3){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    } else if(lambda_cv == -4){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    }
    
    
    sgn<-function(x){
        if(x == 0){
            return(1)
        }
        return(sign(x))
    }
    beta.x = coef(glmnet(X[,-ind], y.std - x*X[,ind], standardize = FALSE), s = lambda_cv, exact = TRUE, x = X[,-ind], y = y.std- x*X[,ind])
    v.x = ( -sum(X[,ind]*(y.std.hat - x*X[,ind] - Z%*%beta.x)) +n*lambda_cv*sgn(x) )/(sigma.hat.std*second_denom_term)
    
    beta.0 = coef(glmnet_object[[ind_to_look]]$glmnet.fit, s = lambda, exact = TRUE, x = X[,-ind], y = y.std)
    v1 = ( -sum(X[,ind]*(y.std.hat - Z%*%beta.0)) -n*lambda )/(sigma.hat.std*second_denom_term)
    v2 = ( -sum(X[,ind]*(y.std.hat - Z%*%beta.0)) +n*lambda )/(sigma.hat.std*second_denom_term)

    denom = 1 - (qhaar(v2,n-ncol(X), lower.tail = TRUE) - qhaar(v1,n-ncol(X), lower.tail = TRUE) )	# i am calulating this outside as this give sthe dropping probability

    uncond_prob = qhaar(v.x,n-ncol(X),lower.tail = TRUE)
    
    if(!adjusted){
        return(uncond_prob)
    }
    
    beta.hat = coef(glmnet_object[[1]]$glmnet.fit, s=lambda, exact = TRUE, x = X, y = y.std)
    

    if(beta.hat[1+ind] == 0){
        cond_prob = -1
    } else{
        numer_1 = qhaar(v.x, n-ncol(X), lower.tail = TRUE)
        numer_2 = 0
        if(v.x>v1){
            numer_2 = qhaar(min(c(v.x,v2)), n- ncol(X), lower.tail = TRUE) - qhaar(v1, n- ncol(X), lower.tail = TRUE)
        }
        cond_prob = (numer_1 - numer_2)/denom
    }

    
    if(!return_both){
        toret = cond_prob
    } else{
        toret = c(uncond_prob, cond_prob)
    }
    if(return_dropprob){
    	toret = list(toret, 1-denom)
    }
    
    return(toret)
}


#--------------------------- l-testing codes ---------------------------------------------

l.test<-function(y,X,ind, lambda=-1, lambda_cv=-1, glmnet_object=NULL, glmnet_object_type = 1, adjusted = FALSE, smoothed = FALSE, return_both = FALSE){
	# y is the response vector
	# X is the data matrix WITHOUT the intercept column
    # ind is the index for which we want to test
    # the null is always gamma = 0
    # lambda is the lambda used for the selection lasso
    # lambda_cv is the one used for cross-validation
	# adjusted = TRUE performs test valid conditionally on \hat \beta_j^{\lambda}\neq 0
	# if adjusted = FALSE, lambda is not used.
    # lambda/lambda_cv = -1 is the default choice (and corresponds to the min rule)
    # lambda/lambda_cv = -2 corresponds to the 1se rule.
    # if glmnet_object_type == 1, then
    # glmnet_object is a glmnet object of length 2. The first one for the entire model while the second one for the model with ind dropped
    # if glmnet_object_type ==2, then
    # glmnet_object is a glmnet object of length 1+p. The first one is for the entire model while the ith one is when X_i is dropped
    # if return_both == TRUE, adjusted will be forced to be TRUE




	if(return_both){
		adjusted = FALSE
	}
	n = nrow(X)
    p = ncol(X)
    
    y.std = std(y)


    Z = cbind(rep(1,n),X[,-ind])
	proj = Z%*%solve(t(Z)%*%Z)%*%t(Z)
	y.std.hat = proj %*%y.std
	sigma.hat.std = sqrt(sum((y.std - y.std.hat)^2))
	V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
	u = rnorm(n-ncol(Z))
	u = u/sqrt(sum(u^2))
	y_temp = as.vector(y.std.hat + sigma.hat.std*V%*%u)
	y_temp.std = std(y_temp)

    if(is.null(glmnet_object)){
        glmnet_object_type = 1
        glmnet_object = list()
        glmnet_object[[1]] = cv.glmnet(X,y.std, standardize=FALSE)
        glmnet_object[[2]] = cv.glmnet(X[,-ind],y_temp.std, standardize=FALSE)
    }
    
    ind_to_look = (glmnet_object_type-1)*(1+ind) + (2-glmnet_object_type)*2

    if(lambda == -1){
        lambda = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda == -2){
        lambda = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda == -3){
        lambda = glmnet_object[[1]]$lambda.1se
    } else if(lambda == -4){
        lambda = glmnet_object[[1]]$lambda.1se
    }
    
    if(lambda_cv == -1){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.min
    } else if(lambda_cv == -2){
        lambda_cv = glmnet_object[[ind_to_look]]$lambda.1se
    } else if(lambda_cv == -3){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    } else if(lambda_cv == -4){
        lambda_cv = glmnet_object[[1]]$lambda.1se
    }

    x = abs(as.numeric(coef(glmnet_object[[1]]$glmnet.fit, s = lambda_cv, exact = TRUE, x = X, y = y.std)[1+ind]) )

    if(x == 0){
    	uncond_pval = 1
    	cond_pval = 1

    	if(!adjusted){
    		if(smoothed){

    			#first calculate u1
    			beta.hat.null = as.vector(coef(glmnet_object[[ind_to_look]]$glmnet.fit, s = lambda_cv, exact = TRUE, y = y_temp.std, x = X[,-ind]))
    			second_denom_term = sqrt(sum((X[,ind] - proj%*%X[,ind])^2))
    			mid = -as.numeric(X[,ind]%*%(y.std.hat - Z%*%beta.hat.null))/(sigma.hat.std * second_denom_term)
    			v1 = mid - n*lambda_cv/(sigma.hat.std * second_denom_term)
    			v2 = mid + n*lambda_cv/(sigma.hat.std * second_denom_term)

    			V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]

    			u1 = sum(X[,ind]*(y.std - y.std.hat))/(sigma.hat.std*second_denom_term)

    			dist_from_mid = abs(mid - u1)
    			w1 = mid - dist_from_mid
    			w2 = mid + dist_from_mid


    			#print(c(v1,w1,w2,v2, u1))
    			left_prob = qhaar(w1, n-ncol(Z), lower.tail = TRUE)
    			right_prob = 1 - qhaar(w2, n-ncol(Z), lower.tail = TRUE)

    			uncond_pval = left_prob + right_prob
    		}
    		if(!return_both){
    			return(uncond_pval)
    		}
    	}
    	beta.selection = abs(as.numeric(coef(glmnet_object[[1]]$glmnet.fit, s = lambda, exact = TRUE, x = X, y = y.std)[1+ind]) )
    	if(beta.selection == 0){
    		cond_pval = -1
    	}
    	if(!return_both){
    		return(cond_pval)
    	} else{
    		return(c(uncond_pval, cond_pval))
    	}
    } else{
    	pval.right = 1-l.cdf_glmnet(x,y,X,ind,lambda = lambda,lambda_cv = lambda_cv, glmnet_object = glmnet_object, glmnet_object_type = glmnet_object_type, adjusted = adjusted, return_both = return_both)
    	pval.left = l.cdf_glmnet(-x,y,X,ind,lambda = lambda,lambda_cv = lambda_cv, glmnet_object = glmnet_object, glmnet_object_type = glmnet_object_type, adjusted = adjusted, return_both = return_both)
    	pval = pval.left + pval.right
    	if(return_both){
    		if(pval.left[2] == -1){
    			pval[2] = -1
    		}
    	}
    	return(pval)
    }
}


l.ci<-function(y,X,ind, gamma_range, lambda=-1, lambda_cv=-1, coverage = 0.95, adjusted = FALSE, display = TRUE, smoothed = TRUE, outer_approx = FALSE, outer_grid.length = 10){
	# gamma_range is the grid of values to invert the l-test.

	gamma_range = sort(gamma_range)
	if(adjusted & (lambda<0)){
		stop('If adjusted, the selection lambda must be supplied')
	}
	g.length = length(gamma_range)
	pvals = vector(length = g.length)
	for(i in 1:g.length){
		y_test = y - gamma_range[i]*X[,ind]
		pvals[i] = l.test(y_test,X, ind, lambda, lambda_cv, adjusted = adjusted, smoothed = smoothed)
		if((i==1) & (pvals[1]==-1)){
			toc = Sys.time()
			if(display){
				print(toc - tic)
			}
			return(NA)	#indicating that there was no selection and hence, no adjusted CI.
		}
	}
	inds = which(pvals > 1-coverage)
	if(length(inds) == 0){
		warning('None of the values from the grid selected. Try a different grid?')
	}
	if(outer_approx){
		if(min(inds)!=1){
			inds1 = min(inds)-1
			left_additional = seq(from = gamma_range[inds1], to = gamma_range[min(inds)], length.out = outer_grid.length)
			pvals.left = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				y_test = y - left_additional[j]*X[,ind]
				pvals.left[j] = l.test(y_test,X, ind, lambda, lambda_cv, adjusted = adjusted, smoothed = smoothed)
			}
			inds.left = which(pvals.left > 1-coverage)
			if(min(inds.left)!=1){
				inds.left = c(min(inds.left)-1, inds.left)
			}
			if(inds.left[1] == Inf){    #which means empty vector
			    ci.left = left_additional[1]
			} else{
			    ci.left = left_additional[inds.left[1]]
			}
		} else{
			ci.left = gamma_range[inds[1]]
			warning('CI hit the lower limit. Try increasing the range?')
		}
		if(ci.left == gamma_range[1]){  
		    warning('CI hit the lower limit. Try increasing the range?')
		} 
		if(max(inds)!=g.length){
			inds2 = max(inds)+1
			right_additional = seq(from = gamma_range[max(inds)], to = gamma_range[inds2], length.out = outer_grid.length)
			pvals.right = vector(length = outer_grid.length)
			for(j in 1:outer_grid.length){
				y_test = y - right_additional[j]*X[,ind]
				pvals.right[j] = l.test(y_test,X, ind, lambda, lambda_cv, adjusted = adjusted, smoothed = smoothed)
			}
			inds.right = which(pvals.right > 1-coverage)
			if(max(inds.right)!= outer_grid.length){
				inds.right = c(inds.right, max(inds.right)+1)
			}
			if(inds.right[1] == -Inf){
			   ci.right = right_additional[1] 
			} else{
			    ci.right = right_additional[inds.right[length(inds.right)]]   
			}
		} else{
			ci.right = gamma_range[inds[length(inds)]]
			warning('CI hit the upper limit. Try increasing the range?')
		}
		if(ci.right == max(gamma_range)){ 
		    warning('CI hit the upper limit. Try increasing the range?')
		}
	} else{
		temp = range( gamma_range[inds] )
		 ci.left = temp[1]
		 ci.right = temp[2]
	}
	return(c(ci.left, ci.right))
}
