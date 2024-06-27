#may need to change these depending on the location of the files.
source('~/l_testing.R')
source('liu_ci.R')
source('~adjusted_ci.R')

# TESTS FOR COMPARING LENGTH AND COVERAGE OF THE VARIOUS CIs
ci_simulations<-function(bx_range, n=60, p = 30, ind = 1,intercept = 0,s=5, rho=0,lambda = 0.02, m=3, coverage= 0.95, r= 2, margin = 10, sigma, normalize_cols = TRUE, adjusted=TRUE, mode = 1, known_sigma = FALSE, extra_counter = '1'){
	require(CVXR)
	require(glmnet)
	require(MASS)
	lambda_flag = lambda
	if(adjusted){	
		temp_table_abs = vector(length = 2)	#coverage and length respectively
		temp_table_new = vector(length = 2)
		temp_table_liu = vector(length = 2)
		temp_table_lee = vector(length = 2)
		for(i in 1:m){
			f<-function(bx){
				X = mvrnorm(n,rep(0,p),toeplitz(rho^(0:(p-1))))
				if(normalize_cols){
					X = apply(X,2,g)
				}

				amplitudes = rep(bx,s)
				random_ind = sample(1:p,s, replace = FALSE)
				beta = rep(0,p)
				beta[random_ind] = amplitudes
				ind = random_ind[1]

 				y = intercept + X%*%beta + rnorm(n)
 				y.std = std(as.vector(y))

 				value_range = seq(from = bx-margin, to = bx+margin, length.out = r)
 				ci_abs = range( l.ci_adjusted(y,X,ind, value_range, lambda = lambda, lambda_cv = -1, coverage = coverage, smoothed = TRUE, display = FALSE, outer_approx = TRUE,outer_grid.length = 20) )
  				ci_new = range( l.ci_adjusted(y,X,ind, value_range, lambda = lambda, lambda_cv = lambda, coverage = coverage, smoothed = TRUE, display = FALSE, outer_approx = TRUE,outer_grid.length = 20) )

  				if(length(ci_abs) == 0){
  				    print('CV lambda length 0')
  				}
  				if(length(ci_new) == 0){
  				    print('same lambda length 0')
  				}
 				ci_liu = liu_full_interval(y,X,ind,lambda, coverage, c(bx-10^5, bx+10^5), sigma)	#the second to last argument is the interval for uniroot to search a root
 				if(ci_abs[1] == -1){
 					temp_table_abs[1] = temp_table_abs[2] = -1
 				}
 				else{

 					temp_table_abs[1] = as.numeric(ci_abs[1]<= bx & bx<= ci_abs[2])
 					temp_table_abs[2] = ci_abs[2] - ci_abs[1] 
 				}

 				if(ci_new[1]==-1){
 					temp_table_new[1] = temp_table_new[2] = -1
 				}
 				else{

 					temp_table_new[1] = as.numeric(ci_new[1] <= bx & bx <= ci_new[2])
 					temp_table_new[2] = ci_new[2] - ci_new[1] 
 				}
  				if(ci_liu[1]==-1){
 					temp_table_liu[1] = temp_table_liu[2] = -1
 				}
 				else{

 					temp_table_liu[1] = as.numeric(ci_liu[1] <= bx & bx <= ci_liu[2])
 					temp_table_liu[2] = ci_liu[2] - ci_liu[1] 
	 			}
	 			return(c(temp_table_abs, temp_table_new, temp_table_liu))
	 		}
 			result_table = t(sapply(bx_range,f))
 			result_table = as.data.frame(cbind(bx_range, result_table) )
 			names(result_table) = c('Amplitude', 'coverage_cv','length_cv', 'coverage_lambda', 'length_lambda', 'coverage_liu', 'length_liu')
 			write.csv(result_table, paste(extra_counter,'CI_results',slurm_ind,i,'.csv'), row.names = FALSE)
 			print(paste('Iteration ',i,' complete.'))
 		}
	}
	else{
		s = length(bz_vector)-1
		table_lci.min = vector(length = 2)	#coverage and length respectively
		table_lci.1se = vector(length = 2)	#coverage and length respectively
		table_tci = vector(length = 2)
		table_tci_target = vector(length = 2)
		for(i in 1:m){
			f<-function(bx){
				X = mvrnorm(n,rep(0,p),toeplitz(rho^(0:(p-1))))
				if(normalize_cols){
					X = apply(X,2,g)
				}

				amplitudes = rep(bx,s)
				random_ind = sample(1:p,s, replace = FALSE)
				beta = rep(0,p)
				beta[random_ind] = amplitudes
				ind = random_ind[1]

 				y = intercept + X%*%beta + rnorm(n)
				
 				value_range = seq(from = bx-margin,to = bx+margin, length.out = r)
 				temp = l.ci(y,X,ind, value_range, lambda = -1, lambda_cv = -1, coverage = coverage, adjusted = FALSE, display = FALSE, outer_approx = TRUE)

 				ci_l.min = range(temp)

 				temp = l.ci(y,X,ind, value_range, lambda = -2, lambda_cv = -2, coverage = coverage, adjusted = FALSE, display = FALSE, outer_approx = TRUE)
 				ci_l.1se = range(temp)

 				table_lci.min[1] =  as.numeric(ci_l.min[1] <= bx & bx <= ci_l.min[2])
 				table_lci.min[2] =  ci_l.min[2] - ci_l.min[1]


 				table_lci.1se[1] =  as.numeric(ci_l.1se[1] <= bx & bx <= ci_l.1se[2])
 				table_lci.1se[2] =  ci_l.1se[2] - ci_l.1se[1]

				indic1 = vector(length = length(value_range))	#this is for the confidence interval that knows the true sign
				indic2 = indic1	#this is for the two-sided confidence interval

				for(i in 1:length(value_range)){
 					val = value_range[i]
 					z = y - val*X[,ind]
 					z.std = std(z)
 		 			fit = lm(z.std~X)
			    	est = summary(fit)$coef[1+ind,1]
			    	se = summary(fit)$coef[1+ind,2]
  					tstat = (est)/se
					indic2[i] = as.numeric( 2*pt(abs(tstat), df = n-p-1, lower.tail = FALSE)> 1-coverage)
 					indic1[i] = as.numeric( pt(tstat, df = n-p-1, lower.tail = (bx<val)) > 1-coverage)
				}
                inds2 = which(indic2 == 1)
	        	if(min(inds2)!=1){
		        	inds2 = c(min(inds2)-1, inds2)
	        	}
	        	if(max(inds2) != length(value_range) ){
	        		inds2 = c(inds2, max(inds2)+1)
	        	}
 	  			ci_t = range(value_range[inds2])
  				table_tci[1] = as.numeric(ci_t[1] <= bx & bx <= ci_t[2])
 				table_tci[2] = ci_t[2] - ci_t[1] 

                inds1 = which(indic1 == 1)
	        	if(min(inds1)!=1){
		        	inds1 = c(min(inds1)-1, inds1)
	        	}
	        	if(max(inds1)!= length(value_range) ){
	        		inds1 = c(inds1, max(inds1)+1)
	        	}
 	  			ci_t = range(value_range[inds1])
 	  			
  				table_tci_target[1] = as.numeric(ci_t[1] <= bx & bx <= ci_t[2])
 				table_tci_target[2] = ci_t[2] - ci_t[1] 

 				return(c(table_lci.min, table_lci.1se, table_tci, table_tci_target))
 			}
  			result_table = t(sapply(bx_range,f))
 			result_table = as.data.frame(cbind(bx_range, result_table) )
 			names(result_table) = c('Amplitude', 'coverage_lci_min','length_lci_min', 'coverage_lci_1se', 'length_lci_1se', 'coverage_t', 'length_t', 'coverage_t_target', 'length_t_target')
 			write.csv(result_table, paste(extra_counter,'CI_results',slurm_ind,i,'.csv'), row.names = FALSE)
 			print(paste('Iteration ',i,' complete.'))
		}
	}
}


#input to the above function is a parameter file that should be in the form of a list of lists.
#each of the individual list should contain two lists. The first list should be the values of c(n,p,ind,intercept, rho, lambda, m, coverage, r, margin, sigma, normalize_cols, adjusted, mode, known_sigma) in this order where they mean:

#n: number of rows of X
#p: the dimension of X
#ind: The index (without considering the intercept) of the coefficient which we are inferrring upon.
#intercept: A logical variable. TRUE if we want to include the intercept.
#rho: Intervariable correlations (in a toeplitz design)
#lambda: The selection lambda (of no use if adjusted = TRUE)
#m: Number of replications in a single run
#coverage: The coverage of the intervals
#r: Length of the grid of hypothesized values on which the \ell-test would be inverted.
#margin: The length of the half interval from which grid of values is selected.
#sigma: Error standard deviation
#normalize_cols = Whether the columns of the design matrix should be normalized
#adjusted: TRUE indicates that we want an adjusted confidence interval
#known_sigma: Whether we should run the \ell-test with a known-sigma. This is only used with the adjusted confidence interval.

#The second list should be the amplitude values to try.


#the parameter file should be stored in an .RData file created using saveRDS function.



slurm_ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(slurm_ind)



par_filename = commandArgs(trailingOnly=TRUE)
par_list = readRDS(par_filename)

par_length = length(par_list)


for(iter in 1:par_length){
	li = par_list[[iter]]
	par_vec = li[[1]]
	bx_range = li[[2]]

	n = par_vec[1]
	p = par_vec[2]
	ind = par_vec[3]
	intercept = par_vec[4]
	s = par_vec[5]
	rho = par_vec[6]
	lambda = par_vec[7]
	m = par_vec[8]
	coverage = par_vec[9]
	r = par_vec[10]
	margin = par_vec[11]
	sigma = par_vec[12]
	normalize_cols = par_vec[13]
	adjusted = par_vec[14]
	mode = par_vec[15]
	known_sigma = par_vec[16]	#is of need only for the adjusted-CI case

	if(adjusted && (lambda == -1)){
		stop('Lambda has to be fixed to obtain adjusted CI')
	}

	ci_simulations(bx_range,n,p,ind, intercept, s, rho, lambda, m, coverage, r, margin, sigma, normalize_cols, adjusted, mode, known_sigma, as.character(iter))
}
