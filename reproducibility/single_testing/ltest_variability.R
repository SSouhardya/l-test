slurm_ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(slurm_ind)

source('~/l_testing.R')


single_test_variability<-function(a){
n = 100
p = 50
m = 100
ind = 1

pvals = matrix(0, nrow = m, ncol = m)
lambda_vals = pvals

for(i in 1:m){
	print(i)
	X = matrix(rnorm(n*p), nrow = n)
	X = apply(X,2,g)
	beta = c(rep(a,5), rep(0,p-5))
	y = as.vector(X%*%beta + rnorm(n))
	Z = cbind(rep(1,n), X[,-1])
 	proj = Z%*%solve(t(Z)%*%Z,t(Z))
 	y.hat = proj%*%y
 	sigma.hat = sqrt(sum(((diag(n)-proj)%*%(y))^2))
 	V = qr.Q(qr(diag(n)-proj))[,1:(n-ncol(Z))]
	for(j in 1:m){
		#print(c(i,j))
 		 u = rnorm(n-ncol(Z))
   		 u = u/sqrt(sum(u^2))
    	 y_new = y.hat + sigma.hat*V%*%u
    	 lambda = cv.glmnet( X,y_new,standardize=FALSE)$lambda.min
    	 lambda_vals[i,j] = lambda
    	 pvals[i,j] = l.test_glmnet(y,X,ind, lambda= lambda, lambda_cv= lambda, adjusted = FALSE, smoothed = TRUE)
	}
}
pvals1 = pvals

u = var(as.vector(pvals))*(m^2-1)
v = sum(apply(pvals,1,var))*(m-1)
return(c(a,v,u))
}


arange = c(1.3,2.3,3.3,4.3,5.3,6.3)

variation_summary = as.data.frame(matrix(0, nrow = length(arange), ncol = 3))
names(variation_summary) = c('Amplitude', 'Within_variation', 'Total_variation')
for(i in 1:length(arange)){
	variation_summary[i,] = single_test_variability(arange[i])
	print(paste('Simulation ',i,'complete'))
}

write.csv(variation_summary,paste0('variation_single_test_',slurm_ind,'.csv'), row.names = FALSE)
