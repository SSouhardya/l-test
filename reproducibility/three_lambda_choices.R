source('single_test.R')

par_number = 7
n = rep(100, par_number)
p = rep(50, par_number)
#p = n/2
s = rep(5,par_number)
a_vals = c(0, 1.3,2.3,3.3,4.3,5.3,6.3)
rho = rep(0,par_number)


times = 1000

run_multiple<-function(times, lambda_choice = -2){
	D = matrix(0, nrow = length(a_vals), ncol = 4)
	L = list(length = length(a_vals))
	D[,1] = a_vals
	for(i in 1:length(a_vals)){
		print(i)
		L[[i]] = compare_unconditional_tests(n[i],p[i],s[i],a_vals[i],a_vals[i],rho[i],lambda_choice,0.05,times)
		D[i,-1] = L[[i]][[1]]
	}
	D = as.data.frame(D)
	names(D) = c('a','LASSO','t_both', 't_one')
	return(D)
}


D1 = run_multiple(times, -2)
print('FIRST DONE')
D2 = run_multiple(times, -5)
print('SECOND DONE')
D3 = run_multiple(times, -8)
print('THIRD DONE')


saveRDS(D1,'unconditional_test_summary_our_lambda.list')
saveRDS(D2,'unconditional_test_summary_hat_lambda.list')
saveRDS(D3,'unconditional_test_summary_full_lambda.list')








