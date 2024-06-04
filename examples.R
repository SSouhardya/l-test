set.seed(1)

n = 100
p = 50
s =  5
A = 2.3

X = matrix(rnorm(n*p), nrow = n)
X = apply(X,2,g)
beta = rep(0,p)
rand_ind = sample(1:p, size = s, replace = FALSE)
j = rand_ind[1]	#index to test
beta[rand_ind] = (1-2*rbinom(s,1,0.5))*A
y = as.numeric(X%*%beta + rnorm(n)) #data



pval_l = l.test(y,X,j)	#l-test for H_j:\beta_j = 0 with lambda chosen using cross-validation

pval_l = l.test(y-2.3*X[,j],X,j) #for testing H_j(2.3):\beta_j = 2.3

pval_l = l.test(y,X,j, lambda_cv = 0.01) #l-test for H_j:\beta_j = 0 with a supplied lambda

pval_l_adjusted = l.test(y,X,j, adjusted = TRUE, lambda = 0.01) #adjusted l-test for H_j:\beta_j = 0 valid conditionally on LASSO selection using penalty 0.01, and the penalty for the test statistic chosen using cross-validation	



gamma_range = seq(from = beta[j]-10, to = beta[j]+10, length.out = 100) #the grid of \gamma values to test on

ci_l = l.ci(y,X,j, gamma_range = gamma_range, coverage = 0.95) #l-CI

ci_l_adjusted = l.ci(y,X,j, gamma_range = gamma_range, coverage = 0.95, lambda = 0.01, adjusted = TRUE) #post-selection l-CI for \beta_j valid conditionally on LASSO with penalty 0.01 selecting the coefficient



