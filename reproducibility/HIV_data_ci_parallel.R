source('~/l_testing.R')

ttest_bh_nointercept<-function(X,y,q, method = 'bh'){
    pvals = summary(lm(y~X-1))$coef[,4]
    if(method == 'bh'){
        rejections = bh(pvals,q)
    } else if(method == 'holm'){
        rejections = holm_bonferroni(pvals,q)
    } else if(method == 'no_mult'){
        rejections = which(pvals<q)
    }
    return(rejections)
}




HIV_ci <- function(X, y, alpha, gamma_length, col_ind, gamma_range = NULL){ 
    library(glmnet)
    ## Log-transform the drug resistance measurements.
    y <- log(y)

    ## Remove patients with missing measurements.
    missing <- is.na(y)
    y <- y[!missing]
    X <- X[!missing,]

    ## Remove predictors that appear less than 3 times.
    X <- X[,colSums(X) >= 3]

    ## Remove duplicate predictors.
    X <- X[,colSums(abs(cor(X)-1) < 1e-4) == 1]


    y.std = std(y)

    ## Get names
    genes <- colnames(X)


    ## Get stats
    
    n <- nrow(X)
    p <- ncol(X)
    coverage = 1-alpha
    
    if(col_ind>p){
        return(c(-1,-1))
    }

    W = apply(X,2,g)

    
    fit = lm(y~W)
    
   
    j = col_ind
    
    ci_t = as.vector(confint(fit, level = coverage)[1+j,])
    mid = mean(ci_t)
    half_width = (ci_t[2] - ci_t[1])/2
    if(sum(is.null(gamma_range))>0){
        gamma_range = seq(from = mid-10*half_width, to = mid + 10*half_width, length.out = gamma_length)
    }
    ci_l = range(l.ci(y,W,j,gamma_range = gamma_range, coverage = coverage, outer_approx = TRUE, outer_grid.length = gamma_length)
    l_length = ci_l[2] - ci_l[1]
    t_length = ci_t[2] - ci_t[1]
    
    return(c(l_length, t_length))
}


slurm_ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(2023)
col_ind = slurm_ind

display = FALSE
alpha = 0.05
gamma_length = 100


load('~/reproducibility/data/HIV_data.RData')

l = matrix(0, nrow = 16, ncol = 2)
l = as.data.frame(l)
names(l) = c('l','t')
index = 1

for (drug_class in names(data)){
  print(drug_class)
  X = data[[drug_class]]$X
  Y = data[[drug_class]]$Y
  col_len = ncol(Y)
  for(i in 1:col_len){
    print(paste0('Index: ',i))
    if(drug_class == 'NRTI' & i == 1 & slurm_ind == 270){
        gamma_range = seq(from = 80, to = 88, length.out = gamma_length/2)
    } else{
        gamma_range = NULL
    }
    l[index,] = HIV_ci(X, Y[[i]], alpha, gamma_length,col_ind,gamma_range)
    index = index + 1
  }
}

write.csv(l, paste0('HIV_ci_',slurm_ind,'.csv'), row.names = FALSE )




