source('~/l_testing.R')


HIV_pvals <- function(X, y){    
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


    W = apply(X,2,g)

    
    p_bundle = vector(length = ncol(W))
    for(i in 1:ncol(W)){
        p_bundle[i] = l.test(y, W, i, smoothed = TRUE)
    }
    pval_l = p_bundle
    pval_t = summary(lm(y~X-1))$coef[,4]
    pval_t_int =  summary(lm(y~X))$coef[-1,4]
    
    dat = data.frame(pval_l, pval_t, pval_t_int)
    
    return(dat)
}



slurm_ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(slurm_ind)


display = FALSE


load('HIV_data.RData')


if(slurm_ind %in% 1:7){
    drug_class = 'PI' # "PI"    "NRTI"  "NNRTI"
    j = slurm_ind
} else if(slurm_ind %in% 8:13){
    drug_class = 'NRTI'
    j = slurm_ind-7
} else{
    drug_class = 'NNRTI'
    j = slurm_ind-13
}

X = data[[drug_class]]$X
Y = data[[drug_class]]$Y
y = Y[[j]]

dat_pvals = HIV_pvals(X,y)

write.csv(dat_pvals, paste0('pvals_hiv_unadjusted_',slurm_ind,'.csv'))
