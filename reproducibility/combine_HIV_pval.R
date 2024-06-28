get_position <- function(x){
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)
}


load('/reproducibility/data/HIV_data.RData')

raw_disc = matrix(0, nrow = 16, ncol = 4)
gene_disc = raw_disc

for(i in 1:16){
    if(i %in% 1:7){
        drug_class = 'PI' # "PI"    "NRTI"  "NNRTI"
        j = i
    } else if(i %in% 8:13){
        drug_class = 'NRTI'
        j = i-7
    } else{
        drug_class = 'NNRTI'
        j = i-13
    }

    X = data[[drug_class]]$X
    Y = data[[drug_class]]$Y
    y = Y[[j]]
    
    
    y <- log(y)
    missing <- is.na(y)
    y <- y[!missing]
    X <- X[!missing,]
    X <- X[,colSums(X) >= 3]
    X <- X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
    genes <- colnames(X)
    
    dat = read.csv(paste0('pvals_hiv_unadjusted_',i,'.csv'))[,-1]
    
    genes_l = genes[which(dat[,1]<0.05)]
    genes_t = genes[which(dat[,2]<0.05)]
    genes_t_int = genes[which(dat[,3]<0.05)]
    
    disc_l = length(unique(get_position(genes_l)))
    disc_t = length(unique(get_position(genes_t)))
    disc_t_int = length(unique(get_position(genes_t_int)))
    disc_total = length(unique(get_position(genes)))
    
    
    raw_disc_l = length(unique((genes_l)))
    raw_disc_t = length(unique((genes_t)))
    raw_disc_t_int = length(unique((genes_t_int)))
    raw_disc_total = length(unique((genes)))
    
    print( c(raw_disc_l, raw_disc_t, raw_disc_t_int, raw_disc_total))
    
    raw_disc[i,] =  c(raw_disc_l, raw_disc_t, raw_disc_t_int, raw_disc_total)
    gene_disc[i,] =  c(disc_l, disc_t, disc_t_int, disc_total)
}

raw_disc = as.data.frame(raw_disc)
names(raw_disc) = c('l_test','t_test','t_test_int','total')

gene_disc = as.data.frame(gene_disc)
names(gene_disc) = c('l_test','t_test','t_test_int','total')

write.csv(raw_disc, 'raw_HIV_discoveries.csv', row.names = FALSE)
write.csv(gene_disc, 'gene_HIV_discoveries.csv', row.names = FALSE)
