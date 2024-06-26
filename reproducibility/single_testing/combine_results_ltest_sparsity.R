total_runs = 500

D = 0

for(i in 1:total_runs){
    l = readRDS(paste0('single_testing_',i,'.list'))
    D = D+l[[1]]
}
D = D/total_runs
D = as.data.frame(D)
names(D) = c('a_vals','LASSO','t_both','t_one')

saveRDS(list(D),'unconditional_test_summary_new.list')