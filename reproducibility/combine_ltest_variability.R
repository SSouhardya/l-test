total = 200
mean_tab = 0
sd_tab = 0
m = 200
for(i in 1:total){
	temp = (read.csv(paste0('variation_single_test_',i,'.csv'), header = TRUE)[,-1])
	temp[,1] = 0.5*log(temp[,1]/(m*(m-1)))
	temp[,2] = 0.5*log(temp[,2]/(m^2-1))
	mean_tab = mean_tab + temp
	sd_tab = sd_tab + temp^2
}
mean_tab = mean_tab/total
sd_tab = sqrt(sd_tab/total - mean_tab^2)/sqrt(total)
ll = list(mean_tab, sd_tab)
saveRDS(ll, 'variation_ltest_log.list')
