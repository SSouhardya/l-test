L = commandArgs(trailingOnly=TRUE) #parameter_filename, total_jobs, amplitude?, correlation/sparisty
L = strsplit(L,',')[[1]]

par_list = readRDS(L[1])

total_jobs = as.numeric(L[2])

instances = length(par_list)

l = list(length = 2)
l[[1]] = NA
l[[2]] = NA

for(i in 1:instances){
	par = par_list[[i]]
	m = par[[1]][8]	#this is number of times the for loop was executed
	adjusted = par[[1]][14]
	known_sigma = par[[1]][16]

    if(adjusted){
	    label = 'adjusted'
	} else{
	    label = 'unadjusted'
	}
	
	if(adjusted){
		bx_length = length(par[[2]])
		na_count_abs = rep(total_jobs*m, bx_length ) 
		na_count_both = na_count_abs
		na_count_liu = na_count_both
		tab_sum = tab_sq = 0

		for(a in 1:total_jobs){
			for(b in 1:m){
				tab = as.data.frame(read.csv(paste(i,'CI_results',a,b,'.csv'), header = TRUE))	

				na_count_abs[which(tab[,2] == -1 )] = na_count_abs[which(tab[,2] == -1 )] - 1
				na_count_both[which(tab[,4] == -1 )] = na_count_both[which(tab[,4] == -1 )] - 1
				na_count_liu[which(tab[,6] == -1 )] = na_count_liu[which(tab[,6] == -1 )] - 1

				tab[which(tab[,2] == -1 ), c(2,3)] = 0
				tab[which(tab[,4] == -1 ), c(4,5)] = 0
				tab[which(tab[,6] == -1 ), c(6,7)] = 0

				tab_sum = tab_sum + tab
				tab_sq = tab_sq + tab^2
		}
	}
	mean_tab = tab_sum
	for(k in 1:bx_length){
		mean_tab[k, c(2,3)] = mean_tab[k, c(2,3)]/na_count_abs[k]
		mean_tab[k, c(4,5)] = mean_tab[k, c(4,5)]/na_count_both[k]
		mean_tab[k, c(6,7)] = mean_tab[k, c(6,7)]/na_count_liu[k]
		#mean_tab[k, c(8,9)] = mean_tab[k, c(8,9)]/na_count_liu[k]
	}
	sd_tab = tab_sq
	for(k in 1:bx_length){
		sd_tab[k, c(2,3)] = sd_tab[k, c(2,3)]/na_count_abs[k]
		sd_tab[k, c(4,5)] = sd_tab[k, c(4,5)]/na_count_both[k]
		sd_tab[k, c(6,7)] = sd_tab[k, c(6,7)]/na_count_liu[k]
		#sd_tab[k, c(8,9)] = sd_tab[k, c(8,9)]/na_count_liu[k]
	}
	sd_tab[,c(2,3)] = sd_tab[,c(2,3)] - mean_tab[,c(2,3)]^2
	sd_tab[,c(4,5)] = sd_tab[,c(4,5)] - mean_tab[,c(4,5)]^2
	sd_tab[,c(6,7)] = sd_tab[,c(6,7)] - mean_tab[,c(6,7)]^2
	#sd_tab[,c(8,9)] = sd_tab[,c(8,9)] - mean_tab[,c(8,9)]^2

	for(k in 1:bx_length){
		sd_tab[k, c(2,3)] = sqrt(sd_tab[k, c(2,3)]/na_count_abs[k])
		sd_tab[k, c(4,5)] =sqrt(sd_tab[k, c(4,5)]/na_count_both[k])
		sd_tab[k, c(6,7)] = sqrt(sd_tab[k, c(6,7)]/na_count_liu[k])
		#sd_tab[k, c(8,9)] = sqrt(sd_tab[k, c(8,9)]/na_count_liu[k])
	}

	mean_tab[,1] = par[[2]]
	sd_tab[,1] = par[[2]]

	names(mean_tab) = c('Amplitude', 'coverage_lmin', 'length_lmin', 'coverage_same_cv', 'length_same_cv',  'coverage_liu', 'length_liu')
	names(sd_tab) = c('Amplitude', 'coverage_lmin', 'length_lmin', 'coverage_same_cv', 'length_same_cv', 'coverage_liu', 'length_liu')

	write.csv(mean_tab, paste(i,'mean_tab.csv'), row.names = FALSE)
	write.csv(sd_tab, paste(i,'sd_tab.csv'), row.names = FALSE)

	write.csv(data.frame(Amplitude = par[[2]], n_lmin = na_count_abs, n_same_cv = na_count_both, n_liu = na_count_liu), paste(i,'counts.csv'), row.names = FALSE)
	}
	else{
		tab_sum = 0
		tab_sq = 0
		#when things are not adjusted
		for(a in 1:total_jobs){
			for(b in 1:m){
				tab = as.data.frame(read.csv(paste(i,'CI_results',a,b,'.csv'), header = TRUE))
				tab_sum = tab_sum + tab
				tab_sq = tab_sq + tab^2
 			}
		}
		mean_tab = tab_sum/(total_jobs*m)
		sd_tab = tab_sq/(total_jobs*m) - mean_tab^2
		sd_tab = sqrt(sd_tab/(m*total_jobs))
		sd_tab[,1] = mean_tab[,1]

		write.csv(mean_tab, paste(i,'mean_tab.csv'), row.names = FALSE)
		write.csv(sd_tab, paste(i,'sd_tab.csv'), row.names = FALSE)
	}
	if(as.numeric(L[3])){
		l[[1]] = mean_tab
		l[[2]] = sd_tab
		if(!known_sigma){
		    saveRDS(l,paste0('CI_summary_',label,i,'.list'))
		 } else{
		    saveRDS(l,paste0('CI_summary_',label,i,'_known_sigma.list'))
		 }
	} else{
		l[[1]] =rbind(l[[1]], mean_tab)
		l[[2]] = rbind(l[[2]], sd_tab)
	}
}

if(!as.numeric(L[3])){
	l[[1]] = l[[1]][-1,]
	l[[2]] = l[[2]][-1,]
	saveRDS(l,paste0('CI_summary_',label,'_',L[4],'.list'))
}
