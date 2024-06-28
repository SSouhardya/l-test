temp = read.csv(paste0('HIV_ci_',1,'.csv'), header = TRUE)
temp = temp[which(temp[,1]!=-1),]
combined = temp
for(i in 2:313){
    temp = read.csv(paste0('HIV_ci_',i,'.csv'), header = TRUE)
    temp = temp[which(temp[,1]!=-1),]
    combined = rbind(combined, temp)
	}
write.csv(combined, 'all_widths.csv', row.names = FALSE)
