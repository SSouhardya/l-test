
#Reproducing Figure 2
mgp = c(1.5,0.5,0)  
oma = c(1.25,0.25,0.25,0.25) 
mfrow = c(1,2) 
mar = c(3, 2.5, 0.1, 0.3) 
cex.axis = 0.84 
cex = 0.4 
mfrow = c(1,3)

dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
par(mgp = mgp, oma = oma, mfrow = mfrow, mar = mar)

plot_uncond('unconditional_test_summary_setting1.list', c(0,1.3,2.3,3.3,4.3,5.3,6.3), 'Amplitude', plot_position = -2, cex.axis = cex.axis*3/2, cex = cex*3/2, pch = 19)
plot_uncond('unconditional_test_summary_sparsity.list', c(1,5,10,15,20,30,40,50), 'Sparsity', plot_position = -2, ylim = c(0,1),  times = 200000, cex = cex*3/2,pch = 19, cex.axis = cex.axis*3/2)
plot_uncond('unconditional_test_summary_correlation.list', c(0,0.1,0.3,0.5,0.7,0.9), 'Correlation', plot_position = -2, ylim = c(0,1), cex.axis = cex.axis*3/2, cex = cex*3/2, pch = 19)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = cex.axis)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
pp = locator(1)
legend('bottom',legend = c('\u2113-test',expression(paste(italic(t),'-test')),expression(paste('1-sided ',italic(t),'-test'))),lwd = rep(2,3), col = c('blue','red','chartreuse3'), horiz = TRUE, bty='n', xpd = TRUE)






#Reproducing Figure 3
mgp = c(1.5,0.5,0)  
oma = c(1.25,0.25,0.25,0.25) 
mfrow = c(1,2)
mar = c(3, 2.5, 1, 0.3) 
cex.axis = 0.84 
cex = 0.4
xrange = c(6,10,16,20,26,30,36,40) 


dev.new(width=6,height=3,unit = 'in',noRStudioGD = TRUE)
par(mgp = mgp, oma = oma, mfrow = c(2,2), mar = mar)

string = paste0('robustness_ltest_size_','t',2,'.list')
plot_uncond(string, xrange, 'Sample size', plot_position = -2,ylim = c(0,0.2), skip_onesided=TRUE, cex.axis = cex.axis*3/(2.5), cex = cex*3/(2.5), pch = 19, main = 'Heavy-tailed error', ylab = 'Type I error')
abline(h = 0.05, col = 'black', lwd = 2, lty = 2)

string = paste0('robustness_ltest_size_','gamma',1,'.list')
plot_uncond(string, xrange, 'Sample size', plot_position = -2,ylim = c(0,0.2), skip_onesided=TRUE, cex.axis = cex.axis*3/(2.5), cex = cex*3/(2.5), pch = 19, main = 'Skewed error', ylab = 'Type I error')
abline(h = 0.05, col = 'black', lwd = 2, lty = 2)

string = paste0('robustness_ltest_size_','hetero1',8,'.list')
plot_uncond(string, xrange, 'Sample size', plot_position = -2, ylim = c(0,0.2), skip_onesided=TRUE,  cex.axis = cex.axis*3/(2.5), cex = cex*3/(2.5), pch = 19, main = 'Heteroskedastic error', ylab = 'Typr I error')
abline(h = 0.05, col = 'black', lwd = 2, lty = 2)

string = paste0('robustness_ltest_size_non_linear',4,'.list')
plot_uncond(string, xrange, 'Sample size', plot_position = -2, ylim = c(0,0.2), skip_onesided=TRUE, cex.axis = cex.axis*3/(2.5), cex = cex*3/(2.5), pch = 19, main = 'Non-linearity', ylab = 'Type I error')
abline(h = 0.05, col = 'black', lwd = 2, lty = 2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = cex.axis)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
pp = locator(1)
legend('bottom',legend = c('\u2113-test',expression(paste(italic(t),'-test'))),lwd = rep(2,2), col = c('blue','red'), horiz = TRUE, bty='n', xpd = TRUE)





#Reproducing Figure 4
mgp = c(1.5,0.5,0)  
oma = c(1.25,0.25,0.25,0.25) 
mfrow = c(1,2) 
mar = c(3, 2.5, 0.1, 0.3) 
cex.axis = 0.84 
cex = 0.4 
mfrow = c(1,3)

dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
par(mgp = mgp, oma = oma, mfrow = mfrow, mar = mar)

plot_CI_unadj('CI_summary_unadjusted1.list','Amplitude',c(0,1.3,2.3,3.3,4.3,5.3,6.3),6, se = 2, only_power = TRUE, cex.axis = cex.axis*3/2, cex = cex*3/2 ,pp=-2)
plot_CI_unadj('CI_summary_unadjusted_sparsity.list','Sparsity',c(1,5,10,15,20,30,40,50),6, se = 2, only_power = TRUE, cex.axis = cex.axis*3/2, cex = cex*3/2 ,pp=-2)
plot_CI_unadj('CI_summary_unadjusted_correlation.list','Correlation',c(0,0.1,0.3,0.5,0.7,0.9),18.5,only_power = TRUE, se = 2, pp = -2, cex.axis = cex.axis*3/2, cex = cex*3/2 )

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = cex.axis)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
pp = locator(1)
legend('bottom',legend =  c(paste0('\u2113','-CI'), expression(paste(italic(t),'-CI')), expression(paste('1-sided ', italic(t),'-CI'))),lwd = rep(2,3), col = c('blue','red','chartreuse3'), horiz = TRUE, bty='n', xpd = TRUE)




#Reproducing Figure 5
dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
plot_CI_adj('CI_summary_adjusted1.list',7,se = 2, cex = 0.4, cex.axis = 0.84*1.2, cex.axis.legend = 0.8)




#Reproducing Figure 6
dev.new(width=4,height=3,unit = 'in',noRStudioGD = TRUE)
cex.axis = 1
cex = 0.3
se = 2
par(mgp = c(1.5,0.5,0), oma = c(0.25,0.25,0.25,0.25), mfrow = c(1,1), mar = c(2.5, 2.5, 0.3, 0.3))
single_var_plot(2,ylim = c(-14,-1),string = 'variation_ltest_log.list', log = TRUE, cex.axis = cex.axis, cex=cex)





#Reproducing Figure 7
times = 1000
cex.axis = 1
cex = 0.3
se = 2

D1 = readRDS('unconditional_test_summary_our_lambda.list')
D2 = readRDS('unconditional_test_summary_hat_lambda.list')
D3 = readRDS('unconditional_test_summary_full_lambda.list')

xrange = c(0,1.3,2.3,3.3,4.3,5.3,6.3)

dev.new(width=4,height=3,unit = 'in',noRStudioGD = TRUE)

par(mgp = c(1.5,0.5,0), oma = c(0.25,0.25,0.25,0.25), mfrow = c(1,1), mar = c(2.5, 2.5, 0.3, 0.3))

plot(0, xlim = range(xrange), ylim = c(0,1), xlab = 'Amplitude', ylab = 'Power', col = 'white', cex.axis = cex.axis, cex.lab = cex.axis)
points(xrange, D1[,2], pch = 19, col = 'blue', type = 'l', lwd = 2, cex = cex)
points(xrange, D2[,2], pch = 19, col = 'brown', type = 'l', lwd = 2, cex = cex)
points(xrange, D3[,2], pch = 19, col = 'orange', type = 'l', lwd = 2, cex = cex)

points(xrange, D1[,2], pch = 19, col = 'blue', type = 'p', lwd = 2, cex = cex)
points(xrange, D2[,2], pch = 19, col = 'brown', type = 'p', lwd = 2, cex = cex)
points(xrange, D3[,2], pch = 19, col = 'orange', type = 'p', lwd = 2, cex = cex)

for(i in 1:nrow(D1)){
    sd = sqrt(D1[i,2]*(1-D1[i,2])/times)
    points(c(xrange[i], xrange[i]), c(D1[i,2] -se*sd, D1[i,2] + se*sd), pch = 19, col ='blue', type = 'l', lwd = 2, )

    sd = sqrt(D2[i,2]*(1-D2[i,2])/times)
    points(c(xrange[i], xrange[i]), c(D2[i,2] -se*sd, D2[i,2] + se*sd), pch = 19, col ='brown', type = 'l', lwd = 2)

    sd = sqrt(D3[i,2]*(1-D3[i,2])/times)
    points(c(xrange[i], xrange[i]), c(D3[i,2] -se*sd, D3[i,2] + se*sd), pch = 19, col ='orange', type = 'l', lwd = 2)
}

legend('bottomright', 
    legend = c('Choice 1', 'Choice 2', 'Choice 3'),
    lwd = c(2,2,2),
    col = c('blue', 'brown','orange' ))




#Reproducing Figure 8
dev.new(width=4,height=3,unit = 'in',noRStudioGD = TRUE)
cex.axis = 1
cex = 0.3
se = 2
par(mgp = c(1.5,0.5,0), oma = c(0.25,0.25,0.25,0.25), mfrow = c(1,1), mar = c(2.5, 2.5, 0.3, 0.3))
plot_uncond('unconditional_test_summary_setting2.list', c(0,1.3,2.3,3.3,4.3,5.3,6.3), 'Amplitude', ylim = c(0,1), plot_position = 'topleft',cex.axis = cex.axis, cex = cex)





#Reproducing figures 9-12
mgp = c(1.5,0.5,0)  
oma = c(1,0.25,0.25,0.25) 
mfrow = c(3,2) 
mar = c(3, 2.5, 2, 0.3)
cex.axis = 0.84
cex = 0.4
dev.new(width=6,height=5,unit = 'in',noRStudioGD = TRUE)
times = 1000
xlab = 'Sample size'


#un-comment either of these settings to draw the respective figure and then copy paste the following codes.

#For Figure 9

#dist_type = 't'
#par_vec = c(30,20,15,10,5,2)
#par_name = 'df'


#For Figure 10

#dist_type = 'gamma'
#par_vec = c(1,2,4,6,8,10)
#par_name = 'shape'

#For Figure 11

#dist_type = 'hetero1'
#par_vec = c(0.01,0.25,0.5,1,4,8)
#par_name = '\u03b7\u00b2'


#For Figure 12

#dist_type = 'non_linear'
#par_vec = c(0.3,0.5,1,2,3,4)
#par_name = '\u03b4'



xrange = c(6,10,16,20,26,30,36,40) 
par(mgp = mgp, oma = oma, mfrow = mfrow, mar = mar)
for(i in 1:length(par_vec)){
    string = paste0('robustness_ltest_size_',dist_type,par_vec[i],'.list')
    plot_uncond(string, xrange, xlab = xlab, ylab = 'Type I error',save_pdf = FALSE, ylim = c(0,0.2), height = 5, width = 5, plot_position = -2, main = paste0(par_name,' = ',par_vec[i]), skip_onesided = TRUE, cex = cex*1.35,pch = 19, cex.axis = cex.axis*1.35)
    abline(h = 0.05, col = 'black', lwd = 2, lty = 2)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = cex.axis)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
pp = locator(1)
legend('bottom',legend = c('\u2113-test',expression(paste(italic(t),'-test'))),lwd = rep(2,2), col = c('blue','red'), horiz = TRUE, bty='n', xpd = TRUE)





#Reproducing Figure 13
dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
plot_CI_unadj('CI_summary_unadjusted1.list','Amplitude',c(0,1.3,2.3,3.3,4.3,5.3,6.3),6,plot_1se = TRUE, se = 2,cex = 0.4, cex.axis = 0.84*1.2, cex.axis.legend=0.84)

#Reproducing Figure 14
dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
plot_CI_unadj('CI_summary_unadjusted2.list','Amplitude',c(0,1.3,2.3,3.3,4.3,5.3,6.3),17,plot_1se = TRUE, se = 2,cex = 0.4, cex.axis = 0.84*1.2, cex.axis.legend=0.84)

#Reproducing Figure 15
dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
plot_CI_unadj('CI_summary_unadjusted_sparsity.list','Sparsity',c(1,5,10,15,20,30,40,50),7,plot_1se = TRUE, se = 2, cex = 0.4, cex.axis = 0.84*1.2, cex.axis.legend=0.84)

#Reproducing Figure 16
dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
plot_CI_unadj('CI_summary_unadjusted_correlation.list','Correlation',c(0,0.1,0.3,0.5,0.7,0.9),18.5,plot_1se = TRUE, se = 2, pp = 'topleft',cex = 0.4, cex.axis = 0.84*1.2, cex.axis.legend=0.84)


#Reproducing Figure 17
dev.new(width=6,height=2.5,unit = 'in',noRStudioGD = TRUE)
plot_CI_adj('CI_summary_adjusted2.list',18, se = 2, cex = 0.4, cex.axis = 0.84*1.2, cex.axis.legend=0.84, string2 = 'CI_summary_adjusted2_known_sigma.list')

