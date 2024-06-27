
plot_uncond<-function(string, xrange, xlab, ylab = 'Power', save_pdf = FALSE, height = 5, width = 5, plot_position = 'topleft', main = NULL, ylim = NULL, skip_onesided = FALSE, se = 1, times = 200, pch = 19, cex = 1, cex.axis = 2, rm_y = FALSE){
	ll = readRDS(string)
	tab = ll[[1]]

    sd = tab
    sd = sqrt(sd*(1-sd)/times)
    sd[,1] = tab[,1]

	if(save_pdf){
			pdf(file = paste0("/Users/Souhardya/Desktop/Presentation_figures/uncond.pdf"),   # The directory you want to save the file in
   		 width = width, # The width of the plot in inches
   		 height = height)
	}
	amp = tab[,1]
    if(is.null(ylim)){
        if(!skip_onesided){
             ylim = range(cbind(tab[,-1]+se*sd[,-1], tab[,-1]-se*sd[,-1]))
        } else{
             ylim = range(cbind(tab[,-c(1,4)]+se*sd[,-c(1,4)], tab[,-c(1,4)]-se*sd[,-c(1,4)]))   
        }
    }
    if(rm_y){
	   plot(0, xlim = range(xrange), ylim = ylim, xlab = xlab, ylab = ylab, col = 'white', main = main,cex.axis = cex.axis, cex.lab = cex.axis, yaxt = 'n')# xaxt = 'n')
    } else{
        plot(0, xlim = range(xrange), ylim = ylim, xlab = xlab, ylab = ylab, col = 'white', main = main,cex.axis = cex.axis, cex.lab = cex.axis)# yaxt = 'n', xaxt = 'n')
    }
    if(plot_position == -1){
        if(!skip_onesided){
            pp = locator(1)
          legend(pp$x, pp$y,
              legend = c('\u2113-test',expression(paste(italic(t),'-test')),expression(paste('1-sided ',italic(t),'-test'))),
             col = c('blue','red','chartreuse3'),
             lwd = c(2,2,2)
             )
        } else{
            pp = locator(1)
          legend(pp$x, pp$y,
              legend = c(paste0('\u2113','-test'), '2-sided t-test'),
             col = c('blue','red'),
             lwd = c(2,2)
             )
        }
    } else if(plot_position == -2){
    }else{
        if(!skip_onesided){
          legend(plot_position,
              legend = c('\u2113-test',expression(paste(italic(t),'-test')),expression(paste('1-sided ',italic(t),'-test'))),
             col = c('blue','red','chartreuse3'),
             lwd = c(2,2,2)
             )
        } else{
          legend(plot_position,
              legend = c(paste0('\u2113','-test'), '2-sided t-test'),
             col = c('blue','red'),
             lwd = c(2,2)
             )
        }
    }

	points(xrange, tab[,2], type = 'l', lwd = 2, col = 'blue',pch = pch, cex = cex)
	points(xrange, tab[,3], type = 'l', lwd = 2, col = 'red',pch = pch, cex = cex)
    if(!skip_onesided){
	   points(xrange, tab[,4], type = 'l', lwd = 2, col = 'chartreuse3',pch = pch, cex = cex)
    }

    points(xrange, tab[,2], type = 'p', lwd = 2, col = 'blue',pch = pch, cex = cex)
    points(xrange, tab[,3], type = 'p', lwd = 2, col = 'red',pch = pch, cex = cex)
    if(!skip_onesided){
       points(xrange, tab[,4], type = 'p', lwd = 2, col = 'chartreuse3',pch = pch, cex = cex)
    }

    for(i in 1:nrow(tab)){
        #print(c(xrange[i],xrange[i]))
        points(c(xrange[i],xrange[i]), c(tab[i,2] - se*sd[i,2],tab[i,2] + se*sd[i,2]), type = 'l', lwd = 2, col = 'blue')
        points(c(xrange[i],xrange[i]), c(tab[i,3] - se*sd[i,3],tab[i,3] + se*sd[i,3]), type = 'l', lwd = 2, col = 'red')
        if(!skip_onesided){
             points(c(xrange[i],xrange[i]), c(tab[i,4] - se*sd[i,4],tab[i,4] + se*sd[i,4]), type = 'l', lwd = 2, col = 'chartreuse3')
        }
    }

	if(save_pdf){
		dev.off()
	}
}







plot_CI_unadj<-function(string, xname, xrange, ylim, save_pdf = FALSE, width = 5, height = 5, only_power = FALSE, plot_1se = FALSE, se = 1, color = c('blue','lightblue','red','chartreuse3'), pp = 'bottomleft', pch = 19, cex = 1, cex.axis = 1, cex.axis.legend = 1){
	bundle = readRDS(string)
	tab = bundle[[1]]
	tab[,1] = xrange
	std_error = bundle[[2]]
	std_error[,1] = xrange
	amplitude = xrange
	if(save_pdf){
			pdf(file = paste0("/Users/Souhardya/Desktop/Presentation_figures/Unadj_CI.pdf"),   # The directory you want to save the file in
    width = width, # The width of the plot in inches
    height = height)
	}
	if(!only_power){
		par(mgp = c(1.5,0.5,0), oma = c(1,0.25,0.25,0.25), mfrow = c(1,2), mar = c(3, 2.5, 0.1, 0.3), cex = cex.axis.legend)
		plot(0, xlim = range(amplitude), ylim = c(0,1), xlab = xname, ylab = 'Coverage', col = 'white',cex.axis = cex.axis, cex.lab = cex.axis)
		#points(amplitude, tab[,2], pch = 19, col = 'red', type = 'b', lwd = 2)
		points(amplitude, tab[,2], pch = pch, cex = cex, col = color[1], type = 'p', lwd = 2)
		points(amplitude, tab[,2], pch = pch, cex = cex, col = color[1], type = 'l', lwd = 2)
		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c( tab[i,2]-se*std_error[i,2],tab[i,2]+se*std_error[i,2]), pch = pch, cex = cex, col = color[1], type = 'l', lwd = 2)
		}
		points(amplitude, tab[,6], pch = pch, cex = cex,  col = color[3], type = 'p', lwd = 2)
		points(amplitude, tab[,6], pch = pch, cex = cex,  col = color[3], type = 'l', lwd = 2)
		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c(tab[i,6]-se*std_error[i,6],tab[i,6]+se*std_error[i,6]), pch = pch, cex = cex, col = color[3], type = 'l', lwd = 2)
		}
		points(amplitude, tab[,8], pch = pch, cex = cex, col = color[4], type = 'p', lwd = 2)
		points(amplitude, tab[,8], pch = pch, cex = cex, col = color[4], type = 'l', lwd = 2)
		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c(tab[i,8]-se*std_error[i,8],tab[i,8]+se*std_error[i,8]), pch = pch, cex = cex, col = color[4], type = 'l', lwd = 2)
		}
		if(plot_1se){
			points(amplitude, tab[,4], pch = pch, cex = cex, col = color[2], type = 'p', lwd = 2)
			points(amplitude, tab[,4], pch = pch, cex = cex, col = color[2], type = 'l', lwd = 2)
			for(i in 1:length(amplitude)){
				points(c(amplitude[i],amplitude[i]), c(tab[i,4]-se*std_error[i,4],tab[i,4]+se*std_error[i,4]), pch = pch, cex = cex,  col = color[2], type = 'l', lwd = 2)
			}
		}
		abline(h = 0.95, col = 'black', lwd = 2, lty =2)
		#p = locator(1)
	} 

	plot(0, xlim = range(amplitude), ylim = c(0,ylim), xlab = xname, ylab = 'Length', col = 'white',cex.axis = cex.axis, cex.lab = cex.axis)
	#points(amplitude, tab[,3], pch = 19, col = 'red', type = 'b', lwd = 2)
	points(amplitude, tab[,3], pch = pch,cex = cex,  col = color[1], type = 'p', lwd = 2)
	points(amplitude, tab[,3], pch = pch,cex = cex,  col = color[1], type = 'l', lwd = 2)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c(tab[i,3]-se*std_error[i,3],tab[i,3]+se*std_error[i,3]), pch = pch, cex = cex, col = color[1], type = 'l', lwd = 2)
	}
	points(amplitude, tab[,7], pch = pch, cex = cex, col = color[3], type = 'p', lwd = 2)
	points(amplitude, tab[,7], pch = pch, cex = cex, col = color[3], type = 'l', lwd = 2)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c(tab[i,7]-se*std_error[i,7],tab[i,7]+se*std_error[i,7]), pch = pch, cex = cex, col = color[3], type = 'l', lwd = 2)
	}
	points(amplitude, tab[,9], pch = pch, cex = cex, col = color[4], type = 'p', lwd = 2)
	points(amplitude, tab[,9], pch = pch, cex = cex, col = color[4], type = 'l', lwd = 2)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c(tab[i,9]-se*std_error[i,9],tab[i,9]+se*std_error[i,9]), pch = pch, cex = cex, col = color[4], type = 'l', lwd = 2)
	}
	if(plot_1se){
		points(amplitude, tab[,5], pch = pch,cex = cex,  col = color[2], type = 'p', lwd = 2)
		points(amplitude, tab[,5], pch = pch,cex = cex,  col = color[2], type = 'l', lwd = 2)
		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c(tab[i,5]-se*std_error[i,5],tab[i,5]+se*std_error[i,5]), pch = pch,cex = cex,  col = color[2], type = 'l', lwd = 2)
		}
	}
	#p = locator(1)

	if(pp != -2){
	if(plot_1se){
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = 0.84)
		plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
		pp = locator(1)
		legend('bottom',legend = c(paste0('\u2113','-CI (min)'), paste('\u2113','-CI (1se)'), expression(paste(italic(t),'-CI')), expression(paste('1-sided ', italic(t),'-CI'))),lwd = c(2,2,2,2), col = color, horiz = TRUE, bty='n', xpd = TRUE)
	} else{
		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = 0.84)
		plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
		pp = locator(1)
		legend('bottom',legend = c(paste0('\u2113','-CI (min)'), expression(paste(italic(t),'-CI')), expression(paste('1-sided ', italic(t),'-CI'))),lwd = c(2,2,2), col = color[-2], horiz = TRUE, bty='n', xpd = TRUE)
	}
	}

	if(save_pdf){
		dev.off()
	}	
}







plot_CI_adj<-function(string, ylim, se = 2, coverage = 0.95, color = c('blue','red','darkgreen'),cex = 1, cex.axis = 1, cex.axis.legend = 1, string2 = NULL){

	bundle = readRDS(string)
	tab = bundle[[1]]
	std_error = bundle[[2]]
	amplitude = tab[,1]

	if(!is.null(string2)){
		bundle1 = readRDS(string2)
		tab1 = bundle1[[1]]
		std_error1 = bundle1[[2]]
	}

	par(mgp = c(1.5,0.5,0), oma = c(1,0.25,0.25,0.25), mfrow = c(1,2), mar = c(3, 2.5, 0.1, 0.3), cex = cex.axis.legend)
	plot(0, xlim = range(amplitude), ylim = c(0,1), xlab = 'Amplitude', ylab = 'Coverage', col = 'white',cex.axis = cex.axis, cex.lab = cex.axis, cex = cex)
	points(amplitude, tab[,2], pch = 19, col = color[1], type = 'l', lwd = 2, cex = cex)
	points(amplitude, tab[,2], pch = 19, col = color[1], type = 'p', lwd = 2, cex = cex)

	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c( tab[i,2]-se*std_error[i,2],tab[i,2]+se*std_error[i,2]), pch = 19, col = color[1], type = 'l', lwd = 2, cex = cex)
	}
	points(amplitude, tab[,4], pch = 19, col = color[2], type = 'l', lwd = 2, cex = cex)
	points(amplitude, tab[,4], pch = 19, col = color[2], type = 'p', lwd = 2, cex = cex)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c( tab[i,4]-se*std_error[i,4],tab[i,4]+se*std_error[i,4]), pch = 19, col = color[2], type = 'l', lwd = 2, cex = cex)
	}
	points(amplitude, tab[,6], pch = 19, col = color[3], type = 'l', lwd = 2, cex = cex)
	points(amplitude, tab[,6], pch = 19, col = color[3], type = 'p', lwd = 2, cex = cex)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c( tab[i,6]-se*std_error[i,6],tab[i,6]+se*std_error[i,6]), pch = 19, col = color[3], type = 'l', lwd = 2, cex = cex)
	}
	if(!is.null(string2)){
		points(amplitude, tab1[,2], pch = 19, col = color[1], type = 'l',  lty=2,lwd = 2, cex = cex)
		points(amplitude, tab1[,2], pch = 19, col = color[1], type = 'p', lwd = 2, cex = cex)

		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c( tab1[i,2]-se*std_error1[i,2],tab1[i,2]+se*std_error1[i,2]), pch = 19, col = color[1], type = 'l', lwd = 2,  cex = cex)
		}
		points(amplitude, tab1[,4], pch = 19, col = color[2], type = 'l', lty=2, lwd = 2, cex = cex)
		points(amplitude, tab1[,4], pch = 19, col = color[2], type = 'p', lwd = 2, cex = cex)
		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c( tab1[i,4]-se*std_error1[i,4],tab1[i,4]+se*std_error1[i,4]), pch = 19, col = color[2], type = 'l', lwd = 2, cex = cex)
		}
	}
	abline(h = coverage, col = 'black', lwd = 2, lty =2)

	plot(0, xlim = range(amplitude), ylim = c(0,ylim), xlab = 'Amplitude', ylab = 'Length', col = 'white',cex.axis = cex.axis, cex.lab = cex.axis, cex = cex)

	points(amplitude, tab[,3], pch = 19, col = color[1], type = 'l', lwd = 2, cex = cex)
	points(amplitude, tab[,3], pch = 19, col = color[1], type = 'p', lwd = 2, cex = cex)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c( tab[i,3]-se*std_error[i,3],tab[i,3]+se*std_error[i,3]), pch = 19, col = color[1], type = 'l', lwd = 2, cex = cex)
	}
	points(amplitude, tab[,5], pch = 19, col = color[2], type = 'l', lwd = 2, cex = cex)
	points(amplitude, tab[,5], pch = 19, col = color[2], type = 'p', lwd = 2, cex = cex)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c( tab[i,5]-se*std_error[i,5],tab[i,5]+se*std_error[i,5]), pch = 19, col = color[2], type = 'l', lwd = 2, cex = cex)
	}
	points(amplitude, tab[,7], pch = 19, col = color[3], type = 'l', lwd = 2, cex = cex)
	points(amplitude, tab[,7], pch = 19, col = color[3], type = 'p', lwd = 2, cex = cex)
	for(i in 1:length(amplitude)){
		points(c(amplitude[i],amplitude[i]), c( tab[i,7]-se*std_error[i,7],tab[i,7]+se*std_error[i,7]), pch = 19, col = color[3], type = 'l', lwd = 2, cex = cex)
	}

	if(!is.null(string2)){
		points(amplitude, tab1[,3], pch = 19, col = color[1], type = 'l', lty=2, lwd = 2, cex = cex)
		points(amplitude, tab1[,3], pch = 19, col = color[1], type = 'p', lwd = 2, cex = cex)

		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c( tab1[i,3]-se*std_error1[i,3],tab1[i,3]+se*std_error1[i,3]), pch = 19, col = color[1], type = 'l', lwd = 2, cex = cex)
		}
		points(amplitude, tab1[,5], pch = 19, col = color[2], type = 'l', lty=2, lwd = 2, cex = cex)
		points(amplitude, tab1[,5], pch = 19, col = color[2], type = 'p', lwd = 2, cex = cex)
		for(i in 1:length(amplitude)){
			points(c(amplitude[i],amplitude[i]), c( tab1[i,5]-se*std_error1[i,5],tab1[i,5]+se*std_error1[i,5]), pch = 19, col = color[2], type = 'l', lwd = 2, cex = cex)
		}
	}

	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = 0.84)
	plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
	pp = locator(1)
	legend('bottom',legend = c('adjusted \u2113-CI', 'adjusted \u2113-CI with CV', 'Liu et al. (2018)'),lwd = c(2,2,2), col = color, horiz = TRUE, bty='n', xpd = TRUE)

}





single_var_plot_cluster<-function(se = 2, ylim = NA, m = 100, log = TRUE, pp='topleft', string = 'variation_sigle_test_non_log_list.list', cex.axis = 1, cex = 1){
	bundle = readRDS(string)
	mean_tab = bundle[[1]]
	sd_tab = bundle[[2]]
	mean_tab = mean_tab[,2:1]
	sd_tab = sd_tab[,2:1]
	if(log){
		mean_tab[,1] = 0.5*(mean_tab[,1] - log(m^2-1))
		mean_tab[,2] = 0.5*(mean_tab[,2] - log(m^2-m))
		sd_tab = 0.5*sd_tab
	}
	xlim = c(min(mean_tab[,1] - se*sd_tab[,1]), max(mean_tab[,1] + se*sd_tab[,1]))
	if(sum(is.na(ylim))){
		ylim = c(min(mean_tab[,2] - se*sd_tab[,2]), max(mean_tab[,2] + se*sd_tab[,2]))
	}

	xlim1 = c(min(c(xlim,ylim)), max(c(xlim,ylim)))	
	ylim1 = xlim1

	xlab = 'log(Overall s.d.)'
	ylab = 'log(Within-replicate s.d.)'

	if(!log){
		xlab = 'Overall standard deviation'
		ylab = 'Within-replicate standard deviation'
	}

	plot(0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, pch = 19, col = 'white',cex.axis = cex.axis, cex.lab = cex.axis)
	points(mean_tab[,1], mean_tab[,2], pch = 19, col = 'cyan4', type = 'p', cex = cex)
	points(mean_tab[,1], mean_tab[,2], pch = 19, col = 'cyan4', type = 'l', lwd = 2, cex = cex)
	for(i in 1:nrow(mean_tab)){
		points(c(mean_tab[i,1]-se*sd_tab[i,1], mean_tab[i,1]+se*sd_tab[i,1] ), c(mean_tab[i,2], mean_tab[i,2]), pch = 19, col = 'cyan4', type = 'l', lwd = 2, cex = cex)
		points(c(mean_tab[i,1], mean_tab[i,1]), c(mean_tab[i,2] - se*sd_tab[i,2], mean_tab[i,2] + se*sd_tab[i,2]), pch = 19, col = 'cyan4', type = 'l', lwd = 2, cex = cex)
	}
	points(c(ylim[1],xlim[2]), c(ylim[1],xlim[2]), lwd =2, col = 'black', lty = 2, type ='l', cex = cex)
	legend(pp,
		legend = c('Variation', 'Identity line'),
		col = c('cyan4', 'black'),
		lwd = c(2,2),
		lty = c(1,2)
	)
}











