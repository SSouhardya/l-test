std<-function(y){
	y = as.vector(y)
	m.y = mean(y)
	sd.y = sqrt(var(y)*(length(y) - 1)/length(y))
	return((y - m.y)/sd.y)
}


g<-function(v){	#normalizes a vector
	return(v/sqrt(sum(v^2))  )
}

qhaar <- function(q,n,lower.tail=TRUE,stoperr=FALSE, known_sigma = FALSE, sigma.hat = 1){
	#cdf for the first (or any) element of a random vector distributed uniformly on S_{n-1}, i.e., the (n-1)-dimensional unit spherical shell in ambient dimension n
	if(known_sigma){
		return(pnorm(q*sigma.hat, lower.tail = lower.tail))
	}
	if( (abs(q)>1) & stoperr){
		p = NA
		stop("impossible haar quantile")
	}else if(q>=1){
		p = 1
	}else if(q<=-1){
		p = 0
	}else if(q==0){
		p = 0.5
	}else if(q<0){
		p = pt(sqrt((n-1)/(1/q^2-1)),n-1,lower.tail=FALSE)
	}else{
		p = 1-pt(sqrt((n-1)/(1/q^2-1)),n-1,lower.tail=FALSE)
	}
	if(!lower.tail){p = 1-p}
	p
}
