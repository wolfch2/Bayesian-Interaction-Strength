#############################################################################
#############################################################################
# Script to calculate per capita attack rates, assuming multispecies Type II functional response, using field-based estimates of species abundances,lab-based size & temperature-specific handling time regression coefficients, % feeding individuals from feeding surveys, and concurrent field temperatures.
# Specific to Tauranga Head Haustrum scobina data
#############################################################################
#############################################################################
# ~~~~~~~~~~~~~~~~
# Define functions
# ~~~~~~~~~~~~~~~~

# estimate handling time using lab-based regression coefficients
EstHTime<-function(htimes,feed.obs,PreySp){
	if(PreySp=='Not Feeding'){est<-NA}
	if(PreySp!='Not Feeding'){	
		coefs<-subset(htimes,Prey==PreySp)
		est <- exp(coefs$logIntC + coefs$logPredSizeC * log(feed.obs$PredSize) + coefs$logPreySizeC * log(feed.obs$PreySize) + coefs$logTempC * log(feed.obs$MeanTemp))
			# /24 to convert hours to day scale units
		est <- est/24
	}
	return(est)
}

# per capita attack rate eqn from Novak and Wootton (2008) Ecology
AttackRateEqn<-function(F_i, F_x, A_x, N_i, h_i){ 
	att <- (F_i * A_x)/((F_x - A_x) * N_i * h_i)
	att[is.na(att)]<-0
	att[is.infinite(att)]<-0
	return(att)
}

# wrapper for attack rate calculation
AttackRateWrap<-function(feed.draw,quad.draw,feed){
	# Attack rate calculation
	PreySpNF <- sort(unique(feed$Prey))
	PreySp <- PreySpNF[which(PreySpNF!='Not Feeding')]

	# Fraction of *all* individuals feeding on prey i
	A<-prop.table(table(factor(feed.draw$Prey,levels=PreySpNF)))
	A<-A[which(names(A)!='Not Feeding')]
	
	# Fraction of *feeding* individuals feeding on prey i
	F<-prop.table(table(factor(feed.draw$Prey,levels=PreySp), exclude='Not Feeding'))
	
	# Identify species 'x', used throughout calculation, as species fed on most frequently.  (The "[1]" ensures that only a single entry is return when there is a tie.)
	Sp_x<-which(F==max(F))[1]

	# Mean handling times for survey
	h<-xtabs(x~Prey,aggregate(feed.draw$HTime, by=list(Prey=factor(feed.draw$Prey,levels=PreySp)), mean,na.rm=TRUE))
	
	# Mean prey abundances
	N<-apply(quad.draw,2,mean,na.rm=TRUE)
	
	AttR<-AttackRateEqn(F, F[Sp_x], A[Sp_x], N, h)
	return(AttR)
}

# random draw from feeding surveys (w/ replacement)
Feed.draw<-function(feed){
	# Random draw of surveys
	surv.samp<-sort(sample(unique(feed$CensusID),replace=TRUE))
	# Pull data from drawn surveys
	feed.draw<-do.call(rbind, lapply(as.list(surv.samp),function(x){subset(feed,CensusID==x)}))
	return(feed.draw)
}
	
# random draw from quadrats (w/ replacement)
Quad.draw<-function(quads){
	quad.samp<-quads[sample(1:nrow(quads),replace=TRUE),]
}

# Wrapper for bootstrapping (B=number of bootstrap replicates)
Boot.AttackRate<-function(B,feed,quads){
	out<-replicate(B,AttackRateWrap(Feed.draw(feed),Quad.draw(quads),feed))
	return(out)
}

nonpara_bs = function(data){
	# data setup
	feed<-data$feed
	temps<-data$temps
	quads<-data$quads

	htimes<-data$htimes_reg
	htimes<-subset(htimes,ConLevel==0.8 & Pred=='Haustrum scobina')[,c(1:3,6:10)]
	htimes<-subset(htimes,Type=='Weighted'|Type=='Mean')
	htimes$Prey[htimes$Prey == "Notoacmea Radialspokes"] = "Notoacmea Radial"

	feed$HTime<-NA
	for(p in unique(feed$Prey)){
		feed[which(feed$Prey==p),]$HTime<-EstHTime(htimes,subset(feed,Prey==p),p)
	}
	
	AttR<-t(Boot.AttackRate(B=nsamples,feed,quads))

	return(AttR)
}
