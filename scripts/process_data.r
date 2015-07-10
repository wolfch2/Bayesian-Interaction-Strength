process_htimes_data = function(htimes){
	htimes = htimes[which(htimes$Pred == "Haustrum scobina"),]
	htimes = htimes[which(htimes$est._._consumed > 80),]

	# combine data for species w/ similar handling times
	htimes[which(htimes$Prey == "Chamaesipho brunnea"),]$Prey = "Chamaesipho columna"

	prey_types = unique(htimes$Prey)[-c(3,9)]
	prey_types = c(prey_types, "Austrolittorina antipodum")
	return(list(htimes=htimes, prey_types=sort(prey_types)))
}

process_survey_data = function(feed, temps, prey_types){
	# used Mark's code here
	colnames(feed)[c(7,9,10,11,13)]<- c('Pred','EstPredSize','ChaemPred','Prey','EstPreySize') 
	temps<-subset(temps,Site=='TaurangaHead')
	J2005<-subset(temps,Monthly=='2005/07')
	J2005$Monthly<-'2005/06'
	temps<-rbind(J2005,temps)
	feed<-subset(feed,Site=='TaurangaHead')
	# Assign each feeding observation its monthly temperature
	feed$Date<-as.Date(feed$Date, "%m/%d/%Y")
	feed$Monthly<-format(as.POSIXct(feed$Date), "%Y/%m")
	feed<-merge(feed,temps,all.x=TRUE)

	feed$Prey[which(feed$Prey == "Notoacmea Radialspokes")] = "Notoacmea Radial"
	
	feeding_status = table(feed$Prey)[prey_types]
	total_surveyed = length(feed$Prey)
	return(list(feed=feed, feeding_status=feeding_status, total_surveyed=total_surveyed))
}

# the comments have escaped! (seriously)
process_abundances_data = function(quads, prey_types){
	quads$Species[which(quads$Species == "Notoacmea Radialspokes")] = "Notoacmea Radial"
	quads<-subset(quads,TideZone=='High' & Site=='TaurangaHead')
	quads<-subset(quads,Transect!=3)
	quads<-subset(quads,Species %in% prey_types)
	quads$SurveyID<-paste(quads$Year, quads$Season, quads$Transect, quads$Quad,sep='.')
	quads2<-reshape(quads,v.names='Density',idvar='SurveyID',timevar='Species', direction='wide',drop=1:7)
	row.names(quads2)<-quads2[,1]
	quads2<-quads2[,-1]
	colnames(quads2)<-sub('Density.','',colnames(quads2))
	quads<-quads2[,order(colnames(quads2))]
	return(list(quads=quads))
}

# setup data for attack rate estimation
process_data = function(raw_data){
	data =  process_htimes_data(raw_data$htimes) # do this first to get prey types
	data = c(data, process_survey_data(raw_data$feed, raw_data$temps, data$prey_types),
		process_abundances_data(raw_data$quads, data$prey_types), list(htimes_reg=raw_data$htimes_reg))
	return(data)	
}

