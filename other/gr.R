require(coda)

setwd("C:/Users/Chris/Dropbox/stat phd/project_code - convergence")

files = list.files("data")

gr = do.call("rbind", lapply(files, function(file){
	print(file)
	data = do.call(rbind, readRDS(paste0("data/", file)))
	out = data.frame(gelman.diag(data, multivariate=FALSE, transform=TRUE, autoburnin=FALSE)$psrf)
	out = data.frame(group=substr(file,1,nchar(file)-4), param=rownames(out), out)
}))

write.csv(gr[order(gr[,4],decreasing=TRUE),], "gr_stat.csv", row.names=FALSE)