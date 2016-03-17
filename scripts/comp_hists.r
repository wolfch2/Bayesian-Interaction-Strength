# histograms for attack rate estimator components posterior dists.

comp_list = readRDS("output/comp_list.RDS") # load components (saved c=1/3 Bayes results here)
colnames(comp_list[[3]]) = colnames(comp_list[[1]])
comp_list$AR = comp_list$ratio_mat / comp_list$mean.h.times / comp_list$N_mat

log_comp_list = sapply(comp_list, log, base=10, simplify=FALSE)
log_comp_list = sapply(log_comp_list, melt, simplify=FALSE)

# restrict plots to 95% intervals
log_comp_list = sapply(log_comp_list, function(x){
	cuts = quantile(x$value, c(.025,.975))
	out = x[(x$value > cuts[1] & x$value < cuts[2]) | x[,2]=="Epopella plicata",]
	if(nrow(out) > 500)
		return(out)
	return(x) # epo h-time
	}, simplify=FALSE)

df = ldply(log_comp_list)[,c(1,3,4)] # convert to data frame for plotting
colnames(df) = c("Component", "Species", "value")

# sort species by median AR
df$Species = factor(df$Species, levels=names(sort(tapply(df[df$Component=="AR",]$value,
	df[df$Component=="AR",]$Species, median), decreasing=TRUE)))
df$Species = factor(df$Species, labels=gsub(" ", "\n", levels(df$Species)))
df$Component = factor(df$Component, levels=c("ratio_mat","mean.h.times","N_mat","AR"),
	labels=c("alpha","h","n","AR"))

# plot labels
lab_fun = function(variable, value){
	if(variable == "Species")
		return(levels(df$Species)[value])
	labels=list(bquote(atop("Feeding Ratio",alpha[i]/alpha[0])),
	bquote(atop("H. Time",eta[i])),
	bquote(atop("Abundance",nu[i])),
	bquote(atop("Attack Rate",xi[i])))
	return(labels)
}

# http://stackoverflow.com/questions/4027004/restrict-y-axis-range-on-ggplotgeom-density
# https://github.com/hadley/ggplot2/wiki/Plotmath
# http://astrostatistics.psu.edu/su07/R/html/grDevices/html/plotmath.html
# reversing log ticks appears hard:
# https://groups.google.com/forum/#!topic/ggplot2/OlGA8Gm9O7w
# instead, change to lightish gray
dev.new(width=12.2, height=12)
print(ggplot(df, aes(x=value)) +
	geom_histogram(aes(y = ..density..), fill="#888888", binwidth=.1, position="identity") +
	facet_grid(Species ~ Component, scales="free_x", space="free_x", labeller=lab_fun) +
	theme_bw() +
	coord_cartesian(ylim=c(0, 6), expand=FALSE) +
	ylab("Probability") +
	scale_y_continuous(breaks=c(0,3,6)) +
	theme(strip.text.y = element_text(size = 10, face = "italic", angle = 0),
		panel.margin = unit(1, "lines"),
		legend.position = "right",
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.title.x = element_blank()) +
	scale_x_continuous(breaks=-10:10,
		labels=trans_format("I", math_format(10^.x))) +
	annotation_logticks(sides="b",short=unit(0.05,"cm"),
		mid=unit(0.1,"cm"),long=unit(0.15, "cm"),size=.25))

savePlot(filename = "output/hist_all.pdf", type = "pdf")
dev.off()

