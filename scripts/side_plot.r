# side-by-side interval and point estimate plots for all methods

prey_types = colnames(results[[1]])
results_quants = lapply(results, apply, 2, quantile, probs = c(.025, .5, .975)) # get median, quantiles
	
# http://stackoverflow.com/questions/13847936/in-r-plotting-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot
results_log_quants = lapply(lapply(results_quants,function(x)x+10^-7),log,base=10) # log transform
# have to include species names in the data (not just as row names)
results_log_quants = sapply(results_log_quants, function(x) data.frame(t(x)), simplify=FALSE)
results_log_quants = sapply(results_log_quants, function(x) data.frame(rownames(x),x), simplify=FALSE)
# variable for species names needs a name (for ldply)
df = ldply(results_log_quants)
colnames(df) = c("Method","Species","lower","mid","upper")
df$Species = gsub(" ", "\n", df$Species)
df$Species = factor(df$Species, levels=rev(unique(df[order(df$mid),"Species"])))

# fix method names for plotting
df$Method = factor(df$Method,
	levels=rev(c("Nonpara BS","Para BS","Bayes (c=0)",
	"Bayes (c=0.333)","Bayes (c=1)")),
	labels=rev(c("Nonparametric Bootstrap ","Parametric Bootstrap ",
	"Bayes (c=0)","Bayes (c=0.333)","Bayes (c=1)")))

# http://stackoverflow.com/questions/11214012/set-only-lower-bound-of-a-limit-for-ggplot
dev.new(width=10, height=12)
print(ggplot(df, aes(x=10^mid,y=as.numeric(Method),color=Method)) +
	geom_errorbar(aes(ymin=as.numeric(Method)-.4, ymax=as.numeric(Method)+.4),
		width=0, size=1) +
	facet_grid(Species ~ .) +
	geom_errorbarh(aes(xmin=10^lower, xmax=10^upper), height=0, size=1) +
	theme_bw() +
	scale_x_log10(breaks=10^(-7:-3),
		labels=trans_format("log10", math_format(10^.x)),limits=c(10^-7,10^-3)) +
	scale_y_continuous(limits=c(.25,NA)) +
	annotation_logticks(sides="b",short=unit(0.05,"cm"),
		mid=unit(0.1,"cm"),long=unit(0.15, "cm"),size=.75) +
	theme(strip.text.y = element_text(size = 10, face = "italic", angle = 0),
		axis.title.y = element_blank(), axis.ticks.y = element_blank(),
		axis.text.y = element_blank(), legend.position = "bottom",
		panel.margin.y = unit(.75, "lines")) +
	scale_color_manual(values=brewer.pal(length(levels(df$Method)),"Set1"),
		limits = rev(levels(df$Method)),
		guide = guide_legend(nrow=2,byrow=TRUE,title=NULL)) +
	xlab(bquote("Attack Rate (" ~ xi[i] ~ ")")))

savePlot(filename = "output/side.pdf", type = "pdf")
dev.off()

