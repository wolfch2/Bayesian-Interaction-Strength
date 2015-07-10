################################################# median 2d -- compare post. median to sample prop. w/ the number we observed not feeding

x_0 = 1629 # number not feeding (from our dataset)
x_i = seq(1,50,by=1) # numbers feeding

# here, the post. dist. of the ratio alpha_i/alpha_0 is F-dist, so we can use qf to get median
c_seq = c(1,1/2,1/3,1/9,0) # values of sparsity parameter to consider
med_mat = matrix(NA,nrow=length(x_i),ncol=length(c_seq),)
# compute difference in logs (log ratio) between sample prop. and post. median
for(i in 1:ncol(med_mat))
	med_mat[,i] = log((x_i + c_seq[i]) / (x_0 + c_seq[i]) * qf(.5,2*(x_i+c_seq[i]),2*(x_0+c_seq[i])),base=10) - log(x_i/x_0,base=10)

colnames(med_mat) = as.character(fractions(c_seq))
med_mat_df = data.frame(x = x_i, med_mat)

df = melt(med_mat_df, id.vars=list("x"))
colnames(df)[2] = "c"
df$c = factor(df$c, levels=levels(df$c), labels=round(c_seq,3)) # note df$c levels and c_seq are aligned (careful)

# http://astrostatistics.psu.edu/su07/R/html/grDevices/html/plotmath.html
# http://vis.supstat.com/2013/04/mathematical-annotation-in-r/
# http://stackoverflow.com/questions/11344561/controlling-line-color-and-line-type-in-ggplot-legend
# http://stackoverflow.com/questions/7034647/save-ggplot-within-a-function
dev.new(width=6, height=6)
print(ggplot(df, aes(x=x,y=value,linetype=c,color=c)) +
	scale_linetype_manual(values = c("dashed","dashed","solid","dotted","dotted")) +
	scale_color_manual(values = c("gray","#666666","black","#666666","gray")) +
	guides(linetype=guide_legend(title="Sparsity (c)"),
	color=guide_legend(title="Sparsity (c)")) +
	geom_line(size=1) +
	theme_bw() +
	coord_cartesian(xlim=c(0,50),ylim=c(-.25,.25)) +
	xlab("Predators Feeding") +
	ylab(expression(paste(	"log"[10], bgroup("(",frac(hat(theta)[Bayes],"X"[i]/"X"[0]),")")))) +
	theme(legend.position = c(1, 1), legend.justification = c(1, 1),
		axis.title.y = element_text(size = 10, face = "italic", angle = 0)))
savePlot("output/med_c.pdf", "pdf")
dev.off()

################################################# median level plots

# now we look at c=1/3 only and consider multiple numbers observed not feeding
n = 100
N = 2089
x_0 = x_1 = 10^(seq(0, log(2089,base=10),length=n)) # linearly spaced on log scale
c = 1/3

gr = expand.grid(x_0,x_1)
x_0 = gr[,1]
x_1 = gr[,2]
z = rep(NA, length(x_0))

for(i in 1:nrow(gr)){
	est = (x_1[i] + c) / (x_0[i] + c) * qf(.5,2*(x_1[i]+c),2*(x_0[i]+c))
	z[i] = log(est) - log(x_1[i] / x_0[i])
}

df = data.frame(x=x_0,y=x_1,z=z)

pts_df = data.frame(cbind(1629, c(265, 185, 3, 2, 2, 1, 1, 1))) # observed survey results

p1 = ggplot(df, aes(x=x,y=y,z=z)) +
	stat_summary2d(bins=100, fun = mean) +
	scale_x_log10(breaks=c(10,50,500,1000),expand=c(0,0)) +
	scale_y_log10(breaks=c(10,50,500,1000),expand=c(0,0)) +
	theme_bw() +
	annotation_logticks() +
	scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlGn")),
		name=expression(paste(	"log"[10], bgroup("(",frac(hat(theta)[Bayes],"X"[i]/"X"[0]),")")))) +
	geom_point(data=pts_df, aes(x=X1,y=X2,z=NULL)) +
	xlab(expression(paste("Predators not feeding (","X"[0], ")"))) +
	ylab(expression(paste("Predators feeding (","X"[i], ")"))) +
	coord_equal()

################################################# median opt c

# like above plot, but now we find the value of c to minimize
# difference between post med. and sample prop.	

# note this is zero for all c when x_0 = x_i
oppy = function(c, x_0, x_i){
	mle = x_i/x_0
	post_med = (x_i + c) / (x_0 + c) * qf(.5,2*(x_i+c),2*(x_0+c))
	post_med_err = log(post_med) - log(mle)
	return(max(abs(post_med_err)))
}

for(i in 1:nrow(gr)){
	if(abs(x_0[i] - x_1[i]) > .0001)
		z[i] = optimize(oppy,c(.3,4),x_0=x_0[i],x_i=x_1[i])$minimum
}

df = data.frame(x=x_0,y=x_1,z=z)

p2 = ggplot(df, aes(x=x,y=y,z=z)) +
	stat_summary2d(bins=100, fun = mean) +
	scale_x_log10(breaks=c(10,50,500,1000),expand=c(0,0)) +
	scale_y_log10(breaks=c(10,50,500,1000),expand=c(0,0)) +
	theme_bw() +
	annotation_logticks() +
	scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlGn")),
		name=expression("c"[opt]),
		limits=c(.3,.334)) +
	geom_point(data=pts_df, aes(x=X1,y=X2,z=NULL)) +
	xlab(expression(paste("Predators not feeding (","X"[0], ")"))) +
	ylab(expression(paste("Predators feeding (","X"[i], ")"))) +
	coord_equal()	

dev.new(width=12,height=6)
grid.arrange(p1, p2, nrow=1)
savePlot("output/c_vs_med.pdf", "pdf")
dev.off()

