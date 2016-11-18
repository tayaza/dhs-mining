xyplot(V1~V2, data = my_data, pch = 20, cex = 1.5, groups = V3,
	auto.key = list(columns = 2),
	xlab = "Distance from DHS center",
	ylab = "SNPs",
	scales = list(tck = c(1,0), x=list(tick.number = 9)),
	#Remove border
	#par.settings = list(superpose.symbol = list(pch = 20, cex = 1.5)),
	par.settings = list(axis.line = list(col = NA), superpose.symbol = list(pch = 20, cex = 1.5)),

	panel = function(x,y,...){
		panel.xyplot(x,y,...);
		panel.axis(side = c("bottom"), outside = TRUE, ticks = TRUE)
}
)