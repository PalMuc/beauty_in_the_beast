library(paletteer)
library(ggplot2)
distfile = "all_distances_levels_ratios_within_Uniplacotomia_lines_no-clades.tsv"
distdata = read.table(distfile, header=TRUE, sep="\t")
boxplotfile = "all_distances_levels_ratios_between_Uniplacotomia_lines_boxplots_no-clades.tsv"
boxplotdist = read.table(boxplotfile, header=TRUE, sep="\t")

levels = c("class","order","family","genus")
genes = c("ATP6", "COX1", "COX2", "COX3", "CYTB", "NAD1", "NAD2", "NAD3", "NAD4", "NAD4L", "NAD5", "NAD6", "concat")

#pdf(file="plots/all_distances_levels_ratios_within_Uniplacotomia_lines_no-clades.pdf", width=8, height=6)
svg(file="plots/all_distances_levels_ratios_within_Uniplacotomia_lines_no-clades.svg", width=8, height=6)

palette = paletteer_d("dichromat::Categorical_12")

par(mar=c(4,5,1,1))
boxplot(avgdistance ~ taxlevel, data = boxplotdist, ylim=c(0,0.3), ylab="Mean uncorrected distance", xlab=NULL, at=c(1,2,3,4), outline=FALSE, boxwex=0.5, frame.plot=TRUE, axes=FALSE, cex.lab=1.5)
axis(1,at=c(1,2,3,4),labels=levels, par(mgp=c(4,1,0)), tick=TRUE, cex.axis=1.3)
axis(2,at=seq(0,0.3,0.1), labels=seq(0.0,0.3,0.1), mgp=c(1,1,0), cex.axis=1.3)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(avgdistance ~ taxlevel, data = boxplotdist, ylim=c(0,0.3), ylab=NULL, xlab=NULL, at=c(1,2,3,4), outline=FALSE, boxwex=0.5, frame.plot=TRUE, axes=FALSE, cex.lab=1.5)
lines(distdata[,2], distdata[,3], type = "p", lwd=1.5, col=palette[1])
lines(distdata[,2], distdata[,4], type = "p", lwd=1.5, col=palette[2])
lines(distdata[,2], distdata[,5], type = "p", lwd=1.5, col=palette[3])
lines(distdata[,2], distdata[,6], type = "p", lwd=1.5, col=palette[4])
lines(distdata[,2], distdata[,7], type = "p", lwd=1.5, col=palette[5])
lines(distdata[,2], distdata[,8], type = "p", lwd=1.5, col=palette[6])
lines(distdata[,2], distdata[,9], type = "p", lwd=1.5, col=palette[7])
lines(distdata[,2], distdata[,10], type = "p", lwd=1.5, col=palette[8])
lines(distdata[,2], distdata[,11], type = "p", lwd=1.5, col=palette[9])
lines(distdata[,2], distdata[,12], type = "p", lwd=1.5, col=palette[10])
lines(distdata[,2], distdata[,13], type = "p", lwd=1.5, col=palette[11])
lines(distdata[,2], distdata[,14], type = "p", lwd=1.5, col=palette[12])
lines(distdata[,2], distdata[,15], type = "p", lwd=1.5, col="white")

legend(4,0.31, legend=genes, fill=palette, cex=1, bty='n')
legend(4,0.155, legend=NA, fill="white", cex=1, bty='n')

dev.off()





