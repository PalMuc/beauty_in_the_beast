infile = "../data/combined-ALL_protein_NAD3.tsv"

levdata = read.table(infile, header=TRUE, sep="\t")

levels = c("genus","clade")

phylanames = unique(levdata[,1])
phylaindex = match(levdata[,1],phylanames)
distvalues = levdata[,3]
barposition = match(levdata[,2],levels) * 3 - 3 + match(levdata[,1],phylanames)
datagroups = split(levdata[,3], barposition)

svg(file="../plots/Figure_S2B.svg", width=8, height=6)

par(mar=c(4,5,1,1))


boxplot(datagroups, ylim=c(0,0.4), ylab="Mean uncorrected distance", at=c(1,2,4,5), xlab=NULL, boxwex=0.7, frame.plot=TRUE, axes=FALSE, col="white", cex.lab=1.7)
axis(1,at=c(1.5,4.5),labels=NA, tick=FALSE)
axis(2,at=seq(0,0.4,0.1), labels=seq(0.0,0.4,0.1), par(mgp=c(2,1,0)), cex.axis=1.3)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(datagroups, ylim=c(0,0.4), ylab=NULL, at=c(1,2,4,5), xlab=NULL, boxwex=0.7, frame.plot=TRUE, axes=FALSE, col="white", cex.lab=1.7)
stripchart(distance ~ points, data = levdata, xlab=NULL, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.35, at=c(1,2,4,5), pch = 19, add = TRUE)
text(1,0.05,labels="within\ngenera", col="black", cex=1.3)
text(2,0.35,labels="between\ngenera", col="black", cex=1.3)
text(4,0.05,labels="within\nclades", col="black", cex=1.3)
text(5,0.35,labels="between\nclades", col="black", cex=1.3)


dev.off()





