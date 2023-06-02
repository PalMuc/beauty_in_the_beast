## INTER DATA:
distclassinter = "inter-CLASS_protein_ALL.tsv"
distorderinter = "inter-ORDER_protein_ALL.tsv"
distfamilyinter = "inter-FAMILY_protein_ALL.tsv"
distgenusinter = "inter-GENUS_protein_ALL.tsv"
distcladeinter = "inter-CLADE_protein_ALL.tsv"

distorderinterPTH = "inter-ORDER_protein_P-TH.tsv"
distfamilyinterPTH = "inter-FAMILY_protein_P-TH.tsv"
distgenusinterPTH = "inter-GENUS_protein_P-TH.tsv"
distcladeinterPTH = "inter-CLADE_protein_P-TH.tsv"

distorderinterTH = "inter-ORDER_protein_T-H.tsv"
distfamilyinterTH = "inter-FAMILY_protein_T-H.tsv"
distgenusinterTH = "inter-GENUS_protein_T-H.tsv"
distcladeinterTH = "inter-CLADE_protein_T-H.tsv"

distclassinterdata = read.table(distclassinter, header=TRUE, sep="\t")
distorderinterdata = read.table(distorderinter, header=TRUE, sep="\t")
distfamilyinterdata = read.table(distfamilyinter, header=TRUE, sep="\t")
distgenusinterdata = read.table(distgenusinter, header=TRUE, sep="\t")
distcladeinterdata = read.table(distcladeinter, header=TRUE, sep="\t")

distorderinterPTHdata = read.table(distorderinterPTH, header=TRUE, sep="\t")
distfamilyinterPTHdata = read.table(distfamilyinterPTH, header=TRUE, sep="\t")
distgenusinterPTHdata = read.table(distgenusinterPTH, header=TRUE, sep="\t")
distcladeinterPTHdata = read.table(distcladeinterPTH, header=TRUE, sep="\t")

distorderinterTHdata = read.table(distorderinterTH, header=TRUE, sep="\t")
distfamilyinterTHdata = read.table(distfamilyinterTH, header=TRUE, sep="\t")
distgenusinterTHdata = read.table(distgenusinterTH, header=TRUE, sep="\t")
distcladeinterTHdata = read.table(distcladeinterTH, header=TRUE, sep="\t")


## INTRA DATA:
distclassintra = "intra-CLASS_protein_ALL.tsv"
distorderintra = "intra-ORDER_protein_ALL.tsv"
distfamilyintra = "intra-FAMILY_protein_ALL.tsv"
distgenusintra = "intra-GENUS_protein_ALL.tsv"
distcladeintra = "intra-CLADE_protein_ALL.tsv"


distclassintradata = read.table(distclassintra, header=TRUE, sep="\t")
distorderintradata = read.table(distorderintra, header=TRUE, sep="\t")
distfamilyintradata = read.table(distfamilyintra, header=TRUE, sep="\t")
distgenusintradata = read.table(distgenusintra, header=TRUE, sep="\t")
distcladeintradata = read.table(distcladeintra, header=TRUE, sep="\t")


groups = c("class","order","family","genus","clade")
colorset = c("#00BFFF", "#ff5555ff", "#32CD32", "#FFA500", "#FF69B4")


#pdf(file="plots/protein_p-distances_graph_combined_color.pdf", width=16, height=24)
svg(file="plots/protein_p-distances_graph_combined_color.svg", width=16, height=24)

par(mar=c(6,6,1,1))
par(mfrow = c(5, 2))  ## set the layout to be 5 by 2


## inter-CLASS plot:
boxplot(distance ~ gene, data = distclassinterdata, ylab="Mean uncorrected distance", xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distclassinterdata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend="between Poly- and Uniplacotoma classes", col=colorset[1], pch=19, cex=2, bty='n')
text(1.5,0.58,labels="A1", col="black", cex=2)
stripchart(distance ~ gene, data = distclassinterdata, vertical = TRUE, offset = TRUE, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 19, add = TRUE, col=colorset[1], cex=1.4)

## intra-CLASS plot:
boxplot(distance ~ gene, data = distclassintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distclassintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend="within Uniplacotoma", col=colorset[1], pch=0, cex=2, bty='n')
text(1.5,0.58,labels="A2", col="black", cex=2)
stripchart(distance ~ gene, data = distclassintradata, vertical = TRUE, offset = TRUE, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 0, add = TRUE, col=colorset[1], cex=1.4, lwd=1.3)

## inter-ORDER plot:
boxplot(distance ~ gene, data = distorderinterdata, ylab="Mean uncorrected distance", xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distorderinterdata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend=c("between Poly- and Uniplacotoma orders","between Uniplacotoma orders"), col=colorset[2], pch=c(19,23), cex=2, bty='n')
text(1.5,0.58,labels="B1", col="black", cex=2)
stripchart(distance ~ gene, data = distorderinterPTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 19, add = TRUE, col=colorset[2], cex=1.4)
stripchart(distance ~ gene, data = distorderinterTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 23, add = TRUE, col=colorset[2], cex=1.4, lwd=1.3)

## intra-ORDER plot:
boxplot(distance ~ gene, data = distorderintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distorderintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend="within Uniplacotoma orders", col=colorset[2], pch=0, cex=2, bty='n')
text(1.5,0.58,labels="B2", col="black", cex=2)
stripchart(distance ~ gene, data = distorderintradata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 0, add = TRUE, col=colorset[2], cex=1.4, lwd=1.3)

## inter-FAMILY plot:
boxplot(distance ~ gene, data = distfamilyinterdata, ylab="Mean uncorrected distance", xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distfamilyinterdata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend=c("between Poly- and Uniplacotoma families","between Uniplacotoma families"), col=colorset[3], pch=c(19,23), cex=2, bty='n')
text(1.5,0.58,labels="C1", col="black", cex=2)
stripchart(distance ~ gene, data = distfamilyinterPTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 19, add = TRUE, col=colorset[3], cex=1.4)
stripchart(distance ~ gene, data = distfamilyinterTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 23, add = TRUE, col=colorset[3], cex=1.4, lwd=1.3)

## intra-FAMILY plot:
boxplot(distance ~ gene, data = distfamilyintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distfamilyintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
text(1.5,0.58,labels="C2", col="black", cex=2)
legend('topright', legend="within Uniplacotoma families", col=colorset[3], pch=0, cex=2, bty='n')
stripchart(distance ~ gene, data = distfamilyintradata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 0, add = TRUE, col=colorset[3], cex=1.4, lwd=1.3)

## inter-GENUS plot:
boxplot(distance ~ gene, data = distgenusinterdata, ylab="Mean uncorrected distance", xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distgenusinterdata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
text(1.5,0.58,labels="D1", col="black", cex=2)
legend('topright', legend=c("between Poly- and Uniplacotoma genera","between Uniplacotoma genera"), col=colorset[4], pch=c(19,23), cex=2, bty='n')
stripchart(distance ~ gene, data = distgenusinterPTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 19, add = TRUE, col=colorset[4], cex=1.4)
stripchart(distance ~ gene, data = distgenusinterTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 23, add = TRUE, col=colorset[4], cex=1.4, lwd=1.3)

## intra-GENUS plot:
boxplot(distance ~ gene, data = distgenusintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distgenusintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
text(1.5,0.58,labels="D2", col="black", cex=2)
legend('topright', legend="within Uniplacotoma genera", col=colorset[4], pch=0, cex=2, bty='n')
stripchart(distance ~ gene, data = distgenusintradata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 0, add = TRUE, col=colorset[4], cex=1.4, lwd=1.3)

## inter-CLADE plot:
boxplot(distance ~ gene, data = distcladeinterdata, ylab="Mean uncorrected distance", xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distcladeinterdata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend=c("between Poly- and Uniplacotoma clades","between Uniplacotoma clades"), col=colorset[2], pch=c(19,23), cex=2, bty='n')
text(1.5,0.58,labels="E1", col="black", cex=2)
stripchart(distance ~ gene, data = distcladeinterPTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 19, add = TRUE, col=colorset[5], cex=1.4)
stripchart(distance ~ gene, data = distcladeinterTHdata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 23, add = TRUE, col=colorset[5], cex=1.4, lwd=1.3)

## intra-CLADE plot:
boxplot(distance ~ gene, data = distcladeintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
axis(1,at=c(2,4,6,8,10,12,14,16,18,20,22,24,26),labels=c("ATP6","COX1","COX2","COX3","CYTB","NAD1","NAD2","NAD3","NAD4","NAD4L","NAD5","NAD6","concat"), mgp=c(4,1,0), tick=TRUE, las=2, cex.axis=1.6)
axis(2,at=seq(0,0.5,0.1), labels=seq(0.0,0.5,0.1), mgp=c(1,1,0), cex.axis=1.6)
rect(par("usr")[1], par("usr")[3],par("usr")[2], par("usr")[4],col="#D3D3D3")
par(new = TRUE)
boxplot(distance ~ gene, data = distcladeintradata, ylab=NULL, xlab=NULL, ylim=c(0,0.6), at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), boxwex=1.5 , col="white", frame.plot=TRUE, axes=FALSE, cex.lab=2, outline=FALSE)
legend('topright', legend="within Uniplacotoma clades", col=colorset[5], pch=0, cex=2, bty='n')
text(1.5,0.58,labels="E2", col="black", cex=2)
stripchart(distance ~ gene, data = distcladeintradata, vertical = TRUE, offset = TRUE, method = "jitter", jitter = 0.5, at=c(2,4,6,8,10,12,14,16,18,20,22,24,26), pch = 0, add = TRUE, col=colorset[5], cex=1.4, lwd=1.3)


dev.off()






