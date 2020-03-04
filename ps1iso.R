library(BiocManager)
install()
library(limma)
library(edgeR)
library(Glimma)
#library(AnnotationDbi)
options(java.parameters = "-Xmx8000m")
library(XLConnect)
library(RColorBrewer)
library(gplots)

cts <- read.delim(file = "quant.featurecounts.counts.unstr.txt", sep = "\t", header = T, stringsAsFactors = F)
head(cts)
dim(cts) # 52459    28

cts <- cts[, c(1,6,7,9,10,12,13,15,16,18,19,21,22)]
head(cts)
mode(cts) # "list"

colnames(cts) <- c("symbol", "ps1_1", "ps1iso_1", "wt_1", "wtiso_1",
                   "ps1_2", "ps1iso_2", "wt_2", "wtiso_2",
                   "ps1_3", "ps1iso_3", "wt_3", "wtiso_3")
head(cts)

cts <- cts[, c(1,2,6,10,3,7,11,4,8,12,5,9,13)]
head(cts)

sum(duplicated(cts$symbol)) # 0
rownames(cts) <- cts$symbol
cts <- cts[, -1]
head(cts)
tail(cts)
dim(cts) # 52459    12

cts[rownames(cts) %in% c("Tfe3", "Tfeb", "App", "Ctsb", "Ctsd", "Acp2", "Hexa", "Psen1", "Psen2", "Mapt", "Akt1"), ]

## Now before making DGElist, make group, batch etc. variables...
group <- as.factor(rep(c("ps1", "ps1iso", "wt", "wtiso"), each = 3))
group

lane <- as.factor(rep(c("L003", "L004", "L005"), 4))
lane

## Make DGEList
y.dge <- DGEList(cts, group = group, genes = rownames(cts))
y.dge
dim(y.dge) # 52459    12

samplenames <- colnames(y.dge)
samplenames

# Adding confounding variables to dge object
y.dge$samples$lane <- lane

y.dge$samples

## Transformation from the raw-scale
cpm <- cpm(y.dge)
head(cpm)
tail(cpm)

lcpm <- cpm(y.dge, log = TRUE) # Default prior.count = 2. This means it adds 2/L to the counts before log2 transformation, where L is mean library size.

## Removing genes that are lowly expressed
table(rowSums(y.dge$counts==0)==12)
# FALSE  TRUE
# 30588 21871

## The smallest lib size is 39.5 million ~ 40 mn. Therefore, cutoff = 10/40 = 0.25
## Or alternatively, the cutoff as discussed in edgeR,
L <- mean(y.dge$samples$lib.size) * 1e-6
M <- median(y.dge$samples$lib.size) * 1e-6
c(L,M) # 42.07802 41.97634
lcpm3.cutoff <- log2(10/M + 2/L) # -1.807123
2^-1.807123 # 0.2857602 (~0.29)

## Let's go with the cut-off of 0.25
keep.exprs_.25 <- rowSums(cpm >= 0.25) >= 3
table(keep.exprs_.25)
# FALSE  TRUE
# 36952 15507
### For cut-off 0.3
keep.exprs_.3 <- rowSums(cpm >= 0.3) >= 3
table(keep.exprs_.3)
# FALSE  TRUE
# 37379 15080
###
keep.exprs_ER <- filterByExpr(y.dge, group = group)
table(keep.exprs_ER)
# FALSE  TRUE
# 36801 15658


## Low-count filtration (threshold cpm = 0.25)
y.filt <- y.dge[keep.exprs_.25,, keep.lib.sizes=FALSE]
y.filt
dim(y.filt) # 15507    12
###
y.filt2 <- y.dge[keep.exprs_.3,, keep.lib.sizes=FALSE]
y.filt2
dim(y.filt2)
###
y.filt3 <- y.dge[keep.exprs_ER,, keep.lib.sizes=FALSE]
y.filt3
dim(y.filt3) # 15658    12

## Low count filtration effect
library(RColorBrewer)
nsamples <- ncol(y.filt)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="A. Raw data", xlab="Log-cpm")
abline(v=log(0.25), lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)

lcpm.filt <- cpm(y.filt, log=TRUE)
plot(density(lcpm.filt[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="B. Filtered data", xlab="Log-cpm")
abline(v=log(0.25), lty=3)
for (i in 2:nsamples){
  den <- density(lcpm.filt[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)
dev.print(device = pdf, "Rplot_lowCountFiltEffect_ps1iso.pdf", width = 12, height = 6)
###
library(RColorBrewer)
nsamples2 <- ncol(y.filt2)
col <- brewer.pal(nsamples2, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="A. Raw data", xlab="Log-cpm")
abline(v=log(0.3), lty=3)
for (i in 2:nsamples2){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)

lcpm.filt2 <- cpm(y.filt2, log=TRUE)
plot(density(lcpm.filt2[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="B. Filtered data", xlab="Log-cpm")
abline(v=log(0.3), lty=3)
for (i in 2:nsamples2){
  den <- density(lcpm.filt2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)
dev.print(device = pdf, "Rplot_lowCountFiltEffect2_ps1iso.pdf", width = 12, height = 6)
###
library(RColorBrewer)
nsamples3 <- ncol(y.filt3)
col <- brewer.pal(nsamples3, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm3.cutoff, lty=3)
for (i in 2:nsamples3){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)

lcpm.filt3 <- cpm(y.filt3, log=TRUE)
plot(density(lcpm.filt3[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm3.cutoff, lty=3)
for (i in 2:nsamples3){
  den <- density(lcpm.filt3[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)
dev.print(device = pdf, "Rplot_lowCountFiltEffect3_ps1iso.pdf", width = 12, height = 6)


## Normalising gene expression distributions
y.norm <- calcNormFactors(y.filt, method = "TMM")
y.norm$samples
###
y.norm2 <- calcNormFactors(y.filt2, method = "TMM")
y.norm2$samples
###
y.norm3 <- calcNormFactors(y.filt3, method = "TMM")
y.norm3$samples


## visualization of normalization effect
par(mfrow=c(1, 2))
par(cex.axis=0.7)
boxplot(lcpm.filt, las=2, col=col, main="A. Example: Unnormalised data", ylab="Log-cpm")
abline(h=median(lcpm.filt), col="red")

lcpm.norm <- cpm(y.norm, log=T)
boxplot(lcpm.norm, las=2, col=col, main="B. Example: Normalised data", ylab="Log-cpm")
abline(h=median(lcpm.norm), col="red")
dev.print(pdf, "Rplot_boxPlot_unNormVSNormData_ps1iso.pdf", height = 6, width = 12)
###
par(mfrow=c(1, 2))
par(cex.axis=0.7)
boxplot(lcpm.filt2, las=2, col=col, main="A. Example: Unnormalised data", ylab="Log-cpm")
abline(h=median(lcpm.filt2), col="red")

lcpm.norm2 <- cpm(y.norm2, log=T)
boxplot(lcpm.norm2, las=2, col=col, main="B. Example: Normalised data", ylab="Log-cpm")
abline(h=median(lcpm.norm2), col="red")
dev.print(pdf, "Rplot_boxPlot_unNormVSNormData2_ps1iso.pdf", height = 6, width = 12)
###
par(mfrow=c(1, 2))
par(cex.axis=0.7)
boxplot(lcpm.filt3, las=2, col=col, main="A. Example: Unnormalised data", ylab="Log-cpm")
abline(h=median(lcpm.filt3), col="red")

lcpm.norm3 <- cpm(y.norm3, log=T)
boxplot(lcpm.norm3, las=2, col=col, main="B. Example: Normalised data", ylab="Log-cpm")
abline(h=median(lcpm.norm3), col="red")
dev.print(pdf, "Rplot_boxPlot_unNormVSNormData3_ps1iso.pdf", height = 6, width = 12)


## Unsupervised clustering of samples
par(mfrow=c(1,1))
colors <- brewer.pal(4, "Dark2")
plotMDS(lcpm.norm, col=colors[group], main="PCA plot (Unsupervised Clustering)")
dev.print(device = pdf, "Rplot_PCA_ps1iso.pdf", height=6, width=8)
###
par(mfrow=c(1,1))
colors <- brewer.pal(4, "Dark2")
plotMDS(lcpm.norm2, col=colors[group], main="PCA plot (Unsupervised Clustering)")
dev.print(device = pdf, "Rplot_PCA2_ps1iso.pdf", height=6, width=8)
###
par(mfrow=c(1,1))
colors <- brewer.pal(4, "Dark2")
plotMDS(lcpm.norm3, col=colors[group], main="PCA plot (Unsupervised Clustering)")
dev.print(device = pdf, "Rplot_PCA3_ps1iso.pdf", height=6, width=8)


## Remove batch effect in PCA
# Create design matrix
design <- model.matrix(~0+group+lane)
design
colnames(design) <- gsub("group", "", colnames(design))

length(colnames(design)) # 6
qr(design)$rank # 6
#  This means the design is VALID. The rank of QR and length of colnames of design matrix should be equal in order to
#  to ensure that the design matrix is of full rank and to avoid warning that certain coefficient(s) is(are) unestimmable.

design_c <- model.matrix(~0+group)
design_c
colnames(design_c) <- gsub("group", "", colnames(design_c))

par(mfrow=c(1,1))
lcpm.norm_c <- removeBatchEffect(lcpm.norm, design = design_c, batch = lane)
plotMDS(lcpm.norm_c, col=colors[group], main="PCA plot (Unsupervised Clustering)-Lane Corrected")
###
par(mfrow=c(1,1))
lcpm.norm2_c <- removeBatchEffect(lcpm.norm2, design = design_c, batch = lane)
plotMDS(lcpm.norm2_c, col=colors[group], main="PCA plot (Unsupervised Clustering)-Lane Corrected")
###
par(mfrow=c(1,1))
lcpm.norm3_c <- removeBatchEffect(lcpm.norm3, design = design_c, batch = lane)
plotMDS(lcpm.norm3_c, col=colors[group], main="PCA plot (Unsupervised Clustering)-Lane Corrected")


## Just using pch to better represent graph
pch <- c(0,1,15,16)
lcpm.norm_5 <- cpm(y.norm, log=T, prior.count = 5) # FOR VISUALIZATION PURPOSE
plotMDS(lcpm.norm_5, col=colors[group], pch=pch[group], main = "PCA plot (Unsupervised Clustering)") # y.norm could also be used in place of lcpm.norm_5
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2)
dev.print(device = pdf, "Rplot_PCA_ps1iso.pdf", height=6, width=8)
###
pch <- c(0,1,15,16)
lcpm.norm2_5 <- cpm(y.norm2, log=T, prior.count = 5) # FOR VISUALIZATION PURPOSE
plotMDS(lcpm.norm2_5, col=colors[group], pch=pch[group], main = "PCA plot (Unsupervised Clustering)") # y.norm could also be used in place of lcpm.norm_5
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2)
dev.print(device = pdf, "Rplot_PCA2_ps1iso.pdf", height=6, width=8)
###
pch <- c(0,1,15,16)
lcpm.norm3_5 <- cpm(y.norm3, log=T, prior.count = 5) # FOR VISUALIZATION PURPOSE
plotMDS(lcpm.norm3_5, col=colors[group], pch=pch[group], main = "PCA plot (Unsupervised Clustering)") # y.norm could also be used in place of lcpm.norm_5
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2)
dev.print(device = pdf, "Rplot_PCA3_ps1iso.pdf", height=6, width=8)


## After removing "lane" effect
lcpm.norm_5_c <- removeBatchEffect(lcpm.norm_5, design = design_c, batch = lane)
plotMDS(lcpm.norm_5_c, col=colors[group], pch=pch[group], main = "PCA plot (Unsupervised Clustering)-Lane Corrected")
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2)
dev.print(device = pdf, "Rplot_PCA_ps1iso_laneCor.pdf", height=6, width=8)
###
lcpm.norm2_5_c <- removeBatchEffect(lcpm.norm2_5, design = design_c, batch = lane)
plotMDS(lcpm.norm2_5_c, col=colors[group], pch=pch[group], main = "PCA plot (Unsupervised Clustering)-Lane Corrected")
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2)
dev.print(device = pdf, "Rplot_PCA2_ps1iso_laneCor.pdf", height=6, width=8)
###
lcpm.norm3_5_c <- removeBatchEffect(lcpm.norm3_5, design = design_c, batch = lane)
plotMDS(lcpm.norm3_5_c, col=colors[group], pch=pch[group], main = "PCA plot (Unsupervised Clustering)-Lane Corrected")
legend("top", legend=levels(group), pch=pch, col=colors, ncol=2)
dev.print(device = pdf, "Rplot_PCA3_ps1iso_laneCor.pdf", height=6, width=8)


## Plotting Lanes
plotMDS(lcpm.norm_5, labels = lane, col = colors[lane], main = "PCA plot- Sequencing Lanes")
###
plotMDS(lcpm.norm2_5, labels = lane, col = colors[lane], main = "PCA plot- Sequencing Lanes")

## Cluster Dendrograms
dist_5 <- dist(as.matrix(t(lcpm.norm_5)))
#head(dist7_dc)
hclust_5 <- hclust(dist_5)
plot(hclust_5)
dev.print(pdf, "Rplot_clusterDendrogram_ps1iso.pdf", height=6, width=8)

dist_5_c <- dist(as.matrix(t(lcpm.norm_5_c)))
#head(dist7_dc_c)
hclust_5_c <- hclust(dist_5_c)
plot(hclust_5_c, main="Cluster Dendrogram- Lane Corrected")
dev.print(pdf, "Rplot_clusterDendrogram_LaneCor_ps1iso.pdf", height=6, width=8)
###
dist2_5 <- dist(as.matrix(t(lcpm.norm2_5)))
#head(dist7_dc)
hclust2_5 <- hclust(dist2_5)
plot(hclust2_5)
dev.print(pdf, "Rplot_clusterDendrogram2_ps1iso.pdf", height=6, width=8)

dist2_5_c <- dist(as.matrix(t(lcpm.norm2_5_c)))
#head(dist7_dc_c)
hclust2_5_c <- hclust(dist2_5_c)
plot(hclust2_5_c, main="Cluster Dendrogram- Lane Corrected")
dev.print(pdf, "Rplot_clusterDendrogram2_LaneCor_ps1iso.pdf", height=6, width=8)
###
dist3_5 <- dist(as.matrix(t(lcpm.norm3_5)))
hclust3_5 <- hclust(dist3_5)
plot(hclust3_5)
dev.print(pdf, "Rplot_clusterDendrogram3_ps1iso.pdf", height=6, width=8)

dist3_5_c <- dist(as.matrix(t(lcpm.norm3_5_c)))
hclust3_5_c <- hclust(dist3_5_c)
plot(hclust3_5_c, main="Cluster Dendrogram- Lane Corrected")
dev.print(pdf, "Rplot_clusterDendrogram3_LaneCor_ps1iso.pdf", height=6, width=8)


## PCA - NEW
pca <- prcomp(t(lcpm.norm_5))
summary(pca)
pca.propVar <- ((pca$sdev^2)/(sum(pca$sdev^2))) * 100
barplot(pca.propVar, xlab = "Principal Components (PC)", ylab = "Proportion of Variation (%)", ylim = c(0,100), main = "Scree plot - PCA")

pca_c <- prcomp(t(lcpm.norm_5_c))
summary(pca_c)
pca.propVar_c <- ((pca_c$sdev^2)/(sum(pca_c$sdev^2))) * 100
barplot(pca.propVar_c, xlab = "Principal Components (PC)", ylab = "Proportion of Variation (%)", ylim = c(0,100),
        main = "Scree plot - PCA (lane Corrected)")
###
pca2 <- prcomp(t(lcpm.norm2_5))
summary(pca2)
pca.propVar2 <- ((pca2$sdev^2)/(sum(pca2$sdev^2))) * 100
barplot(pca.propVar2, xlab = "Principal Components (PC)", ylab = "Proportion of Variation (%)", ylim = c(0,100), main = "Scree plot - PCA")

pca2_c <- prcomp(t(lcpm.norm2_5_c))
summary(pca2_c)
pca.propVar2_c <- ((pca2_c$sdev^2)/(sum(pca2_c$sdev^2))) * 100
barplot(pca.propVar2_c, xlab = "Principal Components (PC)", ylab = "Proportion of Variation (%)", ylim = c(0,100),
        main = "Scree plot - PCA (lane Corrected)")
###
pca3 <- prcomp(t(lcpm.norm3_5))
summary(pca3)
pca.propVar3 <- ((pca3$sdev^2)/(sum(pca3$sdev^2))) * 100
barplot(pca.propVar3, xlab = "Principal Components (PC)", ylab = "Proportion of Variation (%)", ylim = c(0,100), main = "Scree plot - PCA")

pca3_c <- prcomp(t(lcpm.norm3_5_c))
summary(pca3_c)
pca.propVar3_c <- ((pca3_c$sdev^2)/(sum(pca3_c$sdev^2))) * 100
barplot(pca.propVar3_c, xlab = "Principal Components (PC)", ylab = "Proportion of Variation (%)", ylim = c(0,100),
        main = "Scree plot - PCA (lane Corrected)")

library(scatterplot3d)
par(mar=c(4,4,4,4), mfrow=c(1,1)) # mar = margin
scatterplot3d(pca3$x[, 1:3], angle = -40, pch = pch[group], color = colors[group], grid = F, box = F,
              xlab=paste("PC1, ", round(pca.propVar3[1], 2), "%"),
              ylab=paste("PC2, ", round(pca.propVar3[2], 2), "%"),
              zlab=paste("PC3, ", round(pca.propVar3[3], 2), "%"),
              main = "3D-scatterPlot-PCA")
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca3$x[,1:3], grid = c("xy", "xz", "yz"))
legend("topright", legend=levels(group), pch=pch, col=colors)

par(mar=c(4,4,4,4), mfrow=c(1,1)) # mar = margin
scatterplot3d(pca3_c$x[, 1:3], angle = 40, pch = pch[group], color = colors[group], grid = F, box = F,
              xlab=paste("PC1, ", round(pca.propVar3_c[1], 2), "%"),
              ylab=paste("PC2, ", round(pca.propVar3_c[2], 2), "%"),
              zlab=paste("PC3, ", round(pca.propVar3_c[3], 2), "%"),
              main = "3D-scatterPlot-PCA (lane corrected)")
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca3$x[,1:3], grid = c("xy", "xz", "yz"))
legend("topright", legend=levels(group), pch=pch, col=colors)

scatterplot3d(pca3_c$x[, 1:3], angle = -40, pch = pch[group], color = colors[group], grid = F, box = F,
              xlab=paste("PC1, ", round(pca.propVar3_c[1], 2), "%"),
              ylab=paste("PC2, ", round(pca.propVar3_c[2], 2), "%"),
              zlab=paste("PC3, ", round(pca.propVar3_c[3], 2), "%"),
              main = "3D-scatterPlot-PCA (lane corrected)")
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca3$x[,1:3], grid = c("xy", "xz", "yz"))
legend("topright", legend=levels(group), pch=pch, col=colors)


## Contrast matrix
contr.mat <- makeContrasts(
  ps1isoVSps1 = ps1iso - ps1,
  ps1VSwt = ps1 - wt,
  wtisoVSps1 = wtiso - ps1,
  ps1isoVSwt = ps1iso - wt,
  ps1isoVSwtiso = ps1iso - wtiso,
  wtisoVSwt = wtiso - wt,
  levels = colnames(design))
contr.mat

## Removing heteroscedascity from count data
par(mfrow=c(1,2))
v <- voom(y.norm, design = design, plot = T)
vfit.des <- lmFit(v, design = design)
vfit.contr <- contrasts.fit(vfit.des, contrasts = contr.mat)
efit.contr <- eBayes(vfit.contr, robust = T)
plotSA(efit.contr, main="Final model: Mean−variance trend")
dev.print(device = pdf, "Rplot_meanVarTrend_ps1iso.pdf", height = 6, width = 12)
###
par(mfrow=c(1,2))
v2 <- voom(y.norm2, design = design, plot = T)
vfit2.des <- lmFit(v2, design = design)
vfit2.contr <- contrasts.fit(vfit2.des, contrasts = contr.mat)
efit2.contr <- eBayes(vfit2.contr, robust = T)
plotSA(efit2.contr, main="Final model: Mean−variance trend")
dev.print(device = pdf, "Rplot_meanVarTrend2_ps1iso.pdf", height = 6, width = 12)
###
par(mfrow=c(1,2))
v3 <- voom(y.norm3, design = design, plot = T)
vfit3.des <- lmFit(v3, design = design)
vfit3.contr <- contrasts.fit(vfit3.des, contrasts = contr.mat)
efit3.contr <- eBayes(vfit3.contr, robust = T)
plotSA(efit3.contr, main="Final model: Mean−variance trend")
dev.print(device = pdf, "Rplot_meanVarTrend3_ps1iso.pdf", height = 6, width = 12)


## Using "treat" instead of "eBayes" with lfc=1
tfit3.contr <- treat(vfit3.contr, lfc = 1)
plotSA(tfit3.contr, main="Final model: Mean−variance trend (treat)")
plotSA(treat(vfit3.contr, lfc = 1, robust = T), main="Final model: Mean−variance trend (treat_robust)")
###
t2fit3.contr <- treat(vfit3.contr, lfc = log2(1.3)) # log2(1.3) = 0.3785116
plotSA(t2fit3.contr, main="Final model: Mean−variance trend (treat)")
plotSA(treat(vfit3.contr, lfc = log2(1.3), robust = T), main="Final model: Mean−variance trend (treat_robust)")
###
t3fit3.contr <- treat(vfit3.contr, lfc = log2(1.5)) # log2(1.5) = 0.5849625. i.e. 50% change
plotSA(t3fit3.contr, main="Final model: Mean−variance trend (treat)")
plotSA(treat(vfit3.contr, lfc = log2(1.5), robust = T), main="Final model: Mean−variance trend (treat_robust)")
###
t4fit3.contr <- treat(vfit3.contr, lfc = log2(1.1))



dt <- decideTests(efit.contr)
summary(dt)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down          1080    5730       5881       5737          5669       306
# NotSig       13378    3879       3916       3837          3948     14958
# Up            1049    5898       5710       5933          5890       243
###
dt2 <- decideTests(efit2.contr)
summary(dt2)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down          1029    5581       5829       5553          5492       279
# NotSig       12961    3639       3675       3603          3706     14556
# Up            1090    5860       5576       5924          5882       245
###
dt3 <- decideTests(efit3.contr)
summary(dt3)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down          1078    5758       5924       5766          5687       319
# NotSig       13533    3961       4002       3919          4036     15097
# Up            1047    5939       5732       5973          5935       242
###
dt3_t <- decideTests(tfit3.contr)
summary(dt3_t)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down             0    1235       1703       1198          1169         0
# NotSig       15658   12705      12744      12688         12763     15658
# Up               0    1718       1211       1772          1726         0
###
dt3_t2 <- decideTests(t2fit3.contr)
summary(dt3_t2)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down             3    2940       3082       2873          2868         0
# NotSig       15493    9476       9473       9447          9552     15504
# Up              11    3091       2952       3187          3087         3
###
dt3_t3 <- decideTests(t3fit3.contr)
summary(dt3_t3)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down             0    2112       2341       2029          2000         0
# NotSig       15654   11172      11209      11204         11312     15658
# Up               4    2374       2108       2425          2346         0
###
dt3_t4 <- decideTests(t4fit3.contr)
summary(dt3_t4)
###
dt3_p0.01 <- decideTests(efit3.contr, p.value = 0.01)
summary(dt3_p0.01)
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down           505    5174       5333       5168          5105        60
# NotSig       14672    5162       5185       5132          5208     15517
# Up             481    5322       5140       5358          5345        81


###
summary(decideTests(eBayes(vfit3.contr)))   ### WITHOUT ROBUST = TRUE
#        ps1isoVSps1 ps1VSwt wtisoVSps1 ps1isoVSwt ps1isoVSwtiso wtisoVSwt
# Down          1030    5760       5926       5775          5692       267
# NotSig       13623    3961       4004       3896          4022     15171
# Up            1005    5937       5728       5987          5944       220

###
summary(decideTests(treat(vfit3.contr, lfc = log2(1.5))))
summary(decideTests(treat(vfit3.contr, lfc = log2(1.5)), p.value = 0.01))

####### Robust estimation does not increase the false discovery rate when no outliers are present ###### (Phipson & Smyth et al., 2013 and 2016)

write.fit(efit3.contr, dt3, file = "efit3.contr.txt")

sigUP1.2 <- which(dt3[,4] > 0 & dt3[,2] <= 0)
length(sigUP1.2) # 548

sigDOWN1.2 <- which(dt3[,4] < 0 & dt3[,2] >= 0)
length(sigDOWN1.2) # 539

548 + 539 # 1087

sigUP1.2genes <- efit3.contr$genes$genes[sigUP1.2]
write.table(sigUP1.2genes, file = "GO/UniqueSigUPGeneListWithIsoTxt.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

sigDOWN1.2genes <- efit3.contr$genes$genes[sigDOWN1.2]
write.table(sigDOWN1.2genes, file = "GO/UniqueSigDOWNGeneListWithIsoTxt.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

##################################### Four kind of comparisons which make sense (06-04-19) ########################################
## Comparison # 1 (Correction of significantly downregulated genes in PS1KO vs WT to WT level with ISO txt)
# ps1VSwt: LFC < 0 & FDR < 0.05 compared with ps1isoVSps1: LFC > 0 & FDR < 0.05 (At the same time: ps1isoVSwt: FDR >= 0.05)
summary(dt3)

merged.res3 <- merge(res3.ps1VSwt[, c(1,2,6)], res3.ps1isoVSps1[, c(1,2,6)], by="genes")
head(merged.res3)
dim(merged.res3) # 15658     5
merged.res3 <- merge(merged.res3, res3.ps1isoVSwt[, c(1,2,6)], by="genes")
head(merged.res3)
dim(merged.res3) # 15658     7

corr1 <- merged.res3$genes[merged.res3$logFC.x < 0 & merged.res3$adj.P.Val.x < 0.05 & 
                             merged.res3$logFC.y > 0 & merged.res3$adj.P.Val.y < 0.05 & 
                             merged.res3$adj.P.Val >= 0.05]
length(corr1) # 145
corr1
write.table(corr1, file = "correction_of_SigDn_to_WTWithIso_corr1.txt", sep = "\t", 
            row.names = F, col.names = F, quote = F)

## Comparison # 2 (Correction of significantly upregulated genes in PS1KO vs WT to WT level with ISO txt)
# ps1VSwt: LFC > 0 & FDR < 0.05 compared with ps1isoVSps1: LFC < 0 & FDR < 0.05 (At the same time ps1isoVSwt: FDR >= 0.05)

corr2 <- merged.res3$genes[merged.res3$logFC.x > 0 & merged.res3$adj.P.Val.x < 0.05 & 
                             merged.res3$logFC.y < 0 & merged.res3$adj.P.Val.y < 0.05 & 
                             merged.res3$adj.P.Val >= 0.05]
length(corr2) # 159
corr2
write.table(corr2, file = "correction_of_SigUp_to_WTWithIso_corr2.txt", sep = "\t", 
            row.names = F, col.names = F, quote = F)

## Comparison # 3 (No sig downregulation in PS1 vs WT but significantly upregulated with ISO txt but NOT BEYOND WT level.)
# ps1VSwt: LFC < 0 & FDR >= 0.05 compared with ps1isoVSps1: LFC > 0 & FDR < 0.05 (at the same time ps1isoVSwt: FDR >= 0.05)

corr3 <- merged.res3$genes[merged.res3$logFC.x < 0 & merged.res3$adj.P.Val.x >= 0.05 & 
                             merged.res3$logFC.y > 0 & merged.res3$adj.P.Val.y < 0.05 & 
                             merged.res3$adj.P.Val >= 0.05]
length(corr3) # 15
corr3

## Comparison # 4 (No sig upregulation in PS1 vs WT but significantly downregulated with ISO txt but NOT BEYOND WT level.)
# ps1VSwt: LFC > 0 & FDR >= 0.05 compared with ps1isoVSps1: LFC < 0 & FDR < 0.05 (at the same time ps1isoVSwt: FDR >= 0.05)

corr4 <- merged.res3$genes[merged.res3$logFC.x > 0 & merged.res3$adj.P.Val.x >= 0.05 & 
                             merged.res3$logFC.y < 0 & merged.res3$adj.P.Val.y < 0.05 & 
                             merged.res3$adj.P.Val >= 0.05]
length(corr4) # 22
corr4

## Comparison # 5 (1.1) (Significant Correction of significantly downregulated genes in PS1KO vs WT but NOT neccessarily to WT level with ISO txt. Can be ABOVE WT LEVEL)
# ps1VSwt: LFC < 0 & FDR < 0.05 compared with ps1isoVSps1: LFC > 0 & FDR < 0.05 

corr5 <- merged.res3$genes[merged.res3$logFC.x < 0 & merged.res3$adj.P.Val.x < 0.05 & 
                             merged.res3$logFC.y > 0 & merged.res3$adj.P.Val.y < 0.05]
length(corr5) # 471
write.table(corr5, file = "sigUpWithIso_corr5.txt", sep = "\t", 
            row.names = F, col.names = F, quote = F)

## Comparison # 6 (2.1) (Significant Correction of significantly upregulated genes in PS1KO vs WT but NOT neccessarily to WT level with ISO txt. Can be BELOW WT LEVEL)
# ps1VSwt: LFC > 0 & FDR < 0.05 compared with ps1isoVSps1: LFC < 0 & FDR < 0.05 

corr6 <- merged.res3$genes[merged.res3$logFC.x > 0 & merged.res3$adj.P.Val.x < 0.05 & 
                             merged.res3$logFC.y < 0 & merged.res3$adj.P.Val.y < 0.05]
length(corr6) # 586
write.table(corr6, file = "sigDnWithIso_corr5.txt", sep = "\t", 
            row.names = F, col.names = F, quote = F)




###################################################################################################################################

## Venn-diagram
par(mfrow=c(1,1))
vennDiagram(dt3[,c(2,4)], circle.col = c("turquoise", "salmon"))
dev.print(pdf, "Rplot_Venn3_ps1VSwt-ps1isoVSwt.pdf", height=6, width=8)
###
vennDiagram(dt3_t[,c(2,4)], circle.col = c("turquoise", "salmon"))
dev.print(pdf, "Rplot_Venn3_t_ps1VSwt-ps1isoVSwt.pdf", height=6, width=8)
###
vennDiagram(dt3_t2[,c(2,4)], circle.col = c("turquoise", "salmon"))
dev.print(pdf, "Rplot_Venn3_t2_ps1VSwt-ps1isoVSwt.pdf", height=6, width=8)

vennDiagram(dt3[,c(2,1)], circle.col = c("turquoise", "salmon"), include = c("up", "down"))
dev.print(pdf, "Rplot_Venn3_ps1VSwt-ps1isoVSps1.pdf", height=6, width=8)



## Examining individual DE genes from top to bottom
res.ps1VSwt <- topTable(efit.contr, coef = 2, n = Inf)
head(res.ps1VSwt)
write.table(res.ps1VSwt, file = "res.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F) # To not show rows which is a duplicate of "genes" column
#write.table(res.ps1VSwt, file = "res.ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F) # To save df as it is.

res.ps1isoVSwt <- topTable(efit.contr, coef = 4, n = Inf)
head(res.ps1isoVSwt)
write.table(res.ps1isoVSwt, file = "res.ps1isoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res.ps1isoVSps1 <- topTable(efit.contr, coef = 1, n = Inf)
head(res.ps1isoVSps1)
write.table(res.ps1isoVSps1, file = "res.ps1isoVSps1.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res.ps1isoVSwtiso <- topTable(efit.contr, coef = 5, n = Inf)
head(res.ps1isoVSwtiso)
write.table(res.ps1isoVSwtiso, file = "res.ps1isoVSwtiso.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res.wtisoVSwt <- topTable(efit.contr, coef = 6, n = Inf)
head(res.wtisoVSwt)
write.table(res.wtisoVSwt, file = "res.wtisoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
res2.ps1VSwt <- topTable(efit2.contr, coef = 2, n = Inf)
head(res2.ps1VSwt)
write.table(res2.ps1VSwt, file = "res2.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F) # To not show rows which is a duplicate of "genes" column

#write.table(res.ps1VSwt, file = "res.ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F) # To save df as it is.

res2.ps1isoVSwt <- topTable(efit2.contr, coef = 4, n = Inf)
head(res2.ps1isoVSwt)
write.table(res2.ps1isoVSwt, file = "res2.ps1isoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res2.ps1isoVSps1 <- topTable(efit2.contr, coef = 1, n = Inf)
head(res2.ps1isoVSps1)
write.table(res2.ps1isoVSps1, file = "res2.ps1isoVSps1.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res2.ps1isoVSwtiso <- topTable(efit2.contr, coef = 5, n = Inf)
head(res2.ps1isoVSwtiso)
write.table(res2.ps1isoVSwtiso, file = "res2.ps1isoVSwtiso.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res2.wtisoVSwt <- topTable(efit2.contr, coef = 6, n = Inf)
head(res2.wtisoVSwt)
write.table(res2.wtisoVSwt, file = "res2.wtisoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
res3.ps1VSwt <- topTable(efit3.contr, coef = 2, n = Inf)
head(res3.ps1VSwt)
write.table(res3.ps1VSwt, file = "res3.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F) # To not show rows which is a duplicate of "genes" column

#write.table(res.ps1VSwt, file = "res.ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F) # To save df as it is.

res3.ps1isoVSwt <- topTable(efit3.contr, coef = 4, n = Inf)
head(res3.ps1isoVSwt)
write.table(res3.ps1isoVSwt, file = "res3.ps1isoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSps1 <- topTable(efit3.contr, coef = 1, n = Inf)
head(res3.ps1isoVSps1)
write.table(res3.ps1isoVSps1, file = "res3.ps1isoVSps1.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSwtiso <- topTable(efit3.contr, coef = 5, n = Inf)
head(res3.ps1isoVSwtiso)
write.table(res3.ps1isoVSwtiso, file = "res3.ps1isoVSwtiso.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.wtisoVSwt <- topTable(efit3.contr, coef = 6, n = Inf)
head(res3.wtisoVSwt)
write.table(res3.wtisoVSwt, file = "res3.wtisoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
res3_t.ps1VSwt <- topTreat(tfit3.contr, coef = 2, n = Inf)
head(res3_t.ps1VSwt)
write.table(res3_t.ps1VSwt, file = "res3_t.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3_t.ps1isoVSwt <- topTreat(tfit3.contr, coef = 4, n = Inf)
head(res3_t.ps1isoVSwt)
write.table(res3_t.ps1isoVSwt, file = "res3_t.ps1isoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
res3_t2.ps1VSwt <- topTreat(t2fit3.contr, coef = 2, n = Inf)
head(res3_t2.ps1VSwt)
write.table(res3_t2.ps1VSwt, file = "res3_t2.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3_t2.ps1isoVSwt <- topTreat(t2fit3.contr, coef = 4, n = Inf)
head(res3_t2.ps1isoVSwt)
write.table(res3_t2.ps1isoVSwt, file = "res3_t2.ps1isoVSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.wtisoVSps1 <- topTable(efit3.contr, coef = 3, n = Inf)
head(res3.wtisoVSps1)
###
res3_t3.ps1VSwt <- topTreat(t3fit3.contr, coef = 2, n = Inf)
head(res3_t3.ps1VSwt)
write.table(res3_t3.ps1VSwt, file = "res3_t3.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)


## MD plots
par(mfrow=c(1,1))
plotMD(efit.contr, coef = 1, status = dt[,1], main = colnames(efit.contr)[1], xlim = c(-8,13))
# or
plotMD(efit.contr, coef = 1, status = dt[,1], main = "ps1.iso vs ps1", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot_ps1isoVSps1.pdf", height = 6, width = 8)
###
par(mfrow=c(1,1))
plotMD(efit2.contr, coef = 1, status = dt2[,1], main = "ps1.iso vs ps1", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot2_ps1isoVSps1.pdf", height = 6, width = 8)
###
par(mfrow=c(1,1))
plotMD(efit3.contr, coef = 1, status = dt3[,1], main = "ps1.iso vs ps1", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot3_ps1isoVSps1.pdf", height = 6, width = 8)


# Interactive MD Plot
glMDPlot(efit3.contr, coef = 1, status = dt3[,1], main = "ps1.iso vs ps1", counts = lcpm.norm3, groups = group, side.main = "symbol",
         folder = "Glimma-MDPlot-ps1isoVSps1", html = "glMDPlot_ps1isoVSps1", launch = T)

glMDPlot(efit3.contr, coef = 2, status = dt3[,2], main = "ps1 vs wt", counts = lcpm.norm3, groups = group, side.main = "symbol",
         folder = "Glimma-MDPlot-ps1VSwt", html = "glMDPlot_ps1VSwt", launch = T)

glMDPlot(efit3.contr, coef = 4, status = dt3[,4], main = "ps1.iso vs wt", counts = lcpm.norm3, groups = group, side.main = "symbol",
         folder = "Glimma-MDPlot-ps1isoVSwt", html = "glMDPlot_ps1isoVwt", launch = T)


### Other MD plots for various comparisons
plotMD(efit.contr, coef = 2, status = dt[,2], main = "ps1 vs wt", xlim = c(-8,13), values = c(1,-1))
#dev.print(pdf, "Rplot_MDPlot_ps1VSwt.pdf", height = 6, width = 8)
###
plotMD(efit3.contr, coef = 2, status = dt3[,2], main = "ps1 vs wt", xlim = c(-8,13), values = c(1,-1))
dev.print(pdf, "Rplot_MDPlot3_ps1VSwt.pdf", height = 6, width = 8)
###
plotMD(tfit3.contr, coef = 2, status = dt3_t[,2], main = "ps1 vs wt", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot3_t_ps1isoVSps1.pdf", height = 6, width = 8)
###
plotMD(t2fit3.contr, coef = 2, status = dt3_t2[,2], main = "ps1 vs wt", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot3_t2_ps1isoVSps1.pdf", height = 6, width = 8)


plotMD(efit.contr, coef = 4, status = dt[,4], main = "ps1.iso vs wt", xlim = c(-8,13), values = c(1,-1))
#dev.print(pdf, "Rplot_MDPlot_ps1isoVSwt.pdf", height = 6, width = 8)
###
plotMD(efit3.contr, coef = 4, status = dt3[,4], main = "ps1.iso vs wt", xlim = c(-8,13), values = c(1,-1))
dev.print(pdf, "Rplot_MDPlot3_ps1isoVSwt.pdf", height = 6, width = 8)
###
plotMD(tfit3.contr, coef = 4, status = dt3_t[,4], main = "ps1.iso vs wt", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot3_t_ps1isoVSwt.pdf", height = 6, width = 8)
###
plotMD(t2fit3.contr, coef = 4, status = dt3_t2[,4], main = "ps1.iso vs wt", xlim = c(-8,13))
dev.print(pdf, "Rplot_MDPlot3_t2_ps1isoVSwt.pdf", height = 6, width = 8)


plotMD(efit.contr, coef = 5, status = dt[,5], main = "ps1.iso vs wt.iso", xlim = c(-8,13))
#dev.print(pdf, "Rplot_MDPlot_ps1isoVSwtiso.pdf", height = 6, width = 8)


plotMD(efit.contr, coef = 6, status = dt[,6], main = "wt.iso vs wt", xlim = c(-8,13))
#dev.print(pdf, "Rplot_MDPlot_wtisoVSwt.pdf", height = 6, width = 8)


## Heatmap
library(gplots)
par(mfrow=c(1,1))
ps1isoVSps1.topgenes <- res.ps1isoVSps1$genes[1:100]
i <- which(v$genes$genes %in% ps1isoVSps1.topgenes)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap_ps1isoVSps1.pdf", height=10, width=8)
###
library(gplots)
par(mfrow=c(1,1))
ps1isoVSps1.topgenes2 <- res2.ps1isoVSps1$genes[1:100]
i <- which(v2$genes$genes %in% ps1isoVSps1.topgenes2)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v2$E[i,], scale = "row",
          labRow = v2$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap2_ps1isoVSps1.pdf", height=10, width=8)
###
library(gplots)
par(mfrow=c(1,1))
ps1isoVSps1.topgenes3 <- res3.ps1isoVSps1$genes[1:100]
i <- which(v3$genes$genes %in% ps1isoVSps1.topgenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = v3$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap3_ps1isoVSps1.pdf", height=10, width=8)


ps1VSwt.topgenes <- res.ps1VSwt$genes[1:100]
i <- which(v$genes$genes %in% ps1VSwt.topgenes)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v$E[i,], scale = "row",
          labRow = v$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap_ps1VSwt.pdf", height=10, width=8)
###
ps1VSwt.topgenes3 <- res3.ps1VSwt$genes[1:100]
i <- which(v3$genes$genes %in% ps1VSwt.topgenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = v3$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap3_ps1VSwt.pdf", height=10, width=8)


ps1isoVSwt.topgenes3 <- res3.ps1isoVSwt$genes[1:100]
i <- which(v3$genes$genes %in% ps1isoVSwt.topgenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = v3$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap3_ps1isoVSwt.pdf", height=10, width=8)


ps1isoVSwtiso.topgenes3 <- res3.ps1isoVSwtiso$genes[1:100]
i <- which(v3$genes$genes %in% ps1isoVSwtiso.topgenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = v3$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap3_ps1isoVSwtiso.pdf", height=10, width=8)

## Probably the heatmaps for unnecessary comparisons
wtisoVSps1.topgenes3 <- res3.wtisoVSps1$genes[1:100]
i <- which(v3$genes$genes %in% wtisoVSps1.topgenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = v3$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")


wtisoVSwt.topgenes3 <- res3.wtisoVSwt$genes[1:100]
i <- which(v3$genes$genes %in% wtisoVSwt.topgenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = v3$genes$genes[i], labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")


## Heatmap for overall expression profile (all genes rather than top 100)
par(mfrow=c(1,1))
heatmap.2(v3$E, scale = "row",
          labRow = NULL, labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap_ps1VSwt.pdf", height=10, width=8)

### Gene Ontology
## Detected genes in the experiment
genes.detected <- y.dge$counts
head(genes.detected)
dim(genes.detected) # 52459    12

genes.detected <- genes.detected[!(rowSums(genes.detected==0)==12), ]
head(genes.detected)
dim(genes.detected) # 30588    12
write.table(rownames(genes.detected), file = "GO/allGenesDetectedinExpt_NonZeroReadCts.txt", sep = "\t",
            row.names = F, col.names = "symbols", quote = F)

## Highly enriched genes across groups
head(cpm.norm3)
dim(cpm.norm3) # 15658    12
write.table(rownames(cpm.norm3), file = "GO/enrichedGenes_allGroups_15658.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

cpm.norm3_wt <- cpm.norm3[, c(7,8,9)]
table(rowSums(cpm.norm3_wt) > 100)
# FALSE  TRUE
# 9489  6169

cpm.norm3_wt100 <- cpm.norm3_wt[(rowSums(cpm.norm3_wt) > 100), ]
dim(cpm.norm3_wt100) # 6169    3
write.table(rownames(cpm.norm3_wt100), file = "GO/cpm100_in_wt_6169.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)

## All genes regardless of detected or not in illumina chip
all.illum.genes <- y.dge$counts
head(all.illum.genes)
dim(all.illum.genes) # 52459    12
write.table(rownames(all.illum.genes), file = "GO/allIlluminaGenes.txt", sep = "\t",
            row.names = F, col.names = F, quote = F)


head(res3.ps1VSwt)
res3.ps1VSwt_sigUP0.05 <- res3.ps1VSwt[(res3.ps1VSwt$logFC > 0 & res3.ps1VSwt$adj.P.Val <= 0.05), ]
res3.ps1VSwt_sigUP0.05 <- res3.ps1VSwt_sigUP0.05[order(res3.ps1VSwt_sigUP0.05$adj.P.Val), ]
head(res3.ps1VSwt_sigUP0.05)
tail(res3.ps1VSwt_sigUP0.05)
dim(res3.ps1VSwt_sigUP0.05) # 5939    7
write.table(res3.ps1VSwt_sigUP0.05, file = "GO/res3.ps1VSwt_sigUP_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.ps1VSwt_sigDOWN0.05 <- res3.ps1VSwt[(res3.ps1VSwt$logFC < 0 & res3.ps1VSwt$adj.P.Val <= 0.05), ]
res3.ps1VSwt_sigDOWN0.05 <- res3.ps1VSwt_sigDOWN0.05[order(res3.ps1VSwt_sigDOWN0.05$adj.P.Val), ]
head(res3.ps1VSwt_sigDOWN0.05)
tail(res3.ps1VSwt_sigDOWN0.05)
dim(res3.ps1VSwt_sigDOWN0.05) # 5758    7
write.table(res3.ps1VSwt_sigDOWN0.05, file = "GO/res3.ps1VSwt_sigDOWN_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
head(res3.ps1isoVSwt)
res3.ps1isoVSwt_sigUP0.05 <- res3.ps1isoVSwt[(res3.ps1isoVSwt$logFC > 0 & res3.ps1isoVSwt$adj.P.Val <= 0.05), ]
res3.ps1isoVSwt_sigUP0.05 <- res3.ps1isoVSwt_sigUP0.05[order(res3.ps1isoVSwt_sigUP0.05$adj.P.Val), ]
head(res3.ps1isoVSwt_sigUP0.05)
tail(res3.ps1isoVSwt_sigUP0.05)
write.table(res3.ps1isoVSwt_sigUP0.05, file = "GO/res3.ps1isoVSwt_sigUP_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSwt_sigDOWN0.05 <- res3.ps1isoVSwt[(res3.ps1isoVSwt$logFC < 0 & res3.ps1isoVSwt$adj.P.Val <= 0.05), ]
res3.ps1isoVSwt_sigDOWN0.05 <- res3.ps1isoVSwt_sigDOWN0.05[order(res3.ps1isoVSwt_sigDOWN0.05$adj.P.Val), ]
head(res3.ps1isoVSwt_sigDOWN0.05)
tail(res3.ps1isoVSwt_sigDOWN0.05)
dim(res3.ps1isoVSwt_sigDOWN0.05) # 5766    7
write.table(res3.ps1isoVSwt_sigDOWN0.05, file = "GO/res3.ps1isoVSwt_sigDOWN_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
head(res3.ps1isoVSps1)
res3.ps1isoVSps1_sigUP0.05 <- res3.ps1isoVSps1[(res3.ps1isoVSps1$logFC > 0 & res3.ps1isoVSps1$adj.P.Val < 0.05), ]
res3.ps1isoVSps1_sigUP0.05 <- res3.ps1isoVSps1_sigUP0.05[order(res3.ps1isoVSps1_sigUP0.05$adj.P.Val), ]
head(res3.ps1isoVSps1_sigUP0.05)
tail(res3.ps1isoVSps1_sigUP0.05)
dim(res3.ps1isoVSps1_sigUP0.05) # 1047    7
write.table(res3.ps1isoVSps1_sigUP0.05, file = "GO/res3.ps1isoVSps1_sigUP_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSps1_sigDOWN0.05 <- res3.ps1isoVSps1[(res3.ps1isoVSps1$logFC < 0 & res3.ps1isoVSps1$adj.P.Val <= 0.05), ]
res3.ps1isoVSps1_sigDOWN0.05 <- res3.ps1isoVSps1_sigDOWN0.05[order(res3.ps1isoVSps1_sigDOWN0.05$adj.P.Val), ]
head(res3.ps1isoVSps1_sigDOWN0.05)
tail(res3.ps1isoVSps1_sigDOWN0.05)
dim(res3.ps1isoVSps1_sigDOWN0.05) # 1078    7
write.table(res3.ps1isoVSps1_sigDOWN0.05, file = "GO/res3.ps1isoVSps1_sigDOWN_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
###
head(res3_t2.ps1VSwt)
res3_t2.ps1VSwt_sigUP0.05 <- res3_t2.ps1VSwt[(res3_t2.ps1VSwt$logFC > 0 & res3_t2.ps1VSwt$adj.P.Val < 0.05), ]
#res3_t2.ps1VSwt_sigUP0.05 <- res3_t2.ps1VSwt_sigUP0.05[order(res3_t2.ps1VSwt_sigUP0.05$adj.P.Val), ]
head(res3_t2.ps1VSwt_sigUP0.05)
tail(res3_t2.ps1VSwt_sigUP0.05)
dim(res3_t2.ps1VSwt_sigUP0.05) # 3091    6
write.table(res3_t2.ps1VSwt_sigUP0.05, file = "GO/res3_t2.ps1VSwt_sigUP_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)


###
head(res3_t3.ps1VSwt)
res3_t3.ps1VSwt_sigUP0.05 <- res3_t3.ps1VSwt[(res3_t3.ps1VSwt$logFC > 0 & res3_t3.ps1VSwt$adj.P.Val < 0.05), ]
#res3.ps1isoVSwt_sigUP0.05 <- res3.ps1isoVSwt_sigUP0.05[order(res3.ps1isoVSwt_sigUP0.05$adj.P.Val), ]
head(res3_t3.ps1VSwt_sigUP0.05)
tail(res3_t3.ps1VSwt_sigUP0.05)
write.table(res3_t3.ps1VSwt_sigUP0.05, file = "GO/res3_t3.ps1VSwt_sigUP_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

res3_t3.ps1VSwt_sigDN0.05 <- res3_t3.ps1VSwt[(res3_t3.ps1VSwt$logFC < 0 & res3_t3.ps1VSwt$adj.P.Val < 0.05), ]
res3_t3.ps1isoVSwt_sigDN0.05 <- res3_t3.ps1VSwt_sigDN0.05[order(res3_t3.ps1VSwt_sigDN0.05$adj.P.Val), ]
head(res3_t3.ps1VSwt_sigDN0.05)
tail(res3_t3.ps1VSwt_sigDN0.05)
write.table(res3_t3.ps1VSwt_sigDN0.05, file = "GO/res3_t3.ps1VSwt_sigDN_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

## Let's take top ~1100 for both up and down categories from ps1VSwt to feed into enrichr so fairly compare with ps1isoVSps1 since latter comparison has ~1100 in both
## up and down categories.

## Determining top genes based on both logFC and FDR




###########################   For GO analysis (goana and Kegga), we will name y.norm3.1 and so on     ################################################
y.dge
y.filt3
dim(y.filt3) # 15658    12
library(AnnotationDbi)
entrez <- select(org.Mm.eg.db, keys = y.filt3$genes$genes, keytype = "SYMBOL", columns = "ENTREZID")
    # 'select()' returned 1:many mapping between keys and columns
head(entrez)
length(unique(entrez$ENTREZID)) # 14608
length(unique(y.filt$genes$genes)) # 15507

entrez <- entrez[!duplicated(entrez$ENTREZID), ]
dim(entrez) # 14608     2
entrez <- entrez[!is.na(entrez$ENTREZID), ]
dim(entrez) # 14607     2

sum(duplicated(entrez$SYMBOL)) # 1
sum(is.na(entrez$SYMBOL)) # 0

entrez <- entrez[!duplicated(entrez$SYMBOL), ]
dim(entrez) # 14606     2

sum(duplicated(y.filt3.1$genes$ENTREZID))

y.filt3.1 <- y.filt3
y.filt3.1 <- y.filt3.1[rownames(y.filt3.1) %in% entrez$SYMBOL, ]
y.filt3.1
dim(y.filt3.1) # 14606    12

entrez <- entrez[entrez$SYMBOL %in% rownames(y.filt3.1), ]
dim(entrez)

y.filt3.1$genes$entrezID <- entrez$ENTREZID[match(y.filt3.1$genes$genes, entrez$SYMBOL)]
y.filt3.1

#y.filt3.1 <- y.filt3.1[rownames(y.filt3.1) %in% y.filt3.1$genes$SYMBOL, ]
#dim(y.filt3.1) # 14606    12
#y.filt3.1

isTRUE(all.equal(rownames(y.filt3.1), y.filt3.1$genes$genes)) # TRUE

rownames(y.filt3.1) <- y.filt3.1$genes$entrezID
y.filt3.1

lcpm.filt3.1 <- cpm(y.filt3.1, log=TRUE)
plot(density(lcpm.filt3.1[,1]), col=col[1], lwd=2, ylim=c(0,1), las=2,
     main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm3.cutoff, lty=3)
for (i in 2:nsamples3){
  den <- density(lcpm.filt3.1[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n", ncol = 2)

## Normalising gene expression distributions
y.norm3.1 <- calcNormFactors(y.filt3.1, method = "TMM")
y.norm3.1$samples

## visualization of normalization effect
lcpm.norm3.1 <- cpm(y.norm3.1, log=T)
boxplot(lcpm.norm3.1, las=2, col=col, main="B. Example: Normalised data", ylab="Log-cpm")
abline(h=median(lcpm.norm3.1), col="red")

## Unsupervised clustering of samples
par(mfrow=c(1,1))
library(RColorBrewer)
colors <- brewer.pal(4, "Dark2")
plotMDS(lcpm.norm3.1, col=colors[group], main="PCA plot (Unsupervised Clustering)")

## Removing heteroscedascity from count data
par(mfrow=c(1,2))
v3.1 <- voom(y.norm3.1, design = design, plot = T)
vfit3.1.des <- lmFit(v3.1, design = design)
vfit3.1.contr <- contrasts.fit(vfit3.1.des, contrasts = contr.mat)
efit3.1.contr <- eBayes(vfit3.1.contr, robust = T)
plotSA(efit3.1.contr, main="Final model: Mean−variance trend")

## GO
go3.1_ps1VSwt <- goana(efit3.1.contr[,2], species = "Mm")
topGO(go3.1_ps1VSwt, n=30, ontology = "BP")
# try
topGO(goana(efit3.1.contr[,2], species = "Mm", trend = T), n=30, ontology = "BP", sort = "down")
topGO(goana(efit3.1.contr[,2], species = "Mm", trend = T), n=30, ontology = "BP", sort = "up")

go3.1_ps1isoVSwt <- goana(efit3.1.contr[,4], species = "Mm")
topGO(go3.1_ps1isoVSwt, n=75, ontology = "BP")
#topGO(go3.1_ps1isoVSwt, ontology = "BP") [31:50]

topGO(goana(efit3.1.contr[,1], species = "Mm", trend = T), n=40, ontology = "BP")

keg3.1_ps1VSwt <- kegga(efit3.1.contr[,2], species = "Mm")
topKEGG(keg3.1_ps1VSwt, n=40, truncate.path = 60)
topKEGG( kegga(efit3.1.contr[,2], species = "Mm", trend = T), n=40, truncate.path = 60)

keg3.1_ps1isoVSwt <- kegga(efit3.1.contr[,4], species = "Mm")
topKEGG(keg3.1_ps1isoVSwt, n=40, truncate.path = 60)

topKEGG( kegga(efit3.1.contr[,1], species = "Mm", trend = T), n=40, truncate.path = 60)

res3.1.ps1isoVSps1 <- topTable(efit3.1.contr, coef = 1, n=Inf)
head(res3.1.ps1isoVSps1)

res3.1.ps1VSwt <- topTable(efit3.1.contr, coef = 2, n=Inf)
head(res3.1.ps1VSwt)
write.table(res3.1.ps1VSwt, file = "res3.1.ps1VSwt.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

###########################   goana and kegga finished  ###########################
############################### FOR CLUE Connectopedia ###############################
head(res3.1.ps1VSwt)
res3.1.ps1VSwt_sigUP0.05 <- res3.1.ps1VSwt[(res3.1.ps1VSwt$logFC > 0 & res3.1.ps1VSwt$adj.P.Val < 0.05), ]
res3.1.ps1VSwt_sigUP0.05 <- res3.1.ps1VSwt_sigUP0.05[order(res3.1.ps1VSwt_sigUP0.05$adj.P.Val), ]
head(res3.1.ps1VSwt_sigUP0.05)
tail(res3.1.ps1VSwt_sigUP0.05)
dim(res3.1.ps1VSwt_sigUP0.05) # 5729    8
write.table(res3.1.ps1VSwt_sigUP0.05, file = "CLUE/res3.1.ps1VSwt_sigUP_padjCutOff0.05.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)



########################### FOR CLUE Connectopedia completed ##########################
#### Gene set testing
## Gene Set Testing done on 02-26-2019 on geneSets110 (gs100:1-29-18 + gsCLEAR: 12-18-18 (n=10)) 
geneSets110_Hs <- loadWorkbook("~/NKI/Gene Sets/R dir_GeneSets/geneSets110_Hs withNum_121818.xlsx")
geneSets110_Hs <- readWorksheet(geneSets110_Hs, sheet = 1, header = T)
head(geneSets110_Hs)
tail(geneSets110_Hs)
str(geneSets110_Hs)

# Need to convert into mouse symbols (prop cases)
geneSets110_Mm <- data.frame(lapply(geneSets110_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v,1,1)), tolower(substr(v,2,20))), v)
}), stringsAsFactors = F)

class(geneSets110_Mm) # "data.frame"
class(geneSets110_Mm$NfKb.pathway.23.) # "character"
write.table(geneSets110_Mm, file = "geneSets110_Mm withNum_121818.txt", sep = "\t",
            row.names = F, col.names = T, quote = F, na = "")
write.table(geneSets110_Mm, file = "~/NKI/Gene Sets/R dir_GeneSets/geneSets110_Mm withNum_121818.txt", sep = "\t",
            row.names = F, col.names = T, quote = F, na = "")


## convert into ids
idx.all.gs110 <- ids2indices(geneSets110_Mm, id = rownames(v))
str(idx.all.gs110)
###
idx2.all.gs110 <- ids2indices(geneSets110_Mm, id = rownames(v2))
str(idx2.all.gs110)
###
idx3.all.gs110 <- ids2indices(geneSets110_Mm, id = rownames(v3))
str(idx3.all.gs110)


## fry with gs110 (12-18-2018) on 02-26-2019
fry.gs100_ps1isoVSps1 <- fry(v, idx.all.gs110, design, contrast = contr.mat[,1])
fry.gs100_ps1isoVSps1
###
fry2.gs100_ps1isoVSps1 <- fry(v2, idx2.all.gs110, design, contrast = contr.mat[,1])
fry2.gs100_ps1isoVSps1
###
fry3.gs100_ps1isoVSps1 <- fry(v3, idx3.all.gs110, design, contrast = contr.mat[,1])
fry3.gs100_ps1isoVSps1

fry3.gs100_ps1VSwt <- fry(v3, idx3.all.gs110, design, contrast = contr.mat[,2])
fry3.gs100_ps1VSwt

## mroast with gs110 (12-18-2018) on 02-26-2019
mroast.gs100_ps1VSwt <- mroast(v, idx.all.gs110, design, contrast = contr.mat[,2], nrot = 100000)
mroast.gs100_ps1VSwt
write.table(mroast.gs100_ps1VSwt, file = "GSEA/mroast.gs100_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
###
mroast2.gs100_ps1VSwt <- mroast(v2, idx2.all.gs110, design, contrast = contr.mat[,2], nrot = 100000)
mroast2.gs100_ps1VSwt
write.table(mroast2.gs100_ps1VSwt, file = "GSEA/mroast2.gs100_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
###
mroast3.gs100_ps1VSwt <- mroast(v3, idx3.all.gs110, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.gs100_ps1VSwt
write.table(mroast3.gs100_ps1VSwt, file = "GSEA/mroast3.gs100_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


mroast.gs100_ps1isoVSwtiso <- mroast(v, idx.all.gs110, design, contrast = contr.mat[,5], nrot = 100000)
mroast.gs100_ps1isoVSwtiso
write.table(mroast.gs100_ps1isoVSwtiso, file = "GSEA/mroast.gs100_ps1isoVSwtiso.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
###
mroast3.gs100_ps1isoVSwtiso <- mroast(v3, idx3.all.gs110, design, contrast = contr.mat[,5], nrot = 100000)
mroast3.gs100_ps1isoVSwtiso
write.table(mroast3.gs100_ps1isoVSwtiso, file = "GSEA/mroast3.gs100_ps1isoVSwtiso.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


mroast.gs100_ps1isoVSps1 <- mroast(v, idx.all.gs110, design, contrast = contr.mat[,1], nrot = 100000)
mroast.gs100_ps1isoVSps1
write.table(mroast.gs100_ps1isoVSps1, file = "GSEA/mroast.gs100_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
###
mroast3.gs100_ps1isoVSps1 <- mroast(v3, idx3.all.gs110, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.gs100_ps1isoVSps1
write.table(mroast3.gs100_ps1isoVSps1, file = "GSEA/mroast3.gs100_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


mroast.gs100_ps1isoVSwt <- mroast(v, idx.all.gs110, design, contrast = contr.mat[,4], nrot = 100000)
mroast.gs100_ps1isoVSwt
write.table(mroast.gs100_ps1isoVSwt, file = "GSEA/mroast.gs100_ps1isoVSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
###
mroast2.gs100_ps1isoVSwt <- mroast(v2, idx2.all.gs110, design, contrast = contr.mat[,4], nrot = 100000)
mroast2.gs100_ps1isoVSwt
write.table(mroast2.gs100_ps1isoVSwt, file = "GSEA/mroast2.gs100_ps1isoVSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
###
mroast3.gs100_ps1isoVSwt <- mroast(v3, idx3.all.gs110, design, contrast = contr.mat[,4], nrot = 100000)
mroast3.gs100_ps1isoVSwt
write.table(mroast3.gs100_ps1isoVSwt, file = "GSEA/mroast3.gs100_ps1isoVSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


mroast3.gs100_wtisoVSwt <- mroast(v3, idx3.all.gs110, design, contrast = contr.mat[,6], nrot = 100000)
mroast3.gs100_wtisoVSwt
write.table(mroast3.gs100_wtisoVSwt, file = "GSEA/mroast3.gs100_wtisoVSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


## Barcode plots
par(mfrow=c(1,2))
barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$Curated.Induction.23., main = "ps1 vs wt: Curated Induction (n=23)")
barcodeplot(efit3.contr$t[,4], index = idx3.all.gs110$Curated.Induction.23., main = "ps1.iso vs wt: Curated Induction (n=23)")
dev.print(pdf, "GSEA/Rplot_barcode3_curatedInd_ps1VSwt_ps1isoVSwt.pdf", height=6, width=12)
barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$Curated.Induction.23., main = "ps1.iso vs ps1: Curated Induction (n=23)")
dev.print(pdf, "GSEA/Rplot_barcode3_curatedInd_ps1isoVSps1.pdf", height=6, width=8)
dev.print(pdf, "GSEA/Rplot_barcode3_curatedInd_ps1VSwt_ps1isoVSps1.pdf", height=6, width=12)

barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$New.MT.Assoc.All.minus.motors.53., main = "ps1 vs wt: New.MT.Assoc.All.minus.motors (n=53)")
barcodeplot(efit3.contr$t[,4], index = idx3.all.gs110$New.MT.Assoc.All.minus.motors.53., main = "ps1.iso vs wt: New.MT.Assoc.All.minus.motors.53. (n=53)")
dev.print(pdf, "GSEA/Rplot_barcode3_MTminusMotors_ps1VSwt_ps1isoVSwt.pdf", height=6, width=12)

barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$vesicle.mediated.transport.16., main = "ps1 vs wt: Vesicle mediated transport (n=16)")
barcodeplot(efit3.contr$t[,4], index = idx3.all.gs110$vesicle.mediated.transport.16., main = "ps1.iso vs wt: Vesicle mediated transport (n=16)")
dev.print(pdf, "GSEA/Rplot_barcode3_vesclMedTransport_ps1VSwt_ps1isoVSwt.pdf", height=6, width=12)

barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$All.trasport.genes.232., main = "ps1 vs wt: All Transport Genes (n=232)")
barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$All.trasport.genes.232., main = "ps1.iso vs ps1: All Transport Genes (n=232)")
dev.print(pdf, "GSEA/Rplot_barcode3_AllTransport_ps1VSwt_ps1isoVSps1.pdf", height=6, width=12)

barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$protein.transport.40., main = "ps1 vs wt: Protein Transport (n=40)")
barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$protein.transport.40., main = "ps1.iso vs ps1: Protein Transport (n=40)")
dev.print(pdf, "GSEA/Rplot_barcode3_ProtTransport_ps1VSwt_ps1isoVSps1.pdf", height=6, width=12)

barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$Upregulated.upon.TFEB.OE.in.Hela.291., main = "ps1 vs wt: Upregulated.upon.TFEB.OE.in.Hela (n=291)")
barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$Upregulated.upon.TFEB.OE.in.Hela.291., main = "ps1.iso vs ps1: Upregulated.upon.TFEB.OE.in.Hela (n=291)")
dev.print(pdf, "GSEA/Rplot_barcode3_UpwithTFEB.OE_ps1VSwt_ps1isoVSps1.pdf", height=6, width=12)

barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$endosome.to.lysosome.transport.10., main = "ps1 vs wt: endosome.to.lysosome.transport (n=10)")
barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$endosome.to.lysosome.transport.10., main = "ps1.iso vs ps1: endosome.to.lysosome.transport (n=10)")
#dev.print(pdf, "GSEA/Rplot_barcode3_UpwithTFEB.OE_ps1VSwt_ps1isoVSps1.pdf", height=6, width=12)


barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$Lysoplex.All.Autophagy.23., main = "ps1.iso vs ps1: Lysoplex: All Autophagy (n=23)")
dev.print(pdf, "GSEA/Rplot_barcode3_lysoplexAuto_ps1isoVSps1.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$All.Autophagic.genes.147., main = "ps1.iso vs ps1: All Autophagic Genes (n=147)")
dev.print(pdf, "GSEA/Rplot_barcode3_allAutophagy_ps1isoVSps1.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$Autophagy.17., main = "ps1.iso vs ps1: Autophagy (n=17)")
dev.print(pdf, "GSEA/Rplot_barcode3_Autophagy_ps1isoVSps1.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$Autophagy.17., main = "ps1.iso vs ps1: Autophagy (n=17)")
dev.print(pdf, "GSEA/Rplot_barcode3_Autophagy_ps1isoVSps1.pdf", height=6, width=8)

## Checking the difference in GSEA results ps1VSwt versus ps1isoVSwt
head(mroast3.gs100_ps1VSwt)
head(mroast3.gs100_ps1isoVSwt)

ord.mroast3.gs100_ps1VSwt <- mroast3.gs100_ps1VSwt[match(rownames(mroast3.gs100_ps1isoVSwt), rownames(mroast3.gs100_ps1VSwt)), ]
isTRUE(all.equal(rownames(ord.mroast3.gs100_ps1VSwt), rownames(mroast3.gs100_ps1isoVSwt))) # TRUE
###
ord3.1.mroast3.gs100_ps1VSwt <- mroast3.gs100_ps1VSwt[match(rownames(mroast3.gs100_ps1isoVSps1), rownames(mroast3.gs100_ps1VSwt)), ]
isTRUE(all.equal(rownames(ord3.1.mroast3.gs100_ps1VSwt), rownames(mroast3.gs100_ps1isoVSps1))) # TRUE

mroast3.diff <- mroast3.gs100_ps1isoVSwt[mroast3.gs100_ps1isoVSwt$Direction != ord.mroast3.gs100_ps1VSwt$Direction, ] # Direction: Change; FDR: Sig
mroast3.diff

mroast3.diff2 <- mroast3.gs100_ps1isoVSwt[mroast3.gs100_ps1isoVSwt$FDR < 0.1 & ord.mroast3.gs100_ps1VSwt$FDR >= 0.1, ] # Direction: NA; FDR: Sig
mroast3.diff2     ## Consider FDR < 0.05 as significant although 0.1 has been used as cutoff just to reval names of gene sets.

# (03-014-2019) Same direction but FDR sig with iso (Less Robust iso effect)
mroast3.diff3 <- mroast3.gs100_ps1isoVSwt[ord.mroast3.gs100_ps1VSwt$FDR >= 0.1 & ord.mroast3.gs100_ps1VSwt$Direction == mroast3.gs100_ps1isoVSwt$Direction
                                          & mroast3.gs100_ps1isoVSwt$FDR < 0.1, ]
mroast3.diff3
###
mroast3.diff3.1 <- mroast3.gs100_ps1isoVSps1[ord3.1.mroast3.gs100_ps1VSwt$FDR >= 0.1 & ord3.1.mroast3.gs100_ps1VSwt$Direction == mroast3.gs100_ps1isoVSps1$Direction
                                          & mroast3.gs100_ps1isoVSps1$FDR < 0.1, ]
mroast3.diff3.1

# ps1VSwt: NOT sig. BUT Direction changes significantly with iso (Moderately robust iso effect)
mroast3.diff4 <- mroast3.gs100_ps1isoVSwt[ord.mroast3.gs100_ps1VSwt$FDR >= 0.1 & ord.mroast3.gs100_ps1VSwt$Direction != mroast3.gs100_ps1isoVSwt$Direction
                                          & mroast3.gs100_ps1isoVSwt$FDR < 0.1, ]
mroast3.diff4
###
mroast3.diff4.1 <- mroast3.gs100_ps1isoVSps1[ord3.1.mroast3.gs100_ps1VSwt$FDR >= 0.1 & ord3.1.mroast3.gs100_ps1VSwt$Direction != mroast3.gs100_ps1isoVSps1$Direction
                                          & mroast3.gs100_ps1isoVSps1$FDR < 0.1, ]
mroast3.diff4.1

# ps1VSwt: sig. AND Direction changes significantly with iso (MOST ROBUST iso effect)
mroast3.diff5 <- mroast3.gs100_ps1isoVSwt[ord.mroast3.gs100_ps1VSwt$FDR < 0.1 & ord.mroast3.gs100_ps1VSwt$Direction != mroast3.gs100_ps1isoVSwt$Direction
                                          & mroast3.gs100_ps1isoVSwt$FDR < 0.1, ]
mroast3.diff5
###
mroast3.diff5.1 <- mroast3.gs100_ps1isoVSps1[ord3.1.mroast3.gs100_ps1VSwt$FDR < 0.1 & ord3.1.mroast3.gs100_ps1VSwt$Direction != mroast3.gs100_ps1isoVSps1$Direction
                                          & mroast3.gs100_ps1isoVSps1$FDR < 0.1, ]
mroast3.diff5.1

### NEW: 03-18-2019: psVSwt: Any direction and sig.; ps1isoVSwt: Any direction and NOT sig.
mroast.diff6 <- mroast3.gs100_ps1isoVSwt[ord.mroast3.gs100_ps1VSwt$FDR < 0.05 & mroast3.gs100_ps1isoVSwt$FDR >= 0.05, ]
mroast.diff6

## Checking which genes changed direction significantly in DGE matrices of ps1VSwt and ps1isoVSwt
# Making order of genes same in both dataframes
ord.res3.ps1VSwt <- res3.ps1VSwt[match(res3.ps1isoVSwt$genes, res3.ps1VSwt$genes), ]
head(ord.res3.ps1VSwt)
head(res3.ps1isoVSwt)
tail(ord.res3.ps1VSwt)
tail(res3.ps1isoVSwt)

sigChange1 <- res3.ps1isoVSwt[sign(res3.ps1isoVSwt$logFC) != sign(ord.res3.ps1VSwt$logFC) & res3.ps1isoVSwt$adj.P.Val <= 0.05, ]
dim(sigChange1) # 169   7
head(sigChange1)
write.table(sigChange1, file = "sigDirChangeWithIsoTreatment.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

sigUP1 <- res3.ps1isoVSwt[sign(res3.ps1isoVSwt$logFC) != sign(ord.res3.ps1VSwt$logFC) & sign(res3.ps1isoVSwt$logFC) > 0
                          & res3.ps1isoVSwt$adj.P.Val <= 0.05, ]
dim(sigUP1) # 81   7
head(sigUP1)
tail(sigUP1)
write.table(sigUP1, file = "sigUPWithIsoTreatment.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

sigDOWN1 <- res3.ps1isoVSwt[sign(res3.ps1isoVSwt$logFC) != sign(ord.res3.ps1VSwt$logFC) & sign(res3.ps1isoVSwt$logFC) < 0
                            & res3.ps1isoVSwt$adj.P.Val <= 0.05, ]
dim(sigDOWN1) # 88   7
head(sigDOWN1)
tail(sigDOWN1)
write.table(sigDOWN1, file = "sigDOWNWithIsoTreatment.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)


## Camera gene set enrichment analysis with MSigDB i.e. Mm.c2
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata"))
class(Mm.c2) # "list"


# This will load Mm.c2, which is a list of gene sets, each a vector of Entrez Ids.
idx.all.Mm.c2 <- ids2indices(Mm.c2, id = rownames(v3.1))

cam3.1.Mm.c2_ps1VSwt <- camera(v3.1, idx.all.Mm.c2, design, contrast = contr.mat[,2])
head(cam3.1.Mm.c2_ps1VSwt, 10)
write.table(cam3.1.Mm.c2_ps1VSwt, file = "GSEA/cam3.1.Mm.c2_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

cam3.1.Mm.c2_ps1isoVSwt <- camera(v3.1, idx.all.Mm.c2, design, contrast = contr.mat[,4])
head(cam3.1.Mm.c2_ps1isoVSwt, 10)
write.table(cam3.1.Mm.c2_ps1isoVSwt, file = "GSEA/cam3.1.Mm.c2_ps1isoVSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

cam3.1.Mm.c2_ps1isoVSps1 <- camera(v3.1, idx.all.Mm.c2, design, contrast = contr.mat[,1])
head(cam3.1.Mm.c2_ps1isoVSps1, 10)
write.table(cam3.1.Mm.c2_ps1isoVSps1, file = "GSEA/cam3.1.Mm.c2_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


## Some heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$Lysosomal.proteolysis.23., c(7:9,1:3,4:6)], 
          scale = "row", dendrogram = "none", Colv = "dendrogram",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$Lysosomal.proteolysis.23.],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% geneSets110_Mm$Lysosomal.proteolysis.23., c(2,1)], 
          dendrogram = "none", Colv = "dendrogram", scale = "column",
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% geneSets110_Mm$Lysosomal.proteolysis.23.],
          col = mycol3, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)


### Checking a few MSigDb gene sets of interest
# 1. Chloride channel gene sets

## MUST READ https://stackoverflow.com/a/18922914 to read a table with different # of columns.

gs.clcn.msigdb_Hs <- read.table(file = "GSEA/genesets_ChlorideChannel_MSigDb.gmt", sep = "\t", header = F, fill = T,
                                 stringsAsFactors = F)


head(gs.clcn.msigdb_Hs)

gs.clcn.msigdb_Hs <- as.data.frame(t(gs.clcn.msigdb_Hs))
head(gs.clcn.msigdb_Hs)

colnames(gs.clcn.msigdb_Hs) # "V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8"
colnames(gs.clcn.msigdb_Hs) <- as.character(unlist(gs.clcn.msigdb_Hs[1,]))

gs.clcn.msigdb_Hs <- gs.clcn.msigdb_Hs[-c(1,2), ]
rownames(gs.clcn.msigdb_Hs) <- NULL
head(gs.clcn.msigdb_Hs)
tail(gs.clcn.msigdb_Hs)
dim(gs.clcn.msigdb_Hs) # 99  8

### IMP: We ned to convert all characters to prop case
gs.clcn.msigdb_Mm <- data.frame(lapply(gs.clcn.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.clcn.msigdb_Mm)

## convert into ids
idx.all.clcn.msigdb <- ids2indices(gs.clcn.msigdb_Mm, id = rownames(v3))
str(idx.all.clcn.msigdb)

## mroast
mroast3.clcn.msigdb_ps1VSwt <- mroast(v3, idx.all.clcn.msigdb, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.clcn.msigdb_ps1VSwt
write.table(mroast3.clcn.msigdb_ps1VSwt, file = "GSEA/mroast3.clcn.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.clcn.msigdb_ps1isoVSps1 <- mroast(v3, idx.all.clcn.msigdb, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.clcn.msigdb_ps1isoVSps1
write.table(mroast3.clcn.msigdb_ps1isoVSps1, file = "GSEA/mroast3.clcn.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.clcn.msigdb_ps1isoVSwt <- mroast(v3, idx.all.clcn.msigdb, design, contrast = contr.mat[,4], nrot = 100000)
mroast3.clcn.msigdb_ps1isoVSwt

par(mfrow=c(1,2))
barcodeplot(efit3.contr$t[,2], index = idx.all.clcn.msigdb$GO_CHLORIDE_TRANSPORT, main = "", 
            labels = c("UP in PS1KO", "DOWN in PS1KO")) # 	GO:0006821
title(main = "CHLORIDE TRANSPORT (GO:0006821) \nPS1KO vs WT")
barcodeplot(efit3.contr$t[,1], index = idx.all.clcn.msigdb$GO_CHLORIDE_TRANSPORT, main = "", 
            labels = c("UP in PS1KO+ISO", "DOWN in PS1KO+ISO")) # 	GO:0006821
title(main = "CHLORIDE TRANSPORT (GO:0006821) \nPS1KO+ISO vs PS1KO")
dev.print(pdf, "GSEA/Rplot_barcode3_Cl.Trans_ps1VSwt.ps1isoVSps1.pdf", height=6, width=12)  # 	GO:0006821
## Change in orientation
par(mfrow=c(2,1))
barcodeplot(efit3.contr$t[,2], index = idx.all.clcn.msigdb$GO_CHLORIDE_TRANSPORT, main = "ps1 vs wt: GO_CHLORIDE_TRANSPORT (n=46)", 
            labels = c("Down in PS1KO", "Up in PS1KO"))
barcodeplot(efit3.contr$t[,1], index = idx.all.clcn.msigdb$GO_CHLORIDE_TRANSPORT, main = "ps1.iso vs ps1: GO_CHLORIDE_TRANSPORT (n=46)", 
            labels = c("Down in PS1KO+ISO", "Up in PS1KO+ISO"))


## Heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.clcn.msigdb_Mm$BIOCARTA_CFTR_PATHWAY, c(7:9,1:3,4:6)], 
          scale = "row", Colv = T, dendrogram = "column",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.clcn.msigdb_Mm$BIOCARTA_CFTR_PATHWAY],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))



#2. PKA gene sets
gs.pka.msigdb_Hs <- read.table(file = "GSEA/genesets_PKA_MSigDb.gmt", sep = "\t", header = F, fill = T,
                                stringsAsFactors = F)


head(gs.pka.msigdb_Hs)

gs.pka.msigdb_Hs <- as.data.frame(t(gs.pka.msigdb_Hs))
head(gs.pka.msigdb_Hs)

colnames(gs.pka.msigdb_Hs) # "V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8"
colnames(gs.pka.msigdb_Hs) <- as.character(unlist(gs.pka.msigdb_Hs[1,]))

gs.pka.msigdb_Hs <- gs.pka.msigdb_Hs[-c(1,2), ]
rownames(gs.pka.msigdb_Hs) <- NULL
head(gs.pka.msigdb_Hs)
tail(gs.pka.msigdb_Hs)
dim(gs.pka.msigdb_Hs) # 23  4

### IMP: We ned to convert all characters to prop case
gs.pka.msigdb_Mm <- data.frame(lapply(gs.pka.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.pka.msigdb_Mm)

## convert into ids
idx.all.pka.msigdb <- ids2indices(gs.pka.msigdb_Mm, id = rownames(v3))
str(idx.all.pka.msigdb)

## mroast
mroast3.pka.msigdb_ps1VSwt <- mroast(v3, idx.all.pka.msigdb, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.pka.msigdb_ps1VSwt
write.table(mroast3.pka.msigdb_ps1VSwt, file = "GSEA/mroast3.pka.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.pka.msigdb_ps1isoVSps1 <- mroast(v3, idx.all.pka.msigdb, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.pka.msigdb_ps1isoVSps1
write.table(mroast3.pka.msigdb_ps1isoVSps1, file = "GSEA/mroast3.pka.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.pka.msigdb_ps1isoVSwt <- mroast(v3, idx.all.pka.msigdb, design, contrast = contr.mat[,4], nrot = 100000)
mroast3.pka.msigdb_ps1isoVSwt
#write.table(mroast3.pka.msigdb_ps1isoVSps1, file = "GSEA/mroast3.pka.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

res3.ps1VSwt[res3.ps1VSwt$genes %in% gs.pka.msigdb_Mm$BIOCARTA_AGPCR_PATHWAY, ]
res3.ps1isoVSps1[res3.ps1isoVSps1$genes %in% gs.pka.msigdb_Mm$BIOCARTA_AGPCR_PATHWAY, ]
res3.ps1isoVSwt[res3.ps1isoVSwt$genes %in% gs.pka.msigdb_Mm$BIOCARTA_AGPCR_PATHWAY, ]

## Heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.pka.msigdb_Mm$BIOCARTA_AGPCR_PATHWAY, c(7:9,1:3,4:6)], 
          scale = "row", Colv = T, dendrogram = "column",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.pka.msigdb_Mm$BIOCARTA_AGPCR_PATHWAY],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.pka.msigdb_Mm$GO_REGULATION_OF_PROTEIN_KINASE_A_SIGNALING, c(7:9,1:3,4:6)], 
          scale = "row", Colv = T, dendrogram = "column",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.pka.msigdb_Mm$GO_REGULATION_OF_PROTEIN_KINASE_A_SIGNALING],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))





#3. SASSON_PKA genes
gs.sasson_pka.msigdb_Hs <- read.table(file = "GSEA/genesets_SASSON_PKA_MSigDb.gmt", sep = "\t", header = F, fill = T,
                               stringsAsFactors = F)

gs.sasson_pka.msigdb_Hs <- as.data.frame(t(gs.sasson_pka.msigdb_Hs))
head(gs.sasson_pka.msigdb_Hs)

colnames(gs.sasson_pka.msigdb_Hs) # "V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8"
colnames(gs.sasson_pka.msigdb_Hs) <- as.character(unlist(gs.sasson_pka.msigdb_Hs[1,]))

gs.sasson_pka.msigdb_Hs <- gs.sasson_pka.msigdb_Hs[-c(1,2), ]
rownames(gs.sasson_pka.msigdb_Hs) <- NULL
head(gs.sasson_pka.msigdb_Hs)
tail(gs.sasson_pka.msigdb_Hs)
dim(gs.sasson_pka.msigdb_Hs) # 91  5

### IMP: We ned to convert all characters to prop case
gs.sasson_pka.msigdb_Mm <- data.frame(lapply(gs.sasson_pka.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.sasson_pka.msigdb_Mm)

## convert into ids
idx.all.sasson_pka.msigdb <- ids2indices(gs.sasson_pka.msigdb_Mm, id = rownames(v3))
str(idx.all.sasson_pka.msigdb)

## mroast
mroast3.sasson_pka.msigdb_ps1VSwt <- mroast(v3, idx.all.sasson_pka.msigdb, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.sasson_pka.msigdb_ps1VSwt
write.table(mroast3.sasson_pka.msigdb_ps1VSwt, file = "GSEA/mroast3.sasson_pka.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.sasson_pka.msigdb_ps1isoVSps1 <- mroast(v3, idx.all.sasson_pka.msigdb, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.sasson_pka.msigdb_ps1isoVSps1
write.table(mroast3.sasson_pka.msigdb_ps1isoVSps1, file = "GSEA/mroast3.sasson_pka.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

## Heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.sasson_pka.msigdb_Mm$SASSON_RESPONSE_TO_FORSKOLIN_UP, c(7:9,1:3,4:6)], 
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.sasson_pka.msigdb_Mm$SASSON_RESPONSE_TO_FORSKOLIN_UP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.sasson_pka.msigdb_Mm$SASSON_RESPONSE_TO_FORSKOLIN_UP, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", #Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.sasson_pka.msigdb_Mm$SASSON_RESPONSE_TO_FORSKOLIN_UP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)

## Barcode plots
barcodeplot(efit3.contr$t[,2], index = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_UP, index2 = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_DN,
            main = "PS1KO vs WT: SASSON_RESPONSE_TO_FORSKOLIN")
dev.print(pdf, "GSEA/Rplot_barcode3_Sasson.Forsk_ps1VSwt.pdf", height=6, width=8)
# for MS
barcodeplot(efit3.contr$t[,2], index = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_UP, index2 = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_DN,
            main = "PS1KO vs WT: SASSON_RESPONSE_TO_FORSKOLIN", labels = c("", ""))

barcodeplot(efit3.contr$t[,1], index = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_UP, index2 = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_DN,
            main = "PS1KO+ISO vs PS1KO: SASSON_RESPONSE_TO_FORSKOLIN")
dev.print(pdf, "GSEA/Rplot_barcode3_Sasson.Forsk_ps1isoVSps1.pdf", height=6, width=8)
# for MS
barcodeplot(efit3.contr$t[,1], index = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_UP, index2 = idx.all.sasson_pka.msigdb$SASSON_RESPONSE_TO_FORSKOLIN_DN,
            main = "PS1KO+ISO vs PS1KO: SASSON_RESPONSE_TO_FORSKOLIN", labels = c("", ""))

# genas means Genuine Association of Gene Expression Profiles
test1 <- genas(efit3.contr[row.names(efit3.contr) %in% gs.sasson_pka.msigdb_Mm$SASSON_RESPONSE_TO_FORSKOLIN_UP, ], coef = c(2,1), subset = "F", plot = T)
test2 <- genas(efit3.contr[row.names(efit3.contr) %in% gs.sasson_pka.msigdb_Mm$SASSON_RESPONSE_TO_FORSKOLIN_DN, ], coef = c(2,1), subset = "F", plot = T)


#4. Adenylate Cyclase gene sets

## MUST READ https://stackoverflow.com/a/18922914 to read a table with different # of columns.
count.fields("GSEA/genesets_adenylateCyclase_MSigDb.gmt", sep = "\t") # To count maximum # of fields
# 27  66  14  75  70 147  92  90  89  93

# Now read in table
gs.ac.msigdb_Hs <- read.table(file = "GSEA/genesets_adenylateCyclase_MSigDb.gmt", sep = "\t", header = F, fill = T,
                              col.names = paste0("V", seq_len(147)), stringsAsFactors = F)
head(gs.ac.msigdb_Hs)

gs.ac.msigdb_Hs <- as.data.frame(t(gs.ac.msigdb_Hs))
head(gs.ac.msigdb_Hs)
tail(gs.ac.msigdb_Hs)

colnames(gs.ac.msigdb_Hs) # "V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8" "V9" "V10"
colnames(gs.ac.msigdb_Hs) <- as.character(unlist(gs.ac.msigdb_Hs[1,]))

gs.ac.msigdb_Hs <- gs.ac.msigdb_Hs[-c(1,2), ]
rownames(gs.ac.msigdb_Hs) <- NULL
head(gs.ac.msigdb_Hs)
tail(gs.ac.msigdb_Hs)
dim(gs.ac.msigdb_Hs) # 145  10

### IMP: We ned to convert all characters to prop case
gs.ac.msigdb_Mm <- data.frame(lapply(gs.ac.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.ac.msigdb_Mm)

## convert into ids
idx.all.ac.msigdb <- ids2indices(gs.ac.msigdb_Mm, id = rownames(v3))
str(idx.all.ac.msigdb)

## mroast
mroast3.ac.msigdb_ps1VSwt <- mroast(v3, idx.all.ac.msigdb, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.ac.msigdb_ps1VSwt
write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.ac.msigdb_ps1isoVSps1 <- mroast(v3, idx.all.ac.msigdb, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.ac.msigdb_ps1isoVSps1
write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.ac.msigdb_ps1isoVSwt <- mroast(v3, idx.all.ac.msigdb, design, contrast = contr.mat[,4], nrot = 100000)
mroast3.ac.msigdb_ps1isoVSwt

## Heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER, c(7:9,1:3,4:6)], 
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(-log10(efit3.contr$p.value[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER, c(2,1)]), 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)
###
#TESTING
df_test <- merge(res3.ps1VSwt, res3.ps1isoVSps1, by="genes")
head(df_test)
heatmap.2(as.matrix(-log10(df_test[df_test$genes %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER, c(6,12)])), 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = df_test$genes[df_test$genes %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)


heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY, c(7:9,1:3,4:6)], 
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)
###
heatmap.2(-log10(efit3.contr$p.value[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY, c(2,1)]), 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)
###
heatmap.2(as.matrix(-log10(df_test[df_test$genes %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY, c(6,12)])), 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = df_test$genes[df_test$genes %in% gs.ac.msigdb_Mm$GO_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)

#5. mcalpain_biocarta_msigdb
gs.calp.msigdb_Hs <- read.table(file = "GSEA/geneset_mcalpain.BioCarta_MSigDb.gmt", sep = "\t", header = F, fill = T,
                              stringsAsFactors = F)
head(gs.calp.msigdb_Hs)

gs.calp.msigdb_Hs <- as.data.frame(t(gs.calp.msigdb_Hs))
head(gs.calp.msigdb_Hs)
tail(gs.calp.msigdb_Hs)

colnames(gs.calp.msigdb_Hs) # "V1"
colnames(gs.calp.msigdb_Hs) <- as.character(unlist(gs.calp.msigdb_Hs[1,]))

gs.calp.msigdb_Hs <- gs.calp.msigdb_Hs[-c(1,2), ]
rownames(gs.calp.msigdb_Hs) <- NULL
head(gs.calp.msigdb_Hs)
tail(gs.calp.msigdb_Hs)
dim(gs.calp.msigdb_Hs) # NULL

gs.calp.msigdb_Hs <- as.character(gs.calp.msigdb_Hs)
head(gs.calp.msigdb_Hs)
length(gs.calp.msigdb_Hs) # 25

### IMP: We ned to convert all characters to prop case
gs.calp.msigdb_Mm <- as.character(lapply(gs.calp.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}))

head(gs.calp.msigdb_Mm)

## convert into ids
idx.calp.msigdb <- ids2indices(gs.calp.msigdb_Mm, id = rownames(v3))
str(idx.calp.msigdb)

## roast
roast3.calp.msigdb_ps1VSwt <- roast(v3, idx.calp.msigdb, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.calp.msigdb_ps1VSwt
#write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.calp.msigdb_ps1isoVSps1 <- roast(v3, idx.calp.msigdb, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.calp.msigdb_ps1isoVSps1
#write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)


#6. MITF target genes
gs.mitf.Hman_Hs <- loadWorkbook("GSEA/geneset_mitf.targets_Hartman.xlsx")
gs.mitf.Hman_Hs <- readWorksheet(gs.mitf.Hman_Hs, sheet = 1, header = T)
head(gs.mitf.Hman_Hs)
dim(gs.mitf.Hman_Hs) # 78  1

#gs.calp.msigdb_Hs <- as.data.frame(t(gs.calp.msigdb_Hs))
#head(gs.calp.msigdb_Hs)
#tail(gs.calp.msigdb_Hs)

#colnames(gs.calp.msigdb_Hs) # "V1"
#colnames(gs.calp.msigdb_Hs) <- as.character(unlist(gs.calp.msigdb_Hs[1,]))

#gs.calp.msigdb_Hs <- gs.calp.msigdb_Hs[-c(1,2), ]
#rownames(gs.calp.msigdb_Hs) <- NULL
#head(gs.calp.msigdb_Hs)
#tail(gs.calp.msigdb_Hs)
#dim(gs.calp.msigdb_Hs) # NULL

#gs.calp.msigdb_Hs <- as.character(gs.calp.msigdb_Hs)
#head(gs.calp.msigdb_Hs)
#length(gs.calp.msigdb_Hs) # 25

### IMP: We ned to convert all characters to prop case
#gs.mitf.Hman_Mm <- as.data.frame(lapply(gs.mitf.Hman_Hs, function(v) {
#  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
#}), stringsAsFactors = F)

gs.mitf.Hman_Mm <- as.character(lapply(gs.mitf.Hman_Hs$SYMBOLS, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}))

head(gs.mitf.Hman_Mm)
tail(gs.mitf.Hman_Mm)

gs.mitf.Hman_Mm <- gsub("(^\\s+)|(\\s+$)", "", gs.mitf.Hman_Mm)


## convert into ids
idx.mitf.Hman <- ids2indices(gs.mitf.Hman_Mm, id = rownames(v3))
str(idx.mitf.Hman)

## roast
roast3.mitf.Hman_ps1VSwt <- roast(v3, idx.mitf.Hman, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.mitf.Hman_ps1VSwt
#write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.mitf.Hman_ps1isoVSps1 <- roast(v3, idx.mitf.Hman, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.mitf.Hman_ps1isoVSps1
#write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.mitf.Hman_ps1isoVSwt <- roast(v3, idx.mitf.Hman, design, contrast = contr.mat[,4], nrot = 1000000)
roast3.mitf.Hman_ps1isoVSwt

## Check the expression values of these genes for both comparisons
res3.ps1VSwt[match(gs.mitf.Hman_Mm, res3.ps1VSwt$genes, nomatch = 0), ]
res3.ps1isoVSps1[match(gs.mitf.Hman_Mm, res3.ps1isoVSps1$genes, nomatch = 0), ]
res3.ps1isoVSwt[match(gs.mitf.Hman_Mm, res3.ps1isoVSwt$genes, nomatch = 0), ]

## Heatmap
heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% gs.mitf.Hman_Mm, c(2,4)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1)

heatmap.2(efit3.contr$p.value[row.names(efit3.contr) %in% gs.mitf.Hman_Mm, c(2,4)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1)

heatmap.2(efit3.contr$lods[row.names(efit3.contr) %in% gs.mitf.Hman_Mm, c(2,4)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1)

heatmap.2(v3$E[row.names(v3$E) %in% gs.mitf.Hman_Mm, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% gs.mitf.Hman_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10)) # BETTER

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Hman_Mm, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.mitf.Hman_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10)) ## Even BETTER
dev.print(pdf, "Rplot_heatmap_mitf.Hman.pdf", width=8, height=10)

heatmap.2(v3$E[row.names(v3$E) %in% gs.mitf.Hman_Mm, c(1:3,4:6)], scale = "row", Rowv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% gs.mitf.Hman_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

## Finding correlations
corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Hman_Mm, c(7:9,1:3,4:6)], adjust = "BH")
cor_mitf.Hman <- corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Hman_Mm, c(7:9,1:3,4:6)], adjust = "BH")
cor_mitf.Hman$r
cor_mitf.Hman$


#7. Curated Mitf targets_Ennen2017 Clinical Cancer Research

gs.mitf.Ennen_Hs <- loadWorkbook("GSEA/Curated Mitf targets_Ennen2017_ClinCancRes.xlsx")
gs.mitf.Ennen_Hs <- readWorksheet(gs.mitf.Ennen_Hs, sheet = 1, header = T)
head(gs.mitf.Ennen_Hs)
dim(gs.mitf.Ennen_Hs) # 77  2

### IMP: We ned to convert all characters to prop case
gs.mitf.Ennen_Mm <- as.data.frame(lapply(gs.mitf.Ennen_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.mitf.Ennen_Mm)
tail(gs.mitf.Ennen_Mm)

#gs.mitf.Hman_Mm <- gsub("(^\\s+)|(\\s+$)", "", gs.mitf.Hman_Mm)


## convert into ids
idx.mitf.Ennen <- ids2indices(gs.mitf.Ennen_Mm, id = rownames(v3))
str(idx.mitf.Ennen)

## roast
roast3.mitf.Ennen_ps1VSwt <- roast(v3, idx.mitf.Ennen, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.mitf.Ennen_ps1VSwt
#write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.mitf.Ennen_ps1isoVSps1 <- roast(v3, idx.mitf.Ennen, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.mitf.Ennen_ps1isoVSps1
#write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.mitf.Ennen_ps1isoVSwt <- roast(v3, idx.mitf.Ennen, design, contrast = contr.mat[,4], nrot = 1000000)
roast3.mitf.Ennen_ps1isoVSwt

## Check the expression values of these genes for both comparisons
res3.ps1VSwt[match(gs.mitf.Ennen_Mm$Enriched.in.MITF.targets, res3.ps1VSwt$genes, nomatch = 0), ]
res3.ps1isoVSwt[match(gs.mitf.Ennen_Mm$Enriched.in.MITF.targets, res3.ps1isoVSwt$genes, nomatch = 0), ]

## Heatmap
heatmap.2(v3$E[row.names(v3$E) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.targets, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.targets],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.targets, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.targets],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10)) ## Even BETTER
dev.print(pdf, "Rplot_heatmap_enriched.mitf.Ennen..pdf", width=8, height=10)

heatmap.2(v3$E[row.names(v3$E) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10)) ## Even BETTER
dev.print(pdf, "Rplot_heatmap_enriched.mitf.Ennen..pdf", width=8, height=10)

## Find correlations:
library(psych)
cor.test(efit3.contr$p.value[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 2],
         efit3.contr$p.value[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 4]) #Cor 0.3342887, t = 2.4317, df = 47, p-value = 0.0189

cor.test(efit3.contr$p.value[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 2],
         efit3.contr$p.value[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 4])

cor.test(v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 7:9],
         v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 1:3]) # cor 0.8124525

cor.test(v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 7:9],
         v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 4:6]) # cor 0.8106483

corr.test(v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 7:9],
         v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 1:3])

corr.test(v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 7:9],
          v3$E[row.names(v3) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 4:6])

corr.test(efit3.contr$lods[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(2,4,1)])
corr.test(efit3.contr$lods[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(2,4,1)])$p

corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 7:9],
          lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 1:3])

corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 7:9],
          lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, 4:6])

corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(7:9,1:3,4:6)])
corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(7:9,1:3,4:6)])$p # Makes sense

corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.targets, c(7:9,1:3,4:6)])
corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.targets, c(7:9,1:3,4:6)])$p # Makes sense

corr.test(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(7:9,1:3,4:6)])

log10.adjP_mitf_LOW.Ennen <- efit3.contr$p.value[row.names(efit3.contr) %in% gs.mitf.Ennen_Mm$Enriched.in.MITF.low.signature.genes, c(2,4,1)]
log10.adjP_mitf_LOW.Ennen <- -log10(log10.adjP_mitf_LOW.Ennen[,1:3])
head(log10.adjP_mitf_LOW.Ennen)

corr.test(log10.adjP_mitf_LOW.Ennen, adjust = "BH")
corr.test(log10.adjP_mitf_LOW.Ennen, adjust = "BH")$p


#8. ER to Golgi Anterograde Transport from GO: 0006888

gs.ErToGolgi.GO_Mm <- read.table(file = "GO/GO_0006888_ER_to_golgi_anterograde_060319.txt", sep = "\t", header = T, fill = T, row.names = NULL,
                                  stringsAsFactors = F)
head(gs.ErToGolgi.GO_Mm)

gs.ErToGolgi.GO_Mm <- as.data.frame(gs.ErToGolgi.GO_Mm[,2], stringsAsFactors = F)
head(gs.ErToGolgi.GO_Mm)
colnames(gs.ErToGolgi.GO_Mm) <- "Symbol"
dim(gs.ErToGolgi.GO_Mm) # 188   1
class(gs.ErToGolgi.GO_Mm$Symbol) # "character"

sum(duplicated(gs.ErToGolgi.GO_Mm$Symbol)) # 65
gs.ErToGolgi.GO_Mm <- gs.ErToGolgi.GO_Mm[!duplicated(gs.ErToGolgi.GO_Mm$Symbol), ]
length(gs.ErToGolgi.GO_Mm) # 123

sum(is.na(gs.ErToGolgi.GO_Mm)) # 0
sum(gs.ErToGolgi.GO_Mm == "") # 0

## convert into ids
idx.ErToGolgi.GO <- ids2indices(gs.ErToGolgi.GO_Mm, id = rownames(v3))
str(idx.ErToGolgi.GO) # 112

## roast
roast3.ErToGolgi.GO_ps1VSwt <- roast(v3, idx.ErToGolgi.GO, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.ErToGolgi.GO_ps1VSwt
#write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.ErToGolgi.GO_ps1isoVSps1 <- roast(v3, idx.ErToGolgi.GO, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.ErToGolgi.GO_ps1isoVSps1
#write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

## Heatmap
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10)) ## Even BETTER
#dev.print(pdf, "Rplot_heatmap_enriched.mitf.Ennen..pdf", width=8, height=10)

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm, c(1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm],
          labCol = c("PS1KO_1","PS1KO_2","PS1KO_3","PS1KO+ISO_1","PS1KO+ISO_2","PS1KO+ISO_3"),
          col = mycol, trace = "none", density.info = "none", srtCol = 0, cexCol = 0.7, adjCol = c(NA,0),
          margins = c(8,6), lhei = c(2,10))
dev.print(pdf, "Rplot_heatmap_ER_to_Golgi.GO_0006888..pdf", width=8, height=10)


heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm, c(1:3,4:6)], scale = "row", dendrogram = "column",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm],
          labCol = c("PS1KO_1","PS1KO_2","PS1KO_3","PS1KO+ISO_1","PS1KO+ISO_2","PS1KO+ISO_3"),
          col = mycol, trace = "none", density.info = "none", srtCol = 0, cexCol = 0.7, adjCol = c(NA,0),
          margins = c(8,6), lhei = c(2,10)) ## USE THIS
dev.print(pdf, "Rplot_heatmap_ER_to_Golgi_withDendrogram.GO_0006888..pdf", width=8, height=10)


heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.ErToGolgi.GO_Mm, c(2,1)], dendrogram = "none", Colv = "dendrogram", scale = "column",
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.ErToGolgi.GO_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
#dev.print(pdf, "Rplot_heatmap_ER_to_Golgi.GO_0006888..pdf", width=8, height=10)

## Practice heatmap
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm, c(1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.ErToGolgi.GO_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(4,22), lhei = c(1,6), lwid = c(2,4), cexCol = 0.85)
dev.print(pdf, "Rplot_heatmap_ER_to_Golgi.GO_0006888_2.pdf", width=8, height=10)

#9. Goldi to Plasm Membrane Protein Transport from GO: 00430001

gs.GolgiToPM.GO_Mm <- read.table(file = "GO/GO_00430001_golgi_to_PM_anterograde_061719.txt", sep = "\t", header = T, fill = T, row.names = NULL,
                                 stringsAsFactors = F)
head(gs.GolgiToPM.GO_Mm)

gs.GolgiToPM.GO_Mm <- as.data.frame(gs.GolgiToPM.GO_Mm[,2], stringsAsFactors = F)
head(gs.GolgiToPM.GO_Mm)
colnames(gs.GolgiToPM.GO_Mm) <- "Symbol"
dim(gs.GolgiToPM.GO_Mm) # 49   1
class(gs.GolgiToPM.GO_Mm$Symbol) # "character"

sum(duplicated(gs.GolgiToPM.GO_Mm$Symbol)) # 10
gs.GolgiToPM.GO_Mm <- gs.GolgiToPM.GO_Mm[!duplicated(gs.GolgiToPM.GO_Mm$Symbol), ]
length(gs.GolgiToPM.GO_Mm) # 39

sum(is.na(gs.GolgiToPM.GO_Mm)) # 0
sum(gs.GolgiToPM.GO_Mm == "") # 0

## convert into ids
idx.GolgiToPM.GO <- ids2indices(gs.GolgiToPM.GO_Mm, id = rownames(v3))
str(idx.GolgiToPM.GO ) # 1:35

## roast
roast3.GolgiToPM.GO_ps1VSwt <- roast(v3, idx.GolgiToPM.GO, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.GolgiToPM.GO_ps1VSwt
#write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.GolgiToPM.GO_ps1isoVSps1 <- roast(v3, idx.GolgiToPM.GO, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.GolgiToPM.GO_ps1isoVSps1
#write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

## Heatmap
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.GolgiToPM.GO_Mm, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.GolgiToPM.GO_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.GolgiToPM.GO_Mm, c(1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.GolgiToPM.GO_Mm],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

#10. Protein Targetting to Lysosome from MSigDb (GO:0006622)
gs.ProtToLyso.GO <- read.table(file = "GO/geneset_protTargetingToLyso_GO0006622.gmt", sep = "\t", header = F, fill = T,
                               stringsAsFactors = F)


head(gs.ProtToLyso.GO)

gs.ProtToLyso.GO <- as.data.frame(t(gs.ProtToLyso.GO))
head(gs.ProtToLyso.GO)

colnames(gs.ProtToLyso.GO) # "V1" 
colnames(gs.ProtToLyso.GO) <- as.character(unlist(gs.ProtToLyso.GO[1,]))
head(gs.ProtToLyso.GO)

gs.ProtToLyso.GO <- gs.ProtToLyso.GO[-c(1,2), ]
rownames(gs.ProtToLyso.GO) <- NULL
head(gs.ProtToLyso.GO)
tail(gs.ProtToLyso.GO)
gs.ProtToLyso.GO <- as.character(gs.ProtToLyso.GO)
length(gs.ProtToLyso.GO) # 15

### IMP: We ned to convert all characters to prop case
gs.ProtToLyso.GO_Mm <- as.character(lapply(gs.ProtToLyso.GO, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.ProtToLyso.GO_Mm)


## convert into ids
idx.all.protToLyso.GO <- ids2indices(gs.ProtToLyso.GO_Mm, id = rownames(v3))
str(idx.all.protToLyso.GO)

## mroast
roast3.protToLyso.GO_ps1VSwt <- roast(v3, idx.all.protToLyso.GO, design, contrast = contr.mat[,2], nrot = 100000)
roast3.protToLyso.GO_ps1VSwt
#write.table(mroast3.pka.msigdb_ps1VSwt, file = "GSEA/mroast3.pka.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

roast3.protToLyso.GO_ps1isoVSps1 <- roast(v3, idx.all.protToLyso.GO, design, contrast = contr.mat[,1], nrot = 100000)
roast3.protToLyso.GO_ps1isoVSps1
#write.table(mroast3.pka.msigdb_ps1isoVSps1, file = "GSEA/mroast3.pka.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)


#11. Secretory pathways
gs.secret.msigdb_Hs <- loadWorkbook("GSEA/genesets_secretory_mSigDB.xlsx")
gs.secret.msigdb_Hs <- readWorksheet(gs.secret.msigdb_Hs, header = T, sheet = 1)

head(gs.secret.msigdb_Hs)


colnames(gs.secret.msigdb_Hs) # [1] "SECRETORY_PATHWAY_GO.0045045"                 "HALLMARK_PROTEIN_SECRETION"                  
                              # [3] "SECRETION_BY_CELL"                            "GO_PROTON_TRANSPORTING_V_TYPE_ATPASE_COMPLEX"
gs.secret.msigdb_Hs <- gs.secret.msigdb_Hs[-1, ]
head(gs.secret.msigdb_Hs)
tail(gs.secret.msigdb_Hs)
dim(gs.secret.msigdb_Hs) # 116   4


### IMP: We ned to convert all characters to prop case
gs.secret.msigdb_Mm <- data.frame(lapply(gs.secret.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.secret.msigdb_Mm)
tail(gs.secret.msigdb_Mm)

## convert into ids
idx.all.gs.secret.msigdb <- ids2indices(gs.secret.msigdb_Mm, id = rownames(v3))
str(idx.all.gs.secret.msigdb)

## mroast
mroast3.secret.msigdb_ps1VSwt <- mroast(v3, idx.all.gs.secret.msigdb, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.secret.msigdb_ps1VSwt
write.table(mroast3.secret.msigdb_ps1VSwt, file = "GSEA/mroast3.secretory.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.secret.msigdb_ps1isoVSps1 <- mroast(v3, idx.all.gs.secret.msigdb, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.secret.msigdb_ps1isoVSps1
write.table(mroast3.secret.msigdb_ps1isoVSps1, file = "GSEA/mroast3.secretory.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.secret.msigdb_ps1isoVSwt <- mroast(v3, idx.all.gs.secret.msigdb, design, contrast = contr.mat[,4], nrot = 100000)
mroast3.secret.msigdb_ps1isoVSwt
#write.table(mroast3.secret.msigdb_ps1isoVSwt, file = "GSEA/mroast3.secretory.msigdb_ps1isoVSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

mroast3.secret.msigdb_ps1isoVSwtiso <- mroast(v3, idx.all.gs.secret.msigdb, design, contrast = contr.mat[,5], nrot = 100000)
mroast3.secret.msigdb_ps1isoVSwtiso


## Heatmap
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$HALLMARK_PROTEIN_SECRETION, c(7:9,1:3,4:6)], 
          dendrogram = "column", scale = "row", #Colv = "dendrogram"
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$HALLMARK_PROTEIN_SECRETION],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$HALLMARK_PROTEIN_SECRETION, c(1:3,4:6)], 
          scale = "row", Colv = TRUE, dendrogram = "column", Rowv = TRUE, # ColV and RowV are TRUE by default. 
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$HALLMARK_PROTEIN_SECRETION],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.secret.msigdb_Mm$HALLMARK_PROTEIN_SECRETION, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", Rowv = T, # ColV and RowV are TRUE by default. 
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.secret.msigdb_Mm$HALLMARK_PROTEIN_SECRETION],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)



heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$SECRETORY_PATHWAY_GO.0045045, c(1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$SECRETORY_PATHWAY_GO.0045045],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.secret.msigdb_Mm$SECRETORY_PATHWAY_GO.0045045, c(2,1)], dendrogram = "none",
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.secret.msigdb_Mm$SECRETORY_PATHWAY_GO.0045045],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))


heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$SECRETION_BY_CELL, c(1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.secret.msigdb_Mm$SECRETION_BY_CELL],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

## Barcode plot
barcodeplot(efit3.contr$t[,2], index = idx.all.gs.secret.msigdb$HALLMARK_PROTEIN_SECRETION,
            main = "ps1 vs wt: HALLMARK_PROTEIN_SECRETION")
dev.print(pdf, "GSEA/Rplot_barcode3_hallmark_protSecret_ps1VSwt.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,1], index = idx.all.gs.secret.msigdb$HALLMARK_PROTEIN_SECRETION,
            main = "ps1iso vs ps1: HALLMARK_PROTEIN_SECRETION")
dev.print(pdf, "GSEA/Rplot_barcode3_hallmark_protSecret_ps1isoVSps1.pdf", height=6, width=8)


#12. cAMP pathwys from MSigDb
gs.cAMP.msigdb_Hs <- read.table("GSEA/genesets_cAMP_MSigDb.gmx", sep = "\t", 
                                header = T, stringsAsFactors = F)
head(gs.cAMP.msigdb_Hs)
tail(gs.cAMP.msigdb_Hs)

colnames(gs.cAMP.msigdb_Hs) # [1] "CAMP_UP.V1_DN"                                    "CAMP_UP.V1_UP"                                    "GO_CAMP_BIOSYNTHETIC_PROCESS"                    
# [4] "GO_CAMP_METABOLIC_PROCESS"                        "GO_CELLULAR_RESPONSE_TO_CAMP"                     "GO_POSITIVE_REGULATION_OF_CAMP_METABOLIC_PROCESS"
# [7] "GO_REGULATION_OF_CAMP_METABOLIC_PROCESS"          "GO_RESPONSE_TO_CAMP"                              "X"    

gs.cAMP.msigdb_Hs <- gs.cAMP.msigdb_Hs[-1, -9]
head(gs.cAMP.msigdb_Hs)
tail(gs.cAMP.msigdb_Hs)
dim(gs.cAMP.msigdb_Hs) # 200   8


### IMP: We ned to convert all characters to prop case
gs.cAMP.msigdb_Mm <- data.frame(lapply(gs.cAMP.msigdb_Hs, function(v) {
  ifelse(!is.na(v), paste0(toupper(substr(v, 1, 1)), tolower(substr(v, 2, 20))), v)
}), stringsAsFactors = F)

head(gs.cAMP.msigdb_Mm)
tail(gs.cAMP.msigdb_Mm)

## convert into ids
idx.all.cAMP.msigdb <- ids2indices(gs.cAMP.msigdb_Mm, id = rownames(v3))
str(idx.all.cAMP.msigdb)

## mroast
mroast3.cAMP.msigdb_ps1VSwt <- mroast(v3, idx.all.cAMP.msigdb, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.cAMP.msigdb_ps1VSwt
write.table(mroast3.cAMP.msigdb_ps1VSwt, file = "GSEA/mroast3.cAMP.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.cAMP.msigdb_ps1isoVSps1 <- mroast(v3, idx.all.cAMP.msigdb, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.cAMP.msigdb_ps1isoVSps1
write.table(mroast3.cAMP.msigdb_ps1isoVSps1, file = "GSEA/mroast3.cAMP.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.cAMP.msigdb_ps1isoVSwt <- mroast(v3, idx.all.cAMP.msigdb, design, contrast = contr.mat[,4], nrot = 100000)
mroast3.cAMP.msigdb_ps1isoVSwt

## Heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.cAMP.msigdb_Mm$GO_CELLULAR_RESPONSE_TO_CAMP, c(7:9,1:3,4:6)], 
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.cAMP.msigdb_Mm$GO_CELLULAR_RESPONSE_TO_CAMP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.cAMP.msigdb_Mm$GO_CELLULAR_RESPONSE_TO_CAMP, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.cAMP.msigdb_Mm$GO_CELLULAR_RESPONSE_TO_CAMP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)
###
#heatmap.2(-log10(efit3.contr$p.value[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER, c(2,1)]), 
#          scale = "column", Colv = F, dendrogram = "none", Rowv = F,
#          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.ac.msigdb_Mm$G_PROTEIN_SIGNALING_COUPLED_TO_CAMP_NUCLEOTIDE_SECOND_MESSENGER],
#          col = mycol, trace = "none", density.info = "none",
#          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.cAMP.msigdb_Mm$GO_RESPONSE_TO_CAMP, c(7:9,1:3,4:6)], 
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.cAMP.msigdb_Mm$GO_RESPONSE_TO_CAMP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.cAMP.msigdb_Mm$GO_RESPONSE_TO_CAMP, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", #Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.cAMP.msigdb_Mm$GO_RESPONSE_TO_CAMP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)
###
heatmap.2(as.matrix(-log10(df_test[df_test$genes %in% gs.cAMP.msigdb_Mm$GO_RESPONSE_TO_CAMP, c(6,12)])), 
          scale = "column", Colv = F, dendrogram = "none", #Rowv = F,
          labRow = df_test$genes[df_test$genes %in% gs.cAMP.msigdb_Mm$GO_RESPONSE_TO_CAMP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.cAMP.msigdb_Mm$CAMP_UP.V1_UP, c(7:9,1:3,4:6)], 
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% gs.cAMP.msigdb_Mm$CAMP_UP.V1_UP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))
###
heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% gs.cAMP.msigdb_Mm$CAMP_UP.V1_UP, c(2,1)], 
          scale = "column", Colv = F, dendrogram = "none", #Rowv = F,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% gs.cAMP.msigdb_Mm$CAMP_UP.V1_UP],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10), cexCol = 0.7, srtCol = 0)
###

## Barcode plot
barcodeplot(efit3.contr$t[,2], index = idx.all.cAMP.msigdb$GO_CELLULAR_RESPONSE_TO_CAMP, main = "PS1KO vs WT: GO_CELLULAR_RESPONSE_TO_CAMP" )
dev.print(pdf, file="GSEA/Rplot_barcode3_cellRespToCAMP_ps1VSwt.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,1], index = idx.all.cAMP.msigdb$GO_CELLULAR_RESPONSE_TO_CAMP, main = "PS1KO+ISO vs PS1KO: GO_CELLULAR_RESPONSE_TO_CAMP" )
dev.print(pdf, file="GSEA/Rplot_barcode3_cellRespToCAMP_ps1isoVSps1.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,2], index = idx.all.cAMP.msigdb$GO_RESPONSE_TO_CAMP, main = "PS1KO vs WT: GO_RESPONSE_TO_CAMP")
dev.print(pdf, file="GSEA/Rplot_barcode3_respToCAMP_ps1VSwt.pdf", height=6, width=8)

barcodeplot(efit3.contr$t[,1], index = idx.all.cAMP.msigdb$GO_RESPONSE_TO_CAMP, main = "PS1KO+ISO vs PS1KO: GO_RESPONSE_TO_CAMP")
dev.print(pdf, file="GSEA/Rplot_barcode3_respToCAMP_ps1isoVSps1.pdf", height=6, width=8)


#13. Lysosomal Hydrolases (Ju-Hyun)
gs.clcn.JH_Mm<-  c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                      "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")
length(gs.clcn.JH_Mm) # 28

## convert into ids
idx.all.clcn.JH <- ids2indices(gs.clcn.JH_Mm, id = rownames(v3))
str(idx.all.clcn.JH) # 25

## roast
roast3.clcn.JH_ps1VSwt <- roast(v3, idx.all.lysoHydro, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.clcn.JH_ps1VSwt
write.table(roast3.clcn.JH_ps1VSwt, file = "GSEA/roast3.clcn.JH_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

roast3.clcn.JH_ps1isoVSps1 <- roast(v3, idx.all.clcn.JH, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.clcn.JH_ps1isoVSps1
write.table(roast3.clcn.JH_ps1isoVSps1, file = "GSEA/roast3.clcn.JH_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

par(mfrow=c(1,2))
barcodeplot(efit3.contr$t[,2], index = idx.all.clcn.JH$Set1, main="PS1KO vs WT: Chloride Channel Genes")
barcodeplot(efit3.contr$t[,1], index = idx.all.clcn.JH$Set1, main="PS1KO+ISO vs PS1KO: Chloride Channel Genes")
dev.print(pdf, "GSEA/Rplot_barcode3_clcn.JH_ps1VSwt_ps1isoVSwt.pdf", height=6, width=12)
#dev.print(tiff, "GSEA/Rplot_barcode3_clcn_ps1VSwt_ps1isoVSwt.tiff", units="px", res=800, height=4800, width=10400) # 200 MB

#14. THERE IS A COMBINATION OF ALL ABOVE TESTED GENE SETS CALLED "genesets_clcnRelated_070819.xlsx" put together by Sandip on 07/08/19. Let's run that.
geneSets_allClcn.Hs <- loadWorkbook("/Users/sandip_home/NKI/Gene Sets/geneSets_clcnRelated_070819.xlsx")
geneSets_allClcn.Hs <- readWorksheet(geneSets_allClcn.Hs, sheet = 1, header = T)
dim(geneSets_allClcn.Hs) # 201  42
head(geneSets_allClcn.Hs)

geneSets_allClcn.Hs <- geneSets_allClcn.Hs[-1, ]
head(geneSets_allClcn.Hs)
tail(geneSets_allClcn.Hs)
dim(geneSets_allClcn.Hs) # 200  42
class(geneSets_allClcn.Hs) # "data.frame"
lapply(geneSets_allClcn.Hs, class) # All "character"

## Now convert them to prop cases for Mouse. 
geneSets_allClcn.Mm <- as.data.frame(lapply(geneSets_allClcn.Hs, function(v) 
  ifelse(!is.na(v), paste0(toupper(substr(v,1,1)), tolower(substr(v,2,20))), v)
  ))
head(geneSets_allClcn.Mm)
tail(geneSets_allClcn.Mm)
dim(geneSets_allClcn.Mm) # 200  42

write.table(geneSets_allClcn.Mm, file = "~/NKI/Gene Sets/geneSets_clcnRelated_Mm_070819.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F, na = "")

# I manually changed some gene names after saving this file. So let's recall the saved file
geneSets_allClcn.Mm <- loadWorkbook("~/NKI/Gene Sets/geneSets_clcnRelated_Mm_070819_frmtd.xlsx")
geneSets_allClcn.Mm <- readWorksheet(geneSets_allClcn.Mm, sheet = 1, header = T)
head(geneSets_allClcn.Mm)
tail(geneSets_allClcn.Mm)
dim(geneSets_allClcn.Mm) # 200  42

# Convert into ids
idx.allClcn.SJD <- ids2indices(gene.sets = geneSets_allClcn.Mm, id = row.names(v3))
str(idx.allClcn.SJD)

# mroast
mroast3.allClcn_ps1VSwt <- mroast(v3, idx.allClcn.SJD, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.allClcn_ps1VSwt
write.table(mroast3.allClcn_ps1VSwt, file = "GSEA/mroast3.allClcn.SJD_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

mroast3.allClcn_ps1isoVSps1 <- mroast(v3, idx.allClcn.SJD, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.allClcn_ps1isoVSps1
write.table(mroast3.allClcn_ps1isoVSps1, file = "GSEA/mroast3.allClcn.SJD_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


#15. forskolin regulated transcription in neurons. (Benito E, Valor LM, Jimenez-Minchan M, Huber W and Barco A. Comparative transcriptomics identifies CREB 
# as a primary hub of activity-driven neuronal gene expression. J Neurosci. 31(50): 18237-50)
gs.fsk.Valor_Mm <- read.table(file = "GSEA/forskolin_regulated_transcription_NADTranscriptomics.txt", sep = "\t",
                              header = T, fill = T, row.names = NULL, quote = "", stringsAsFactors = F)
head(gs.fsk.Valor_Mm)
tail(gs.fsk.Valor_Mm)
dim(gs.fsk.Valor_Mm) # 518   6

length(unique(gs.fsk.Valor_Mm$Gene.symbol)) # 501
sum(is.na(gs.fsk.Valor_Mm$Gene.symbol)) # 13
sum(duplicated(gs.fsk.Valor_Mm$Gene.symbol)) # 17

gs.fsk.Valor_Mm <- gs.fsk.Valor_Mm[!duplicated(gs.fsk.Valor_Mm$Gene.symbol), ]
dim(gs.fsk.Valor_Mm) # 501   6

sum(gs.fsk.Valor_Mm$Gene.symbol=="") # NA
sum(is.na(gs.fsk.Valor_Mm$Gene.symbol)) # 1

gs.fsk.Valor_Mm <- gs.fsk.Valor_Mm[!(is.na(gs.fsk.Valor_Mm$Gene.symbol)), ]
dim(gs.fsk.Valor_Mm) # 500   6


# For fsk activated genes
gs.fskUP.valor_Mm <- gs.fsk.Valor_Mm[gs.fsk.Valor_Mm$Fold.change > 0, 2]
length(gs.fskUP.valor_Mm) # 316
length(unique(gs.fskUP.valor_Mm)) # 316
#sum(duplicated(gs.fskUP.valor_Mm)) # 10

#gs.fskUP.valor_Mm <- gs.fskUP.valor_Mm[!duplicated(gs.fskUP.valor_Mm)]
#length(gs.fskUP.valor_Mm) # 317

# For fsk repressed genes
gs.fskDN.valor_Mm <- gs.fsk.Valor_Mm[gs.fsk.Valor_Mm$Fold.change < 0, 2]
length(gs.fskDN.valor_Mm) # 184
length(unique(gs.fskDN.valor_Mm)) # 184
#sum(duplicated(gs.fskDN.valor_Mm)) # 142 (really?)

# Let's combine both genesets (need to make them of same lengths)
gs.fskDN.valor_Mm <- c(gs.fskDN.valor_Mm, rep(NA, 316-184))
length(gs.fskDN.valor_Mm) # 316

gs.fsk.Valor_Mm <- data.frame("fsk.up"=gs.fskUP.valor_Mm, "fsk.down"=gs.fskDN.valor_Mm, stringsAsFactors = F)
head(gs.fsk.Valor_Mm)
tail(gs.fsk.Valor_Mm)
dim(gs.fsk.Valor_Mm) # 316   2

write.table(gs.fsk.Valor_Mm, file = "GSEA/gs.fsk.valor_Mm.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


## convert into ids
idx.fsk.valor <- ids2indices(gs.fsk.Valor_Mm, id = rownames(v3))
str(idx.fsk.valor ) # 1:256; 1:119

## roast
mroast3.fsk.valor_ps1VSwt <- mroast(v3, idx.fsk.valor, design, contrast = contr.mat[,2], nrot = 1000000) # nrot=10^7 took long time
mroast3.fsk.valor_ps1VSwt
#write.table(mroast3.ac.msigdb_ps1VSwt, file = "GSEA/mroast3.ac.msigdb_ps1VSwt.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)

mroast3.fsk.valor_ps1isoVSps1 <- mroast(v3, idx.fsk.valor, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.fsk.valor_ps1isoVSps1
#write.table(mroast3.ac.msigdb_ps1isoVSps1, file = "GSEA/mroast3.ac.msigdb_ps1isoVSps1.txt", sep = "\t",
#            row.names = T, col.names = NA, quote = F)


#16. Er/Golgi related GO categories (On Ju-Hyun's request) 09/13/19
geneSets_ERtoGolgi.JH_Mm <- loadWorkbook("GSEA/ER.golgi_related_pathways/ER.golgi_related_pathways_GO.xlsx")
geneSets_ERtoGolgi.JH_Mm <- readWorksheet(geneSets_ERtoGolgi.JH_Mm, sheet = 1, header = T)
dim(geneSets_ERtoGolgi.JH_Mm) # 460   8
head(geneSets_ERtoGolgi.JH_Mm)

#geneSets_ERtoGolgi.JH_Mm <- geneSets_ERtoGolgi.JH_Mm[lapply(geneSets_ERtoGolgi.JH_Mm, function(x) 
#  ifelse(!is.na(x), !duplicated(x), x)), ]
sapply(geneSets_ERtoGolgi.JH_Mm, function(x) length(unique(x)))
rapply(geneSets_ERtoGolgi.JH_Mm, function(x) length(unique(x)))
# rapply(geneSets_ERtoGolgi.JH_Mm, function(x) length(!is.na(unique(x)))) # Same result as above
head(sapply(geneSets_ERtoGolgi.JH_Mm, function(x) unique(x)))
head(rapply(geneSets_ERtoGolgi.JH_Mm, function(x) unique(x)))
class(sapply(geneSets_ERtoGolgi.JH_Mm, function(x) unique(x))) # list
class(rapply(geneSets_ERtoGolgi.JH_Mm, function(x) unique(x))) # "character"

df <- sapply(geneSets_ERtoGolgi.JH_Mm, function(x) unique(x))
head(df) # list of 8

df2 <- sapply(df, '[', 
              seq(max(sapply(df, length))))
head(df2)
tail(df2)
class(df2) # matrix

df2 <- as.data.frame(df2)
head(df2)
tail(df2)
dim(df2) # 262   8
lapply(df2, class) # factor

df2 <- data.frame(lapply(df2, as.character), stringsAsFactors = F)
lapply(df2, class) # All character now

geneSets_ERtoGolgi.JH_Mm <- df2

## convert into ids
idx.all.ErToGolgi.JH <- ids2indices(gene.sets = geneSets_ERtoGolgi.JH_Mm, id = row.names(v3))
str(idx.all.ErToGolgi.JH)

## mroast
mroast3.ErToGolgi.JH_ps1VSwt <- mroast(v3, idx.all.ErToGolgi.JH, design, contrast = contr.mat[,2], nrot = 100000)
mroast3.ErToGolgi.JH_ps1VSwt

mroast3.ErToGolgi.JH_ps1isoVSps1 <- mroast(v3, idx.all.ErToGolgi.JH, design, contrast = contr.mat[,1], nrot = 100000)
mroast3.ErToGolgi.JH_ps1isoVSps1

## Heatmaps
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets_ERtoGolgi.JH_Mm$Golgi.vesicle.transport_GO0048193_262, 1:9], scale = "row",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% geneSets_ERtoGolgi.JH_Mm$Golgi.vesicle.transport_GO0048193_262], labCol = samplenames[1:9],
          dendrogram = "column",
          col = mycol, trace = "none", density.info =  "none", cexCol = 0.75, srtCol = 0, adjCol = c(0.5, NA),
          margins = c(4,6), lhei = c(2,10)) #  

heatmap.2(efit3.contr$t[row.names(efit3.contr) %in% geneSets_ERtoGolgi.JH_Mm$Golgi.vesicle.transport_GO0048193_262, c(2,1)], scale = "column",
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% geneSets_ERtoGolgi.JH_Mm$Golgi.vesicle.transport_GO0048193_262],# labCol = samplenames[1:9],
          col = mycol, trace = "none", density.info =  "none", cexCol = 0.75, srtCol = 0, adjCol = c(0.5, NA),
          margins = c(4,6), lhei = c(2,10))

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets_ERtoGolgi.JH_Mm$Regulation.of.ER.to.Golgi.vesicle.mediated.transport_GO0060628_17, 1:9], scale = "row",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% geneSets_ERtoGolgi.JH_Mm$Regulation.of.ER.to.Golgi.vesicle.mediated.transport_GO0060628_17], labCol = samplenames[1:9],
          dendrogram = "column",
          col = mycol, trace = "none", density.info =  "none", cexCol = 0.75, srtCol = 0, adjCol = c(0.5, NA),
          margins = c(4,6), lhei = c(2,10))

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets_ERtoGolgi.JH_Mm$Post.Golgi.vesicle.mediated.transport_GO0006892_122, 1:9], scale = "row",
          labRow = row.names(lcpm.norm3_5_c)[row.names(lcpm.norm3_5_c) %in% geneSets_ERtoGolgi.JH_Mm$Post.Golgi.vesicle.mediated.transport_GO0006892_122], labCol = samplenames[1:9],
          dendrogram = "column",
          col = mycol, trace = "none", density.info =  "none", cexCol = 0.75, srtCol = 0, adjCol = c(0.5, NA),
          margins = c(4,6), lhei = c(2,10))


#########################################################################################################################################################

### Using msigdbr package
library(msigdbr)

# Retrieve mouse C2 (curated) gene sets
Mm.c2_msigdb <- msigdbr(species = "Mus musculus", category = "C2")
class(Mm.c2_msigdb) # "tbl_df"     "tbl"        "data.frame"
head(Mm.c2_msigdb)

Mm.c2_msigdb_list = Mm.c2_msigdb %>% split(x = .$gene_symbol, f = .$gs_name)
Mm.c2_msigdb_list[1:3]


idx.all.msigdb <- ids2indices(Mm.c2_msigdb_list, id = rownames(v3))
cam3.msigdb_ps1VSwt <- camera(v3, idx.all.msigdb, design, contrast = contr.mat[,2])
head(cam3.msigdb_ps1VSwt,15)
write.table(cam3.msigdb_ps1VSwt, file = "GSEA/cam3.msigdb_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

cam3.msigdb_ps1isoVSwt <- camera(v3, idx.all.msigdb, design, contrast = contr.mat[,4])
head(cam3.msigdb_ps1isoVSwt,15)
write.table(cam3.msigdb_ps1isoVSwt, file = "GSEA/cam3.msigdb_ps1isoVSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

cam3.msigdb_ps1isoVSps1 <- camera(v3, idx.all.msigdb, design, contrast = contr.mat[,1])
head(cam3.msigdb_ps1isoVSps1,15)
write.table(cam3.msigdb_ps1isoVSps1, file = "GSEA/cam3.msigdb_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


library(fgsea)
fgsea.msigdb_ps1VSwt <- fgsea(pathways = Mm.c2_msigdb_list, stats = efit3.contr$t[ ,2], minSize = 10, maxSize = 1000, nperm = 100000)
fgsea.msigdb_ps1VSwt[order(pval), ][1:50, ]

fgsea.msigdb_ps1isoVSwt <- fgsea(pathways = Mm.c2_msigdb_list, stats = efit3.contr$t[ ,4], minSize = 10, maxSize = 1000, nperm = 100000)
fgsea.msigdb_ps1isoVSwt[order(pval), ][1:50, ]

fgsea.msigdb_ps1isoVSps1 <- fgsea(pathways = Mm.c2_msigdb_list, stats = efit3.contr$t[ ,1], minSize = 10, maxSize = 1000, nperm = 100000)
fgsea.msigdb_ps1isoVSps1[order(pval), ][1:50, ]

library(reactome.db)
reactome3.1.msigdb_ps1VSwt <- reactomePathways(names(efit3.1.contr$t[,2]))
reactome3.1.msigdb_ps1VSwt <- fgsea(reactome3.1.msigdb_ps1VSwt, efit3.1.contr$t[,2], nperm=10000, maxSize=1000)
head(reactome3.1.msigdb_ps1VSwt[order(pval), ], 20)

reactome3.1.msigdb_ps1isoVSwt <- reactomePathways(names(efit3.1.contr$t[,4]))
reactome3.1.msigdb_ps1isoVSwt <- fgsea(reactome3.1.msigdb_ps1isoVSwt, efit3.1.contr$t[,4], nperm=10000, maxSize=1000)
head(reactome3.1.msigdb_ps1isoVSwt[order(pval), ], 30)

reactome3.1.msigdb_ps1isoVSps1 <- reactomePathways(names(efit3.1.contr$t[,1]))
reactome3.1.msigdb_ps1isoVSps1 <- fgsea(reactome3.1.msigdb_ps1isoVSps1, efit3.1.contr$t[,1], nperm=10000, maxSize=1000)
head(reactome3.1.msigdb_ps1isoVSps1[order(pval), ], 20)


#Mm.c2_msigdb <- as.data.frame(Mm.c2_msigdb)
#class(Mm.c2_msigdb) # "data.frame"
#head(Mm.c2_msigdb)

## Checking curated induction genes; How individual genes change after iso txt
res3.ps1isoVSwt_ind <- res3.ps1isoVSwt[res3.ps1isoVSwt$genes %in% geneSets110_Mm$Curated.Induction.23., ]
dim(res3.ps1isoVSwt_ind) # 23  7
res3.ps1isoVSwt_ind

res3.ps1VSwt_ind <- res3.ps1VSwt[res3.ps1VSwt$genes %in% geneSets110_Mm$Curated.Induction.23., ]
res3.ps1VSwt_ind <- res3.ps1VSwt_ind[match(res3.ps1isoVSwt_ind$genes, res3.ps1VSwt_ind$genes), ]
dim(res3.ps1VSwt_ind) # 23  7
res3.ps1VSwt_ind



######################################################################   goseq   ######################################################################
library(goseq)
# Format the DE genes into a vector suitable for use with goseq
genes3.ps1VSwt <- as.integer(p.adjust(res3.ps1VSwt$P.Value[res3.ps1VSwt$logFC!=0], method="BH") < 0.05)
names(genes3.ps1VSwt) <- row.names(res3.ps1VSwt[res3.ps1VSwt$logFC != 0, ])
table(genes3.ps1VSwt)
#     0     1
#  3961 11697
###
genes3.ps1VSwt_sigUP <- as.integer(p.adjust(res3.ps1VSwt$P.Value[res3.ps1VSwt$logFC > 0], method = "BH") < 0.05)
names(genes3.ps1VSwt_sigUP) <- row.names(res3.ps1VSwt[res3.ps1VSwt$logFC > 0, ])
table(genes3.ps1VSwt_sigUP)
#     0     1
#   2003 5939
###
genes3.ps1VSwt_sigDN <- as.integer(p.adjust(res3.ps1VSwt$P.Value[res3.ps1VSwt$logFC < 0], method = "BH") < 0.05)
names(genes3.ps1VSwt_sigDN) <- row.names(res3.ps1VSwt[res3.ps1VSwt$logFC < 0, ])
table(genes3.ps1VSwt_sigDN)
#     0     1
#  1958  5758


# Determining Genome and Gene ID
head(supportedOrganisms())
supportedOrganisms()[supportedOrganisms()$Genome == "mm9", ] ## mm10 is not supported yet.

## GO analysis
# Fitting the Probability Weighting Function (PWF)----We will let goseq automatically fetch this data from its databases.
pwf3.ps1VSwt <- nullp(genes3.ps1VSwt, "mm9", "geneSymbol")
###
pwf3.ps1VSwt_sigUP <- nullp(genes3.ps1VSwt_sigUP, "mm9", "geneSymbol")
###
pwf3.ps1VSwt_sigDN <- nullp(genes3.ps1VSwt_sigDN, "mm9", "geneSymbol")


# Using the Wallenius approximation (recommended contrary to random sampling/approximation)
GO3.wall.ps1VSwt <- goseq(pwf = pwf3.ps1VSwt, "mm9", "geneSymbol")
head(GO3.wall.ps1VSwt, 30)

GO3.BP.wall.ps1VSwt <- goseq(pwf = pwf3.ps1VSwt, "mm9", "geneSymbol", test.cats = c("GO:BP"))
head(GO3.BP.wall.ps1VSwt, 30)
GO3.BP.wall.ps1isoVSwt[151:200,]

Kegg3.wall.ps1VSwt <- goseq(pwf = pwf3.ps1VSwt, "mm9", "geneSymbol", test.cats = "KEGG")
head(Kegg3.wall.ps1VSwt, 30)
###
GO3.wall.ps1VSwt_sigUP <- goseq(pwf = pwf3.ps1VSwt_sigUP, "mm9", "geneSymbol")
head(GO3.wall.ps1VSwt_sigUP, 30)

GO3.BP.wall.ps1VSwt_sigUP <- goseq(pwf = pwf3.ps1VSwt_sigUP, "mm9", "geneSymbol", test.cats = c("GO:BP"))
GO3.BP.wall.ps1VSwt_sigUP[30:60,]
#GO3.BP.wall.ps1isoVSwt[151:200,]

Kegg3.wall.ps1VSwt_sigUP <- goseq(pwf = pwf3.ps1VSwt_sigUP, "mm9", "geneSymbol", test.cats = "KEGG")
head(Kegg3.wall.ps1VSwt_sigUP, 30)
###
GO3.wall.ps1VSwt_sigDN <- goseq(pwf = pwf3.ps1VSwt_sigDN, "mm9", "geneSymbol")
head(GO3.wall.ps1VSwt_sigDN, 30)

GO3.BP.wall.ps1VSwt_sigDN <- goseq(pwf = pwf3.ps1VSwt_sigDN, "mm9", "geneSymbol", test.cats = c("GO:BP"))
GO3.BP.wall.ps1VSwt_sigDN[30:60,]
GO3.BP.wall.ps1VSwt_sigDN[61:100,]
GO3.BP.wall.ps1VSwt_sigDN[101:150,]
GO3.BP.wall.ps1VSwt_sigDN[151:200,]
GO3.BP.wall.ps1VSwt_sigDN[201:300,]
GO3.BP.wall.ps1VSwt_sigDN[301:400,]

Kegg3.wall.ps1VSwt_sigDN <- goseq(pwf = pwf3.ps1VSwt_sigDN, "mm9", "geneSymbol", test.cats = "KEGG")
head(Kegg3.wall.ps1VSwt_sigDN, 30)

#install("KEGGREST")
library(KEGGREST)
head(keggList("mmu"))

keggGet("mmu00562")[[1]]$NAME

keggids_ps1VSwt_sigUP <- lapply(Kegg3.wall.ps1VSwt_sigUP$category, function(x) 
  paste0("mmu", x))
head(keggids_ps1VSwt_sigUP)

#keggGet(keggids[1:10])[[rep(1,10)]]$NAME
sapply(keggGet(keggids_ps1VSwt_sigUP[11:20]), "[[", "NAME")  ## Ref: https://support.bioconductor.org/p/66777/#66817
###
keggids_ps1VSwt_sigDN <- lapply(Kegg3.wall.ps1VSwt_sigDN$category, function(x) 
  paste0("mmu", x))
head(keggids_ps1VSwt_sigDN)

#keggGet(keggids[1:10])[[rep(1,10)]]$NAME
sapply(keggGet(keggids_ps1VSwt_sigDN), "[[", "NAME")
  


library(GO.db)
for (go in Kegg3.wall.ps1VSwt$category[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

test <- as.list(org.Mm.egPATH)
head(test)

genes3.ps1isoVSwt <- as.integer(p.adjust(res3.ps1isoVSwt$P.Value[res3.ps1isoVSwt$logFC!=0], method="BH") < 0.05)
names(genes3.ps1isoVSwt) <- row.names(res3.ps1isoVSwt[res3.ps1isoVSwt$logFC != 0, ])
table(genes3.ps1isoVSwt)
#     0     1
#  3919 11739

## GO analysis
# Fitting the Probability Weighting Function (PWF)----We will let goseq automatically fetch this data from its databases.
pwf3.ps1isoVSwt <- nullp(genes3.ps1isoVSwt, "mm9", "geneSymbol")

# Using the Wallenius approximation (recommended contrary to random sampling/approximation)
GO3.wall.ps1isoVSwt <- goseq(pwf = pwf3.ps1isoVSwt, "mm9", "geneSymbol")
head(GO3.wall.ps1isoVSwt, 30)

GO3.BP.wall.ps1isoVSwt <- goseq(pwf = pwf3.ps1isoVSwt, "mm9", "geneSymbol", test.cats = c("GO:BP"))
head(GO3.BP.wall.ps1isoVSwt, 30)
GO3.BP.wall.ps1isoVSwt[31:60,]


genes3.ps1isoVSps1 <- as.integer(p.adjust(res3.ps1isoVSps1$P.Value[res3.ps1isoVSps1$logFC!=0], method="BH") < 0.05)
names(genes3.ps1isoVSps1) <- row.names(res3.ps1isoVSps1[res3.ps1isoVSps1$logFC != 0, ])
table(genes3.ps1isoVSps1)
#     0     1
# 13533  2125

## GO analysis
# Fitting the Probability Weighting Function (PWF)----We will let goseq automatically fetch this data from its databases.
pwf3.ps1isoVSps1 <- nullp(genes3.ps1isoVSps1, "mm9", "geneSymbol")

# Using the Wallenius approximation (recommended contrary to random sampling/approximation)
GO3.wall.ps1isoVSps1 <- goseq(pwf = pwf3.ps1isoVSps1, "mm9", "geneSymbol")
head(GO3.wall.ps1isoVSps1, 30)

GO3.BP.wall.ps1isoVSps1 <- goseq(pwf = pwf3.ps1isoVSps1, "mm9", "geneSymbol", test.cats = c("GO:BP"))
head(GO3.BP.wall.ps1isoVSps1, 30)
GO3.BP.wall.ps1isoVSps1[31:60,]


##################################################################  goseq complete #####################################################################

#################################### EdgeR Quasi Likelyhood (QL) ############################################
## Dispersion estimation
y3.estDisp <- estimateDisp(y.norm3, design = design, robust = T) # This returns a DGEList
plotBCV(y3.estDisp)

fit3.QL <- glmQLFit(y3.estDisp, design, robust = T)
plotQLDisp(fit3.QL)
summary(fit3.QL$df.prior)

## Differential expression analysis
res3_QLF.ps1VSwt <- glmQLFTest(fit3.QL, contrast = contr.mat[,2])
res3_QL.ps1VSwt <- topTags(res3_QLF.ps1VSwt, n=Inf)
res3_QL.tbl.ps1VSwt <- res3_QL.ps1VSwt$table
head(res3_QL.tbl.ps1VSwt)

res3_QLF.ps1isoVSwt <- glmQLFTest(fit3.QL, contrast = contr.mat[,4])
res3_QL.ps1isoVSwt <- topTags(res3_QLF.ps1isoVSwt, n=Inf)
res3_QL.tbl.ps1isoVSwt <- res3_QL.ps1isoVSwt$table
head(res3_QL.tbl.ps1isoVSwt)

res3_QLF.ps1isoVSps1 <- glmQLFTest(fit3.QL, contrast = contr.mat[,1])
res3_QL.ps1isoVSps1 <- topTags(res3_QLF.ps1isoVSps1, n=Inf)
res3_QL.tbl.ps1isoVSps1 <- res3_QL.ps1isoVSps1$table
head(res3_QL.tbl.ps1isoVSps1)


dt3_QL.ps1VSwt <- decideTestsDGE(res3_QLF.ps1VSwt)
summary(dt3_QL.ps1VSwt)  ## Almost exactly SAME numbers as in limma-voom
#        1*ps1 -1*wt
# Down          5752
# NotSig        3955
# Up            5951

dt3_QL.ps1isoVSwt <- decideTestsDGE(res3_QLF.ps1isoVSwt)
summary(dt3_QL.ps1isoVSwt)  ## Almost exactly SAME numbers as in limma-voom
#        1*ps1iso -1*wt
# Down             5762
# NotSig           3920
# Up               5976

dt3_QL.ps1isoVSps1 <- decideTestsDGE(res3_QLF.ps1isoVSps1)
summary(dt3_QL.ps1isoVSps1)  ## Almost exactly SAME numbers as in limma-voom
#        -1*ps1 1*ps1iso
# Down              1100
# NotSig           13483
# Up                1075

res3_QL.tbl.ps1isoVSps1[res3_QL.tbl.ps1isoVSps1$genes %in% c("Mitf","Ostm1","Clcn5","Clcn7","Clcn2","Clcn3"), ]

######################################################      DONE      ######################################################

###################################################### FOR MANUSCRIPT ######################################################
## Pie chart (In excel)

## Venn-diagram
par(mfrow=c(1,1))
vennDiagram(dt3[,c(1,2,5,6)], circle.col = c("blue", "red", "green", "orange"))
dev.print(pdf, "Rplot_Venn3_4way.pdf", height=6, width=8)

vennDiagram(dt3[,c(1,6)], circle.col = c("green","orange"), names = c("PS1KO+ISO vs PS1KO", "WT+ISO vs WT"), cex = c(1,1,1), mar = c(0.05,4,0.05,4))
dev.print(pdf, "Rplot_venn3_ps1isoVSps1-wtisoVSwt.pdf", height=6, width=8)

vennDiagram(dt3[,c(2,1)], circle.col = c("turquoise", "salmon"), include = c("up", "down"), cex = c(1,1,1),
            names = c("PS1KO vs WT", "PS1KO+ISO vs PS1KO"), mar = c(0.05,4,0.05,4))
dev.print(pdf, "Rplot_Venn3_ps1VSwt-ps1isoVSps1.pdf", height=6, width=8)

vennDiagram(dt3[,c(2,1)], circle.col = c("blue", "red", "green", "orange"))
vennDiagram(dt3[,c(2,1,6)], circle.col = c("blue", "red", "green", "orange"), include = c("up", "down"), 
            names = c("PS1KO vs WT", "PS1KO+ISO vs PS1KO", "WT+ISO vs WT"), cex = c(1,0.5,1), mar = c(0.05,4,0.05,4))
dev.print(pdf, "Rplot_Venn3_3way.pdf", height=6, width=8)

summary(dt3_t)
vennDiagram(dt3_p0.01[,c(1,6)], circle.col = c("green","orange"), names = c("PS1KO+ISO vs PS1KO", "WT+ISO vs WT"), cex = c(1,1,1), mar = c(0.1,2,0.1,2))
#dev.print(pdf, "Rplot_venn3_ps1isoVSps1-wtisoVSwt.pdf", height=6, width=8)




library(VennDiagram)
#venn.diagram(list("PS1KO UP"=res3.ps1))

## Heatmap (with gplots)
library(ggplot2)

library(RColorBrewer)
par(mfrow=c(1,1))
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[row.names(v3$E) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% geneSets110_Mm$protein.transport.40.],
          col = mycol, trace = "none", density.info = "none",
          margins = c(8,6), lhei = c(2,10))

mycol2 <- colorpanel(100, "blue", "red")
heatmap.2(v3$E[row.names(v3$E) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", Rowv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% geneSets110_Mm$protein.transport.40.],
          col = mycol2, trace = "none", density.info = "none",
          margins = c(8,10))

heatmap.2(v3$E[row.names(v3$E) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          labRow = v3$genes$genes[v3$genes$genes %in% geneSets110_Mm$protein.transport.40.],
          col = mycol2, trace = "none", density.info = "none",
          margins = c(8,10))

mycol3 <- colorpanel(1000, "green", "black", "red")
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          Colv = "dendrogram", dendrogram = "none", scale = "row", Rowv = "dendrogram", 
          col = mycol, trace = "none", density.info = "none")

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          Colv = "dendrogram", dendrogram = "none", Rowv = T, #scale = "row"
          col = mycol, trace = "none", density.info = "none") # Bad

mycol3 <- colorpanel(1000, "green", "black", "red")
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          Colv = FALSE, dendrogram = "none", # scale = "row"
          col = mycol, trace = "none", density.info = "none")

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          scale = "row", Colv = TRUE, dendrogram = "column",
          col = mycol, trace = "none", density.info = "none")

mycol4 <- colorpanel(1000, "white", "red")
heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          col = mycol4, trace = "none", density.info = "none")

heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(1,2)],
          scale = "row", Colv = TRUE, dendrogram = "none",
          col = mycol, trace = "none", density.info = "none")

heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0,
          sepwidth = c(0.1,0.01), sepcolor = "black", colsep = 1:ncol(efit3.contr$t[row.names(efit3.contr$t),]), rowsep = 1:nrow(efit3.contr$t[row.names(efit3.contr$t),]))
            ## USE THIS

#### OR ####
nrow(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(1,2)]) # 31
ncol(nrow(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(1,2)])) # 2
heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0,
          sepwidth = c(0.0001,0.005), sepcolor = "black", colsep = 0:2, rowsep = 0:31)  ## SAME AS ABOVE WITH BORDERS AROUND CELLS


heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = bluered(256), trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0,
          sepwidth = c(0.0001,0.005), sepcolor = "black", colsep = 0:2, rowsep = 0:31) # SAME AS ABOVE

heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% geneSets110_Mm$protein.transport.40.],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,5), cexCol = 1, srtCol = 0, #margins = c(8,6)
          sepwidth = c(0.0001,0.005), sepcolor = "black", colsep = 0:2, rowsep = 0:31)
############## Example of layout ##################
x <- pmin(3, pmax(-3, stats::rnorm(50)))
y <- pmin(3, pmax(-3, stats::rnorm(50)))
xhist <- hist(x, breaks = seq(-3,3,0.5), plot = FALSE)
yhist <- hist(y, breaks = seq(-3,3,0.5), plot = FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(-3, 3)
yrange <- c(-3, 3)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(3,1), c(1,3), TRUE)
layout.show(nf)

par(mar = c(3,3,1,1))
plot(x, y, xlim = xrange, ylim = yrange, xlab = "", ylab = "")
par(mar = c(0,3,1,1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)
par(mar = c(3,0,1,1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE)
#######################################################

head(df_test)
df_test2 <- df_test
df_test2$lfc_FDR.x <- sign(df_test2$logFC.x) * -log10(df_test2$adj.P.Val.x)
head(df_test2)
head(df_test2)
df_test2$lfc_FDR.y <- sign(df_test2$logFC.y) * -log10(df_test2$adj.P.Val.y)
head(df_test2)
dim(df_test2) # 15658    15

heatmap.2(as.matrix(df_test2[df_test2$genes %in% geneSets110_Mm$protein.transport.40., c(14,15)]),
          Colv = T, dendrogram = "none", Rowv = T,
          labRow = df_test2$genes[df_test2$genes %in% geneSets110_Mm$protein.transport.40.],
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0)  ## Maybe USE THIS?

heatmap.2(efit3.contr$coefficients[row.names(efit3.contr$coefficients) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = T, dendrogram = "none", Rowv = T,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% geneSets110_Mm$protein.transport.40.],
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0) ## Or use this.

mycol5 = colorpanel(100, "blue", "white", "red")
mycol6 = colorpanel(10000, "blue", "white", "red")
heatmap.2(efit3.contr$coefficients[row.names(efit3.contr$coefficients) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = T, dendrogram = "none", Rowv = T,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% geneSets110_Mm$protein.transport.40.],
          col = mycol6, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0) # mycol is better than mycol5

heatmap.2(efit3.contr$coefficients[row.names(efit3.contr$coefficients) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = T, dendrogram = "none", Rowv = T,
          labRow = row.names(efit3.contr)[row.names(efit3.contr) %in% geneSets110_Mm$protein.transport.40.],
          col = bluered(256), trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0)

mycol5 <- colorpanel(100, "blue", "white", "red")
heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0)  ## SAME RESULT. USED in MS.

heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = "dendrogram", dendrogram = "none", 
          col = mycol, trace = "none", density.info = "none",
          lhei = c(0.5,2.5), lwid = c(1,2), margins = c(10,20),cexCol = 1, srtCol = 0, adjCol = c(0.5,1)) # margins = c(8,6)  ## USE THIS IN MS.
#dev.print()

heatmap.2(efit3.contr$lods[row.names(efit3.contr$lods) %in% geneSets110_Mm$protein.transport.40., c(2,1)],
          Colv = "dendrogram", dendrogram = "none", scale = "column",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0)

heatmap.2(v3$E[row.names(v3$E) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          scale = "row", dendrogram = "column", Colv = v3$weights[rownames(v3$weights) %in% geneSets110_Mm$protein.transport.40., ],
          col = mycol, trace = "none", density.info = "none")

heatmap.2(v3$E[row.names(v3$E) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          scale = "row", dendrogram = "none", Colv = "dendrogram",
          col = mycol, trace = "none", density.info = "none")

heatmap.2(v3$E[row.names(v3$E) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          dendrogram = "none", Colv = "dendrogram",
          col = mycol, trace = "none", density.info = "none")

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$protein.transport.40., c(7:9,1:3,4:6)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none")


heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$All.trasport.genes.232., c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1, srtCol = 0)  ## USE THIS

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$All.trasport.genes.232., c(7:9,1:3,4:6)],
          Colv = "dendrogram", dendrogram = "none", scale = "row",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6))


# For Curated induction gene set heatmap
library(gplots)
par(mfrow=c(1,1))
heatmap.2(lcpm.norm3_5[row.names(lcpm.norm3_5) %in% geneSets110_Mm$Curated.Induction.23., c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6))  ## Not bad

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% geneSets110_Mm$Curated.Induction.23., c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6)) ## Better

heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% geneSets110_Mm$Curated.Induction.23., c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6))


## For comparison of GO terms
df_go.fdr <- loadWorkbook("GO/GO_excel_formatted/GO_BP2018_enrichr_AllSigUP_ps1isoVSps1.xlsx")
df_go.fdr <- readWorksheet(df_go.fdr, sheet = 2, startRow = 35,endRow = 50, header = T)
head(df_go.fdr)

# Select only required columns
df_go.fdr <- df_go.fdr[, c(1,6,12)]
head(df_go.fdr)
row.names(df_go.fdr) <- df_go.fdr$Col1
df_go.fdr <- df_go.fdr[, -1]
head(df_go.fdr)
colnames(df_go.fdr) <- c("FDR:PS1isoVSps1", "FDR:PS1vsWT")
df_go.fdr
lapply(df_go.fdr, class) # Both numeric
class(df_go.fdr[, c(2:1)]) # "data.frame"
class(as.matrix(df_go.fdr[, c(2:1)])) # "matrix"

# Heatmap
heatmap.2(as.matrix(df_go.fdr[, c(2:1)]), scale = "none", dendrogram = "none", 
          Rowv = "dendrogram", Colv = "dendrogram",
          labRow = row.names(df_go.fdr),
          col = mycol4, trace = "none", density.info = "none", 
          lhei = c(1,3), lwid = c(1,1),#, margins = c(4,4)), 
          cexCol = 0.7, srtCol = 0, adjCol = c(NA, 0.5))
dev.print(pdf, "ms_figures/heatmap_GO_comparisons.pdf", height=6, width=8)

## Comparison of transcriptomic profiles (BD and HS Fibro)
# Import DEG df from HS Fibro
df_dge.ps1mutVSwt_Fibro.Hs <- read.table("~/Bioinformatics/PS1mut_PS2mut_Fibro_Hs/R_dir_ps1mut_ps2mut_fibro_Hs/res.ps1_mutVSwt.txt", sep = "\t", 
                                         header = T, stringsAsFactors = F, quote = "")
head(df_dge.ps1mutVSwt_Fibro.Hs)
df_dge.ps1mutVSwt_Fibro.Hs <- df_dge.ps1mutVSwt_Fibro.Hs[df_dge.ps1mutVSwt_Fibro.Hs$adj.P.Val < 0.05, ]
dim(df_dge.ps1mutVSwt_Fibro.Hs)
# 10660     7
head(res3.ps1VSwt)
sum(duplicated(res3.ps1VSwt$genes)) # 0

df_dge.ps1VSwt_BD <- res3.ps1VSwt
df_dge.ps1VSwt_BD$genes <- lapply(df_dge.ps1VSwt_BD$genes, function(v) 
  toupper(v))
head(df_dge.ps1VSwt_BD)

df_dge.ps1VSwt_BD$genes <- as.character(df_dge.ps1VSwt_BD$genes)

df_dge.ps1VSwt_BD <- df_dge.ps1VSwt_BD[df_dge.ps1VSwt_BD$adj.P.Val < 0.05, ]
dim(df_dge.ps1VSwt_BD)
# 11697     8

## Keeping same genes in both df.
df_dge.ps1VSwt_BD <- df_dge.ps1VSwt_BD[df_dge.ps1VSwt_BD$genes %in% df_dge.ps1mutVSwt_Fibro.Hs$genes, ]
dim(df_dge.ps1VSwt_BD)
# 6018    8

df_dge.ps1mutVSwt_Fibro.Hs <- df_dge.ps1mutVSwt_Fibro.Hs[df_dge.ps1mutVSwt_Fibro.Hs$genes %in% df_dge.ps1VSwt_BD$genes, ]
dim(df_dge.ps1mutVSwt_Fibro.Hs)
# 6017    8

## Putting them in equal order of genes
sum(duplicated(df_dge.ps1mutVSwt_Fibro.Hs$genes)) # 0
sum(duplicated(df_dge.ps1VSwt_BD$genes)) # 1

df_dge.ps1VSwt_BD <- df_dge.ps1VSwt_BD[!duplicated(df_dge.ps1VSwt_BD$genes), ]

df_dge.ps1VSwt_BD<- df_dge.ps1VSwt_BD[match(df_dge.ps1mutVSwt_Fibro.Hs$genes, df_dge.ps1VSwt_BD$genes), ]

head(df_dge.ps1mutVSwt_Fibro.Hs)
head(df_dge.ps1VSwt_BD)

tail()

cor.test(df_dge.ps1VSwt_BD$logFC, df_dge.ps1mutVSwt_Fibro.Hs$logFC)

plot(df_dge.ps1VSwt_BD$logFC, df_dge.ps1mutVSwt_Fibro.Hs$logFC)

## Combine both df

df_combd <- merge(df_dge.ps1mutVSwt_Fibro.Hs, df_dge.ps1VSwt_BD, by = "genes")
head(df_combd)

library(ggplot2)
ggplot(data = df_combd, aes(x=sign(logFC.x), y=sign(logFC.y)))+
         geom_point(color = "red") + 
         geom_smooth(method = "lm")

reg1 <- lm(logFC.x*-log10(adj.P.Val.x) ~ logFC.y*-log10(adj.P.Val.y), data = df_combd)
summary(reg1)

with(df_combd, plot(logFC.x*-log10(adj.P.Val.x), logFC.y*-log10(adj.P.Val.y)))
abline(reg1)
###
df2_dge.ps1mutVSwt_Fibro.Hs <- read.table("~/Bioinformatics/PS1mut_PS2mut_Fibro_Hs/R_dir_ps1mut_ps2mut_fibro_Hs/res.ps1_mutVSwt.txt", sep = "\t", 
                                         header = T, stringsAsFactors = F, quote = "")
head(df2_dge.ps1mutVSwt_Fibro.Hs)
dim(df2_dge.ps1mutVSwt_Fibro.Hs) # 19120     7

#f_dge.ps1mutVSwt_Fibro.Hs <- df_dge.ps1mutVSwt_Fibro.Hs[df_dge.ps1mutVSwt_Fibro.Hs$adj.P.Val < 0.05, ]
#dim(df_dge.ps1mutVSwt_Fibro.Hs)
# 10660     7
head(res3.ps1VSwt)
sum(duplicated(res3.ps1VSwt$genes)) # 0

df2_dge.ps1VSwt_BD <- res3.ps1VSwt
df2_dge.ps1VSwt_BD$genes <- lapply(df2_dge.ps1VSwt_BD$genes, function(v) 
  toupper(v))
head(df2_dge.ps1VSwt_BD)
dim(df2_dge.ps1VSwt_BD) # 15658     7

df2_dge.ps1VSwt_BD$genes <- as.character(df2_dge.ps1VSwt_BD$genes)

#df_dge.ps1VSwt_BD <- df_dge.ps1VSwt_BD[df_dge.ps1VSwt_BD$adj.P.Val < 0.05, ]
#dim(df_dge.ps1VSwt_BD)
# 11697     8

## Keeping same genes in both df.
df2_dge.ps1mutVSwt_Fibro.Hs <- df2_dge.ps1mutVSwt_Fibro.Hs[df2_dge.ps1mutVSwt_Fibro.Hs$genes %in% df2_dge.ps1VSwt_BD$genes, ]
dim(df2_dge.ps1mutVSwt_Fibro.Hs)
# 11098     7

df2_dge.ps1VSwt_BD <- df2_dge.ps1VSwt_BD[df2_dge.ps1VSwt_BD$genes %in% df2_dge.ps1mutVSwt_Fibro.Hs$genes, ]
dim(df2_dge.ps1VSwt_BD)
# 11099     7

## Putting them in equal order of genes
sum(duplicated(df2_dge.ps1mutVSwt_Fibro.Hs$genes)) # 0
sum(duplicated(df2_dge.ps1VSwt_BD$genes)) # 1

df2_dge.ps1VSwt_BD <- df2_dge.ps1VSwt_BD[!duplicated(df2_dge.ps1VSwt_BD$genes), ]

df2_dge.ps1VSwt_BD<- df2_dge.ps1VSwt_BD[match(df2_dge.ps1mutVSwt_Fibro.Hs$genes, df2_dge.ps1VSwt_BD$genes), ]

head(df2_dge.ps1mutVSwt_Fibro.Hs)
head(df2_dge.ps1VSwt_BD)


cor.test(df2_dge.ps1VSwt_BD$logFC[df2_dge.ps1VSwt_BD$logFC < 0][1:6000], df2_dge.ps1mutVSwt_Fibro.Hs$logFC[df2_dge.ps1mutVSwt_Fibro.Hs$logFC < 0][1:6000])

plot(df_dge.ps1VSwt_BD$logFC[df2_dge.ps1VSwt_BD$logFC < 0][1:6000], df_dge.ps1mutVSwt_Fibro.Hs$logFC[df2_dge.ps1mutVSwt_Fibro.Hs$logFC < 0][1:6000])
abline()

## Combine both df

df2_combd <- merge(df2_dge.ps1mutVSwt_Fibro.Hs, df2_dge.ps1VSwt_BD, by = "genes")
head(df2_combd)

df3_combd <- df2_combd[df2_combd$logFC.x < 0, ]
dim(df3_combd) # 5480   13
library(ggplot2)
library(dplyr)
ggplot(data = df3_combd, aes(x=logFC.x, y=logFC.y))+
  geom_point(color = "red") + 
  geom_smooth(method = "lm")

reg2 <- lm(logFC.x*-log10(adj.P.Val.x) ~ logFC.y*-log10(adj.P.Val.y), data = df_combd)
summary(reg2)

with(df_combd, plot(logFC.x*-log10(adj.P.Val.x), logFC.y*-log10(adj.P.Val.y)))
abline(reg1)
###
reg3 <- lm(logFC.x*-log10(adj.P.Val.x) ~ logFC.y*-log10(adj.P.Val.y), data = df3_combd)
summary(reg3)

with(df3_combd, plot(logFC.x*-log10(adj.P.Val.x), logFC.y*-log10(adj.P.Val.y)))
abline(reg3)


## For whole genome
library(gplots)
par(mfrow=c(1,1))
ps1isoVSps1.sigGenes3 <- res3.ps1isoVSps1$genes[res3.ps1isoVSps1$adj.P.Val < 0.05]
i <- which(v3$genes$genes %in% ps1isoVSps1.sigGenes3)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v3$E[i,], scale = "row",
          labRow = "", labCol = samplenames,
          col = mycol, trace = "none", density.info =  "none",
          margins = c(8,6), lhei = c(2,10), dendrogram = "column")
dev.print(pdf, "Rplot_heatmap_ps1isoVSps1.pdf", height=10, width=8)


labCol = colnames(v3$E)[1:9]

## Transport related genes from Ju-Hyun's template
gs.transport_Mm <- c("Cdkn1a","Cdkn2c","Bub1","Casp8ap2","Ccnd1","Ccne1","Cdc6","Cdk1","Cenpe","Cenpv","Chek1","Dbf4","Ercc6l",
                       "Gmnn","Haus3","Hells","Hmga2","Id4","Kif11","Kif20b","Krt7","Mad2l1","Mcm7","Mcm8","Mlf1","Myb","Nae1",
                       "Nedd9","Npm1","Nuf2","Nup37","Nup43","Phgdh","Pin1","Pola1","Rad51","Ran","Ranbp1","Rcc1","Rcc2","Rif1",
                        "Rprm","Seh1l","Ska1","Skp2","Smc4","Tipin","Tpd52l1")
class(gs.transport_Mm) # "character"

idx.gs.transport <- ids2indices(gs.transport_Mm, id = rownames(v3))
str(idx.gs.transport) # Rprm gene is not present in V3.

roast3.gs.transport_ps1VSwt <- roast(v3, idx.gs.transport, design, contrast = contr.mat[,2], nrot = 1000000)
roast3.gs.transport_ps1VSwt

roast3.gs.transport_ps1isoVSps1 <- roast(v3, idx.gs.transport, design, contrast = contr.mat[,1], nrot = 1000000)
roast3.gs.transport_ps1isoVSps1

mycol3 <- colorpanel(1000, "green", "black", "red")
heatmap.2(v3$E[row.names(v3$E) %in% gs.transport_Mm, c(7:9,1:3,4:6)], scale = "row", Colv = "dendrogram", dendrogram = "none",
          col = mycol3, trace = "none", density.info = "none")

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.transport_Mm, c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none")

heatmap.2(lcpm.norm3_5_c[row.names(lcpm.norm3_5_c) %in% gs.transport_Mm, c(7:9,1:3,4:6)],
          scale = "row", Colv = "dendrogram", dendrogram = "none",
          col = mycol4, trace = "none", density.info = "none")

heatmap.2(efit3.contr$t[row.names(efit3.contr$t) %in% gs.transport_Mm, c(2,1)],
          Colv = "dendrogram", dendrogram = "none",
          col = mycol, trace = "none", density.info = "none",
          lhei = c(2,10), margins = c(8,6), cexCol = 1)  ## USE THIS

## Histograms
barplot(matrix(c(res3.ps1isoVSps1$logFC[res3.ps1isoVSps1$genes %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb")],
                 res3.ps1VSwt$genes[res3.ps1VSwt$genes %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb")]), nrow = 2, ncol = 5),
        names.arg = c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb"),
        beside = TRUE)

barplot(res3.ps1isoVSps1$logFC[res3.ps1isoVSps1$genes %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb")],
        names.arg = res3.ps1isoVSps1$genes[res3.ps1isoVSps1$genes %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb")],
        ylim = c(-0.2,0.4))


## PCA
par(mar=c(4,4,4,4), mfrow=c(1,1)) # mar = margin
scatterplot3d(pca3_c$x[, 1:3], angle = 40, pch = pch[group], color = colors[group], grid = F, box = F,
              xlab=paste("PC1, ", round(pca.propVar3_c[1], 2), "%"),
              ylab=paste("PC2, ", round(pca.propVar3_c[2], 2), "%"),
              zlab=paste("PC3, ", round(pca.propVar3_c[3], 2), "%"),
              main = "Principal Component Analysis (3-Dimensional)")
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca3$x[,1:3], grid = c("xy", "xz", "yz"))
legend("topright", legend=c("PS1KO","PS1KO+ISO","WT","WT+ISO"), pch=pch, col=colors)
dev.print(pdf, "ms_figures/Rplot_PCA_3d_laneCor.pdf", height=6, width=8)

## Barcode plots
par(mfrow=c(1,2))
barcodeplot(efit3.contr$t[,2], index = idx.all.clcn.msigdb$GO_CHLORIDE_TRANSPORT, main = "", 
            labels = c("UP in PS1KO", "DOWN in PS1KO")) # 	GO:0006821
title(main = "CHLORIDE TRANSPORT (GO:0006821) \nPS1KO vs WT")
barcodeplot(efit3.contr$t[,1], index = idx.all.clcn.msigdb$GO_CHLORIDE_TRANSPORT, main = "", 
            labels = c("UP in PS1KO+ISO", "DOWN in PS1KO+ISO")) # 	GO:0006821
title(main = "CHLORIDE TRANSPORT (GO:0006821) \nPS1KO+ISO vs PS1KO")
dev.print(pdf, "GSEA/Rplot_barcode3_Cl.Trans_ps1VSwt.ps1isoVSps1.pdf", height=6, width=12)  # 	GO:0006821


par(mfrow=c(1,2))
barcodeplot(efit3.contr$t[,2], index = idx3.all.gs110$All.trasport.genes.232., main = "", 
            labels = c("UP in PS1KO", "DOWN in PS1KO"))
title(main = "All Transport Genes (n=232) \nPS1KO vs WT")
barcodeplot(efit3.contr$t[,1], index = idx3.all.gs110$All.trasport.genes.232., main = "", 
            labels = c("UP in PS1KO+ISO", "DOWN in PS1KO+ISO"))
title(main = "All Transport Genes (n=232) \nPS1KO+ISO vs PS1KO")
dev.print(pdf, "GSEA/Rplot_barcode3_AllTransport_ps1VSwt_ps1isoVSwt.pdf", height=6, width=12) 

## For prism
lcpm.norm3_5[rownames(lcpm.norm3_5) %in% c("Lamp1","Mfsd8","Npc1","Pqlc2","Slc29a3","Slc36a1","Slc38a9","Tpcn2"), ]

lcpm.norm3_5[rownames(lcpm.norm3_5) %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb"), ]

lcpm.norm3_5[rownames(lcpm.norm3_5) %in% geneSets110_Mm$Curated.Induction.23., ]


cpm.norm3 <- cpm(y.norm3)
head(cpm.norm3)

cpm.norm3[rownames(cpm.norm3) %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb"), ]
cpm.norm3[rownames(cpm.norm3) %in% geneSets110_Mm$Curated.Induction.23., ]

res3.ps1VSwt[res3.ps1VSwt$genes %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb"), ]
res3.ps1isoVSps1[res3.ps1isoVSps1$genes %in% c("Mitf","Foxo1","Foxo3","Tfe3","Tfeb"), ]

res3.ps1VSwt[res3.ps1VSwt$genes %in% c("Lamp1","Mfsd8","Npc1","Pqlc2","Slc29a3","Slc36a1","Slc38a9","Tpcn2"), ]
res3.ps1isoVSps1[res3.ps1isoVSps1$genes %in% c("Lamp1","Mfsd8","Npc1","Pqlc2","Slc29a3","Slc36a1","Slc38a9","Tpcn2"), ]

res3.ps1VSwt[res3.ps1VSwt$genes %in% c('Ostm1','Clcn7'), ]
res3.ps1isoVSps1[res3.ps1isoVSps1$genes %in% c('Ostm1','Clcn7'), ]

res3.ps1VSwt[res3.ps1VSwt$genes %in% geneSets110_Mm$Curated.Induction.23., ]
res3.ps1isoVSps1[res3.ps1isoVSps1$genes %in% geneSets110_Mm$Curated.Induction.23., ]


## Volcano plots
library(calibrate)
par(mfrow=c(1,2))

#with(res3.ps1isoVSwt, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-13,15), ylim=c(0,26)))
#with(subset(res3.ps1isoVSwt, res3.ps1isoVSwt$genes %in%
#              c("Cpq","Ctsa","Ctsb","Ctsc","Ctsd","Ctsf","Ctsh","Ctsk","Ctsl","Ctso","Ctss","Ctsw","Ctsz","Ncstn","Prcp","Scpep1",
#                "Tpp1","Cstb","Cst3")),
#     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1VSwt, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in%
              c("Cpq","Ctsa","Ctsb","Ctsc","Ctsd","Ctsf","Ctsh","Ctsk","Ctsl","Ctso","Ctss","Ctsw","Ctsz","Ncstn","Prcp","Scpep1",
                "Tpp1","Cstb","Cst3")),
     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot")) #xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in%
              c("Cpq","Ctsa","Ctsb","Ctsc","Ctsd","Ctsf","Ctsh","Ctsk","Ctsl","Ctso","Ctss","Ctsw","Ctsz","Ncstn","Prcp","Scpep1",
                "Tpp1","Cstb","Cst3")),
     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))

dev.print(pdf, "Rplot_volcano3_LysoGenes-labs_ps1VSwt-ps1isoVSps1.pdf", height=6, width=12)

with(res3.ps1VSwt, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% geneSets110_Mm$Lysosomal.proteolysis.23.),
     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot")) #xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% geneSets110_Mm$Lysosomal.proteolysis.23.),
     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))


with(res3.ps1VSwt, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in%
              c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")),
     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(P.Value), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in%
              c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")),
     textxy(logFC, -log10(P.Value), labs = genes, cex = 0.5, offset = 0.6, col="red"))

## Use -log10(adj.P.val)
with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in%
              c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in%
              c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))


with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-2.5,2.5), ylim=c(0,12.5)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in%
              c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot", xlim=c(-2.5,2.5), ylim=c(0,4)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in%
              c("Clcn1","Clcn2","Clcn3","Clcn6","Clcn7","Clcn5","Clcn6","Kcne2","Lamp1","Lamp2","Limp2","Mcoln1","Mfsd1","Mfsd8",
                "Npc1","Npc2","Ostm1","Pqlc2","Ramp2","Slc11a2","Slc15a3","Slc29a3","Slc36a1","Slc38a9","Tcirg1","Tmem9","Tpcn1","Tpcn2")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

## For vesicle mediate transport genes (Matteo's pipeline)
with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% geneSets110_Mm$vesicle.mediated.transport.16.),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% geneSets110_Mm$vesicle.mediated.transport.16.),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

## For ER to Golgi Anterograde Transport from GO: 0006888
with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% gs.ErToGolgi.GO_Mm),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% gs.ErToGolgi.GO_Mm),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))


with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% c("Arf1","Sorl1","Stx18","Tbc1d20")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% c("Arf1","Sorl1","Stx18","Tbc1d20")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

## For Golgi to PM Protein Transport from GO: 00430001
with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% gs.GolgiToPM.GO_Mm),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% gs.GolgiToPM.GO_Mm),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

## GO_INTRACELLULAR_CALCIUM_ACTIVATED_CHLORIDE_CHANNEL_ACTIVITY 
with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% gs.clcn.msigdb_Mm$GO_INTRACELLULAR_CALCIUM_ACTIVATED_CHLORIDE_CHANNEL_ACTIVITY),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% gs.clcn.msigdb_Mm$GO_INTRACELLULAR_CALCIUM_ACTIVATED_CHLORIDE_CHANNEL_ACTIVITY),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))



## For Snx genes Snx7, 12, 14, 19, 21
with(res3.ps1VSwt, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1VSwt, res3.ps1VSwt$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 1, offset = 0.6, col="red"))

with(res3.ps1isoVSps1, plot(logFC, -log10(adj.P.Val), pch=20, cex=0.35, col="grey", main="Volcano plot"))# xlim=c(-13,15), ylim=c(0,26)))
with(subset(res3.ps1isoVSps1, res3.ps1isoVSps1$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21")),
     textxy(logFC, -log10(adj.P.Val), labs = genes, cex = 0.5, offset = 0.6, col="red"))
## Using ggplot2
library(ggplot2)
library(dplyr)
library(ggrepel)
res3.ps1VSwt_volc = mutate(res3.ps1VSwt, SIGNIFICANCE=ifelse(res3.ps1VSwt$adj.P.Val<0.05 & res3.ps1VSwt$logFC<0, "Downregulated", 
                                                             ifelse(res3.ps1VSwt$adj.P.Val<0.05 & res3.ps1VSwt$logFC>0, "Upregulated", "No Change")))
head(res3.ps1VSwt_volc)

p_ps1VSwt = ggplot(res3.ps1VSwt_volc, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=SIGNIFICANCE)) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_point(data = res3.ps1VSwt_volc[res3.ps1VSwt_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], color = "white") +

p_ps1VSwt

p_ps1VSwt = ggplot(res3.ps1VSwt_volc, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=SIGNIFICANCE), size=0.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_point(data = res3.ps1VSwt_volc[res3.ps1VSwt_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], color = "white") +
  geom_text_repel(data = res3.ps1VSwt_volc[res3.ps1VSwt_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], aes(label=genes), box.padding = 1.5, color="red")
p_ps1VSwt

p2_ps1VSwt = ggplot(res3.ps1VSwt_volc, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(col=SIGNIFICANCE), size=0.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_point(data = res3.ps1VSwt_volc[res3.ps1VSwt_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], color = "white") +
  geom_text_repel(data = res3.ps1VSwt_volc[res3.ps1VSwt_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], 
                  aes(label=genes), 
                  box.padding = 1.5, 
                  color="red",
                  arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"), 
                  force = 10, 
                  nudge_x = -10, 
                  direction = "x") +
  theme_classic()
p2_ps1VSwt


res3.ps1isoVSps1_volc = mutate(res3.ps1isoVSps1, SIGNIFICANCE=ifelse(res3.ps1isoVSps1$adj.P.Val<0.05 & res3.ps1isoVSps1$logFC<0, "Downregulated", 
                                                             ifelse(res3.ps1isoVSps1$adj.P.Val<0.05 & res3.ps1isoVSps1$logFC>0, "Upregulated", "No Change")))
head(res3.ps1isoVSps1_volc)
p2_ps1isoVSps1 = ggplot(res3.ps1isoVSps1_volc, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=SIGNIFICANCE), size=0.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_point(data = res3.ps1isoVSps1_volc[res3.ps1isoVSps1_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], color = "white") +
  geom_text_repel(data = res3.ps1isoVSps1_volc[res3.ps1isoVSps1_volc$genes %in% c("Snx7","Snx12","Snx14","Snx19","Snx21"), ], 
                  aes(label=genes), 
                  box.padding = 1.5, 
                  color="red",
                  arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"), 
                  force = 5, 
                  nudge_x = 1, 
                  direction = "x") +
  theme_classic()
p2_ps1isoVSps1

## Supplemental Table
res3.ps1VSwt_protTrans <- res3.ps1VSwt[row.names(res3.ps1VSwt) %in% geneSets110_Mm$protein.transport.40., ]
res3.ps1VSwt_protTrans
#res3.ps1VSwt_protTrans$genes <- res3.ps1VSwt_protTrans$genes[order(geneSets110_Mm$protein.transport.40.[!is.na(geneSets110_Mm$protein.transport.40.)])]
res3.ps1VSwt_protTrans <- res3.ps1VSwt_protTrans[match(geneSets110_Mm$protein.transport.40., row.names(res3.ps1VSwt_protTrans), nomatch=0), ]
res3.ps1VSwt_protTrans
write.table(res3.ps1VSwt_protTrans, file = "ms_figures/res3.ps1VSwt_protTrans.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


res3.ps1isoVSps1_protTrans <- res3.ps1isoVSps1[row.names(res3.ps1isoVSps1) %in% geneSets110_Mm$protein.transport.40., ]
res3.ps1isoVSps1_protTrans
res3.ps1isoVSps1_protTrans <- res3.ps1isoVSps1_protTrans[match(geneSets110_Mm$protein.transport.40., row.names(res3.ps1isoVSps1_protTrans))[1:40], ]
res3.ps1isoVSps1_protTrans
write.table(res3.ps1isoVSps1_protTrans, file = "ms_figures/res3.ps1isoVSps1_protTrans.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


res3.ps1VSwt_endoToLys <- res3.ps1VSwt[row.names(res3.ps1VSwt) %in% geneSets110_Mm$endosome.to.lysosome.transport.10., ]
res3.ps1VSwt_endoToLys
res3.ps1VSwt_endoToLys <- res3.ps1VSwt_endoToLys[match(geneSets110_Mm$endosome.to.lysosome.transport.10., row.names(res3.ps1VSwt_endoToLys))[1:10], ]
res3.ps1VSwt_endoToLys
write.table(res3.ps1VSwt_endoToLys, file = "ms_figures/res3.ps1VSwt_endoToLys.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSps1_endoToLys <- res3.ps1isoVSps1[row.names(res3.ps1isoVSps1) %in% geneSets110_Mm$endosome.to.lysosome.transport.10., ]
res3.ps1isoVSps1_endoToLys
res3.ps1isoVSps1_endoToLys <- res3.ps1isoVSps1_endoToLys[match(geneSets110_Mm$endosome.to.lysosome.transport.10., row.names(res3.ps1isoVSps1_endoToLys))[1:10], ]
res3.ps1isoVSps1_endoToLys
write.table(res3.ps1isoVSps1_endoToLys, file = "ms_figures/res3.ps1isoVSps1_endoToLys.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


res3.ps1VSwt_vesMedTrans <- res3.ps1VSwt[row.names(res3.ps1VSwt) %in% geneSets110_Mm$vesicle.mediated.transport.16., ]
res3.ps1VSwt_vesMedTrans
res3.ps1VSwt_vesMedTrans <- res3.ps1VSwt_vesMedTrans[match(geneSets110_Mm$vesicle.mediated.transport.16., row.names(res3.ps1VSwt_vesMedTrans))[1:16], ]
res3.ps1VSwt_vesMedTrans
write.table(res3.ps1VSwt_vesMedTrans, file = "ms_figures/res3.ps1VSwt_vesMedTrans.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSps1_vesMedTrans <- res3.ps1isoVSps1[row.names(res3.ps1isoVSps1) %in% geneSets110_Mm$vesicle.mediated.transport.16., ]
res3.ps1isoVSps1_vesMedTrans
res3.ps1isoVSps1_vesMedTrans <- res3.ps1isoVSps1_vesMedTrans[match(geneSets110_Mm$vesicle.mediated.transport.16., row.names(res3.ps1isoVSps1_vesMedTrans))[1:16], ]
res3.ps1isoVSps1_vesMedTrans
write.table(res3.ps1isoVSps1_vesMedTrans, file = "ms_figures/res3.ps1isoVSps1_vesMedTrans.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


res3.ps1VSwt_allTrans <- res3.ps1VSwt[row.names(res3.ps1VSwt) %in% geneSets110_Mm$All.trasport.genes.232., ]
head(res3.ps1VSwt_allTrans)
res3.ps1VSwt_allTrans <- res3.ps1VSwt_allTrans[match(geneSets110_Mm$All.trasport.genes.232., row.names(res3.ps1VSwt_allTrans))[1:232], ]
head(res3.ps1VSwt_allTrans)
write.table(res3.ps1VSwt_allTrans, file = "ms_figures/res3.ps1VSwt_allTrans.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSps1_allTrans <- res3.ps1isoVSps1[row.names(res3.ps1isoVSps1) %in% geneSets110_Mm$All.trasport.genes.232., ]
head(res3.ps1isoVSps1_allTrans)
res3.ps1isoVSps1_allTrans <- res3.ps1isoVSps1_allTrans[match(geneSets110_Mm$All.trasport.genes.232., row.names(res3.ps1isoVSps1_allTrans))[1:232], ]
head(res3.ps1isoVSps1_allTrans)
write.table(res3.ps1isoVSps1_allTrans, file = "ms_figures/res3.ps1isoVSps1_allTrans.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


res3.ps1VSwt_golgiRecy.Sec <- res3.ps1VSwt[row.names(res3.ps1VSwt) %in% geneSets110_Mm$X.Golgi..recycling.comp..And.secretory..24., ]
res3.ps1VSwt_golgiRecy.Sec
res3.ps1VSwt_golgiRecy.Sec <- res3.ps1VSwt_golgiRecy.Sec[match(geneSets110_Mm$X.Golgi..recycling.comp..And.secretory..24., row.names(res3.ps1VSwt_golgiRecy.Sec))[1:24], ]
res3.ps1VSwt_golgiRecy.Sec
write.table(res3.ps1VSwt_golgiRecy.Sec, file = "ms_figures/res3.ps1VSwt_golgiRecy.Sec.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)

res3.ps1isoVSps1_golgiRecy.Sec <- res3.ps1isoVSps1[row.names(res3.ps1isoVSps1) %in% geneSets110_Mm$X.Golgi..recycling.comp..And.secretory..24., ]
res3.ps1isoVSps1_golgiRecy.Sec
res3.ps1isoVSps1_golgiRecy.Sec <- res3.ps1isoVSps1_golgiRecy.Sec[match(geneSets110_Mm$X.Golgi..recycling.comp..And.secretory..24., row.names(res3.ps1isoVSps1_golgiRecy.Sec))[1:24], ]
res3.ps1isoVSps1_golgiRecy.Sec
write.table(res3.ps1isoVSps1_golgiRecy.Sec, file = "ms_figures/res3.ps1isoVSps1_golgiRecy.Sec.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)


## Tableau
# Preparing df for GO_CHLORIDE_TRANSPORT gene set
mroast3.clcn.msigdb_ps1isoVSps1
head(gs.clcn.msigdb_Mm)

df.clcnTran_ps1isoVSps1 <- efit3.contr$coefficients[rownames(efit3.contr) %in% gs.clcn.msigdb_Mm$GO_CHLORIDE_TRANSPORT, 1]
length(df.clcnTran_ps1isoVSps1) # 46
head(df.clcnTran_ps1isoVSps1)
df.clcnTran_ps1isoVSps1 <- as.data.frame(df.clcnTran_ps1isoVSps1)
head(df.clcnTran_ps1isoVSps1)
colnames(df.clcnTran_ps1isoVSps1) <- "log2FC"
head(df.clcnTran_ps1isoVSps1)

# Adding t-value as a column
df.clcnTran_ps1isoVSps1$t.value <- efit3.contr$t[match(rownames(df.clcnTran_ps1isoVSps1), rownames(efit3.contr$t)), 1]
head(df.clcnTran_ps1isoVSps1)
tail(df.clcnTran_ps1isoVSps1)
dim(df.clcnTran_ps1isoVSps1) # 46  2

write.table(df.clcnTran_ps1isoVSps1, file = "Tableau/df.clcnTran_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

df.clcnTran_ps1VSwt <- efit3.contr$coefficients[rownames(efit3.contr) %in% gs.clcn.msigdb_Mm$GO_CHLORIDE_TRANSPORT, 2]
length(df.clcnTran_ps1VSwt) # 46
head(df.clcnTran_ps1VSwt)
df.clcnTran_ps1VSwt <- as.data.frame(df.clcnTran_ps1VSwt)
head(df.clcnTran_ps1VSwt)
colnames(df.clcnTran_ps1VSwt) <- "log2FC"
head(df.clcnTran_ps1VSwt)

# Adding t-value as a column
df.clcnTran_ps1VSwt$t.value <- efit3.contr$t[match(rownames(df.clcnTran_ps1VSwt), rownames(efit3.contr$t)), 2]
head(df.clcnTran_ps1VSwt)
tail(df.clcnTran_ps1VSwt)
dim(df.clcnTran_ps1VSwt) # 46  2

write.table(df.clcnTran_ps1VSwt, file = "Tableau/df.clcnTran_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

## For GEO
head(cpm.norm3)
tail(cpm.norm3)
dim(cpm.norm3) # 15658    12
write.table(cpm.norm3, file = "~/Bioinformatics/PS1KO Iso/MS/GEO/countsPerMillion_filt_TMMnorm.txt", sep = "\t", 
            row.names = T, col.names = NA, quote = F)



# Preparing df for Curated induction gene set

df.curInd_ps1isoVSps1 <- efit3.contr$coefficients[rownames(efit3.contr) %in% geneSets110_Mm$Curated.Induction.23., 1]
length(df.curInd_ps1isoVSps1) # 23
head(df.curInd_ps1isoVSps1)
df.curInd_ps1isoVSps1 <- as.data.frame(df.curInd_ps1isoVSps1)
head(df.curInd_ps1isoVSps1)
colnames(df.curInd_ps1isoVSps1) <- "log2FC"
head(df.curInd_ps1isoVSps1)

# Adding t-value as a column
df.curInd_ps1isoVSps1$t.value <- efit3.contr$t[match(rownames(df.curInd_ps1isoVSps1), rownames(efit3.contr$t)), 1]
head(df.curInd_ps1isoVSps1)
tail(df.curInd_ps1isoVSps1)
dim(df.curInd_ps1isoVSps1) # 23  2

write.table(df.curInd_ps1isoVSps1, file = "Tableau/df.curInd_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

df.curInd_ps1VSwt <- efit3.contr$coefficients[rownames(efit3.contr) %in% geneSets110_Mm$Curated.Induction.23., 2]
length(df.curInd_ps1VSwt) # 23
head(df.curInd_ps1VSwt)
df.curInd_ps1VSwt <- as.data.frame(df.curInd_ps1VSwt)
head(df.curInd_ps1VSwt)
colnames(df.curInd_ps1VSwt) <- "log2FC"
head(df.curInd_ps1VSwt)

# Adding t-value as a column
df.curInd_ps1VSwt$t.value <- efit3.contr$t[match(rownames(df.curInd_ps1VSwt), rownames(efit3.contr$t)), 2]
head(df.curInd_ps1VSwt)
tail(df.curInd_ps1VSwt)
dim(df.curInd_ps1VSwt) # 23  2

write.table(df.curInd_ps1VSwt, file = "Tableau/df.curInd_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)


# Preparing df for "Curated induction "All Transport Genes" gene set

df.allTrans_ps1isoVSps1 <- efit3.contr$coefficients[rownames(efit3.contr) %in% geneSets110_Mm$All.trasport.genes.232., 1]
length(df.allTrans_ps1isoVSps1) # 194
head(df.allTrans_ps1isoVSps1)
df.allTrans_ps1isoVSps1<- as.data.frame(df.allTrans_ps1isoVSps1)
head(df.allTrans_ps1isoVSps1)
colnames(df.allTrans_ps1isoVSps1) <- "log2FC"
head(df.allTrans_ps1isoVSps1)

# Adding t-value as a column
df.allTrans_ps1isoVSps1$t.value<- efit3.contr$t[match(rownames(df.allTrans_ps1isoVSps1), rownames(efit3.contr$t)), 1]
head(df.allTrans_ps1isoVSps1)
tail(df.allTrans_ps1isoVSps1)
dim(df.allTrans_ps1isoVSps1) # 194  2

write.table(df.allTrans_ps1isoVSps1, file = "Tableau/df.allTrans_ps1isoVSps1.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

df.allTrans_ps1VSwt <- efit3.contr$coefficients[rownames(efit3.contr) %in% geneSets110_Mm$All.trasport.genes.232., 2]
length(df.allTrans_ps1VSwt) # 194
head(df.allTrans_ps1VSwt)
df.allTrans_ps1VSwt<- as.data.frame(df.allTrans_ps1VSwt)
head(df.allTrans_ps1VSwt)
colnames(df.allTrans_ps1VSwt) <- "log2FC"
head(df.allTrans_ps1VSwt)

# Adding t-value as a column
df.allTrans_ps1VSwt$t.value<- efit3.contr$t[match(rownames(df.allTrans_ps1VSwt), rownames(efit3.contr$t)), 2]
head(df.allTrans_ps1VSwt)
tail(df.allTrans_ps1VSwt)
dim(df.allTrans_ps1VSwt) # 194  2

write.table(df.allTrans_ps1VSwt, file = "Tableau/df.allTrans_ps1VSwt.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)

############################################################################################################################
library(biomaRt)
listMarts()
ensembl <- useMart("ensembl")
head(listDatasets(ensembl))
head(searchDatasets(mart = ensembl, pattern = "musculus"))

ensembl_dataSet <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
head(searchAttributes(ensembl_dataSet, pattern = "entrez"))
head(searchFilters(ensembl_dataSet, pattern = "symbol"))

df_genes.detected <- getBM(attributes = c("hgnc_symbol","entrezgene"), 
                           filters = "mgi_symbol", 
                           values = row.names(genes.detected), 
                           mart = ensembl_dataSet)
head(df_genes.detected)
dim(df_genes.detected) # 17459     2
length(unique(df_genes.detected$entrezgene)) # 17389


df_genes.detected3 <- getBM(attributes = c("mgi_symbol","entrezgene", "ensembl_gene_id"), 
                           filters = "mgi_symbol", 
                           values = row.names(genes.detected), 
                           mart = ensembl_dataSet)
head(df_genes.detected3)
dim(df_genes.detected3) # 29045     3
length(unique(df_genes.detected3$entrezgene)) # 17389
length(na.omit(unique(df_genes.detected3$entrezgene))) # 17388

length(na.omit(unique(df_genes.detected3$ensembl_gene_id))) # 28966
write.table(na.omit(unique(df_genes.detected3$ensembl_gene_id)), file = "GO/allGenesDetectedinExpt_ENSEMBL_NonZeroReadCts.txt", 
            row.names = F, col.names = F, quote = F)


library(AnnotationDbi)
library(org.Mm.eg.db)

df_genes.detected2 <- select(org.Mm.eg.db, keys = row.names(genes.detected), keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID", "ENSEMBL"))
                      # 'select()' returned 1:many mapping between keys and columns
head(df_genes.detected2)
length(unique(df_genes.detected2$ENTREZID)) # 24167
length(na.omit(unique(df_genes.detected2$ENTREZID))) # 24166
write.table(na.omit(unique(df_genes.detected2$ENTREZID)), file = "GO/allGenesDetectedinExpt_ENTREZ_NonZeroReadCts.txt", 
            row.names = F, col.names = F, quote = F)

length(unique(df_genes.detected2$ENSEMBL)) # 18719
length(na.omit(unique(df_genes.detected2$ENSEMBL))) # 18718
#write.table(na.omit(unique(df_genes.detected2$ENTREZID)), file = "GO/allGenesDetectedinExpt_ENTREZ_NonZeroReadCts.txt", 
#            row.names = F, col.names = F, quote = F)



############################################################################################################################



############################################################################################################################
## Scatter plots/correlation
library(psych)
par(mfrow=c(1,1))
head(iris)
tail(iris)
pairs.panels(efit3.contr$coefficients,
             method = "pearson",
             hist.col = "lightgreen",
             density = TRUE,
             ellipses = TRUE)

### Perform enrichment analysis for Biological Process (BP)
## Note that the argument is keytype instead of keyType in Bioconductor 3.5
library(clusterProfiler)
library(org.Mm.eg.db)
enrichGO_ps1VSwt<- enrichGO(
  gene = rownames(res3.ps1VSwt)[res3.ps1VSwt$adj.P.Val < 0.05],
  OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
  universe = rownames(res3.ps1VSwt))
head(enrichGO_ps1VSwt,20)


head(gofilter(enrichGO_ps1VSwt, level = 4), 20)
gofilter(enrichGO_ps1VSwt, level = 5)[1:30, ]
gofilter(enrichGO_ps1VSwt, level = 3)[1:20, ]
###
enrichGO_ps1isoVSps1 <- enrichGO(
  gene = rownames(res3.ps1isoVSps1)[res3.ps1isoVSps1$adj.P.Val < 0.001],
  OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
  universe = rownames(res3.ps1isoVSps1))
head(enrichGO_ps1isoVSps1,20)

head(gofilter(enrichGO_ps1isoVSps1, level = 4), 20)
gofilter(enrichGO_ps1isoVSps1, level = 3)[1:30, ]


## Visualize enrichment results
dotplot(enrichGO_ps1isoVSps1, showCategory=20, font.size = 10)
barplot(enrichGO_ps1isoVSps1, showCategory = 20, font.size = 7)

#install("RDAVIDWebService")
library(RDAVIDWebService)
david<-DAVIDWebService(email="sdarji@nki.rfmh.org",
                       url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
david

getTimeOut(david) # [1] 30000

##Increment the time out interval

setTimeOut(david, 500000)

getTimeOut(david) # [1] 500000

getHttpProtocolVersion(david)
#setHttpProtocolVersion(david, "HTTP/1.0")
#getHttpProtocolVersion(david)

enrich_david <- enrichDAVID(gene = rownames(res3.1.ps1isoVSps1)[res3.1.ps1isoVSps1$adj.P.Val < 0.0001],
                            species = "Mm", idType = "ENTREZ_GENE_ID",
                            pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
                            universe = rownames(res3.1.ps1isoVSps1),
                            david.user = "sdarji@nki.rfmh.org") # certificate error
#authenticate("sdarji@nki.rfmh.org")
enrich_kegg <- enrichKEGG(gene = rownames(res3.ps1isoVSps1)[res3.ps1isoVSps1$adj.P.Val < 0.001],
                          organism = "mmu", keyType = "ncbi_geneid",
                          pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
                          universe = rownames(res3.ps1isoVSps1),
                          use_internal_data = FALSE)




y.dge$counts[rownames(y.dge) %in% c("Ctsk"), ]
#    ps1_1    ps1_2    ps1_3 ps1iso_1 ps1iso_2 ps1iso_3     wt_1     wt_2     wt_3  wtiso_1  wtiso_2  wtiso_3
#        7        9        9       11        3       11       26       20       24       28       35       29

