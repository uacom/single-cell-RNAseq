qiu = DGEList(counts = fc$counts, genes = fc$annotation[, c("GeneID", "Length")])

qiu$samples

# simplify samplenames and include group as factor in qiu DGEList object
samplenames = c("FAD1190", "FAD1193", "FAD1219", "wt1191", "wt1213", "wt1214")
colnames(qiu) = samplenames
group = as.factor(c("FAD", "FAD", "FAD", "wt", "wt", "wt"))
qiu$samples$group = group

# organizing gene annotations
# add second dataframe of genes to DGEList object 
GeneID = fc$annotation$GeneID
Chr = fc$annotation$Chr
genes = data.frame(GeneID, Chr)
head(genes)

# remove duplicated gene IDs by keeping only first occurrence of each gene ID
genes = genes[!duplicated(genes$GeneID),]
qiu$genes = genes
head(genes)

rm(GeneID, Chr)

dim(qiu)
#convert raw counts to CPM and lCPM
cpm <- cpm(qiu)
lcpm <- cpm(qiu, log=TRUE)

# L is avg library size (which is in the 20 millions) in millions
L = mean(qiu$samples$lib.size)*1e-6
M = median(qiu$samples$lib.size) * 1e-6

c(L, M)
summary(lcpm)

# graph of raw data
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(qiu)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.62), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# Remove lowly expressed genes

# find the number of genes that are not expressed in all samples
table(rowSums(qiu$counts==0)==6)
# automatically filter genes using filterByExpr function in edgeR package
keep.exprs = filterByExpr(qiu, group = group)
qiu = qiu[keep.exprs,, keep.lib.sizes = FALSE]
dim(qiu)

# graph of filtered data

lcpm <- cpm(qiu, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.62), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# Normalize gene expression distributions vs trimmed mean of M-values
qiu = calcNormFactors(qiu, method = "TMM")
qiu$samples$norm.factors

# Unsupervised clustering of samples: MDS plot using limma
lcpm = cpm(qiu, log=T)
par(mfrow=c(1,1))
col.group = group
levels(col.group)=brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="Sample groups")

# Differential expression analysis

#create design matrix and contrasts

# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
# from what I've read, when creating the design matrix, we can in/exclude the intercept term. If group is a factor (explanatory variable that is categorical in nature)
# then the two design matrices are different, if group is a covariate (explanatory variable that is quantitative) then they are the same

# I've set up the design matrix without a y-int for group
# the design matrix is parameterized for a means model: first col of the design matrix is parameterized for the mean expression of the FAD group (1 is present when the associated sample belongs to the FAD group, 0 if otherwise)
design = model.matrix(~0+group)
colnames(design) = gsub("group", "", colnames(design))
design

# set up contrasts
# subtract first parameter estimate (mean expression for FAD) from the second parameter estimate (mean expression for wt)
contr.matrix = makeContrasts(
  wt_vs_FAD=wt-FAD,
  levels = colnames(design)
)

contr.matrix

# remove heteroscedascity from count data

par(mfrow=c(1,2))
# show mean-variance relationship estimated by voom
v = voom(qiu, design, plot = T)
v

# fitting linear models for comparisons of interest
# use limma package's lmFit and contrasts.fit functions to conduct linear modeling
# Use eBayes function to conduct empirical Bayes moderation to obtain more precise estimates of gene-wise variability
vfit = lmFit(v, design) 
vfit = contrasts.fit(vfit, contrasts=contr.matrix)
efit = eBayes(vfit)
# plot after eBayes
plotSA(efit, main= "Final model: Mean-variance trend")

# Examining the number of DE genes
summary(decideTests(efit, p.value=0.05))


