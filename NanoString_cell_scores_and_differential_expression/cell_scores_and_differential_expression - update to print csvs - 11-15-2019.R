######### citn07 analysis: look at immune cell gene expressions vs. time and cohort, and perform differential expression analyses

rm(list=ls())
library(ggplot2)
library(pheatmap)
source("source_volcplot.r")

# coloring function
fadecols = function(cols,fade=.5)
{
  fcols = c()
  for(i in 1:length(cols))
  {
    tmp = as.vector(col2rgb(cols[i])/256)
    fcols[i] = rgb(tmp[1],tmp[2],tmp[3],fade)
  }
  
  return(fcols)
}

#### load normalized data and sample annotation: ####
load("data/CITN07 data ready for analysis October2018.RData")

# use the log-transformed data:
norm = log2norm; rm(log2norm)  #log2(pmax(norm, 1))
#### screen out low-signal genes:
low.expression = names(which(max.counts < 25))
norm = norm[, setdiff(colnames(norm), low.expression)]


##### reformat the gene annotation matrix:
gannot = gannot[match(colnames(norm), rownames(gannot)), ]
gannot = replace(gannot, is.na(gannot), 0)
rownames(gannot) = colnames(norm)

# load additional gene sets:
genesets = read.csv("data/gene sets.csv", stringsAsFactors = F)
gannot2 = matrix(0, ncol(norm), length(unique(genesets$sig)), dimnames = list(colnames(norm), unique(genesets$sig)))
# fill it in:
for (m in unique(genesets$sig)){
  modgenes = genesets$gene[genesets$sig == m]
  gannot2[is.element(rownames(gannot2), modgenes), m] = 1
}
# remove modules with too few genes to be useful:
gannot2 = gannot2[, colSums(gannot2) >= 3]

gannot = cbind(gannot, gannot2)
# add mini-modules split out from big gene sets:
genesets.small = read.csv("data/gene sets - smaller.csv", stringsAsFactors = F)
genesets.small = genesets.small[, !is.na(genesets.small[1, ])]
gannot.small = matrix(0, nrow(gannot), ncol(genesets.small), 
                       dimnames = list(rownames(gannot), colnames(genesets.small)))
for (set in colnames(genesets.small)) {
  gannot.small[setdiff(genesets.small[, set], c(NA,"")), set] = 1
}
gannot = cbind(gannot, gannot.small)

######################################
#### calculate cell type signatures: #### 
######################################
cellmarkers = read.csv("data/cell type markers pbmc.csv",row.names = 1)
# function to calc cell scores from data and marker list:
calc.cellscores = function(e,markers)
{
  cellscores = c()
  cells = unique(markers[,"Cell.Type"])
  for(cell in cells)
  {
    tempgenes = intersect(colnames(e),rownames(markers)[markers[,"Cell.Type"]==cell])
    lostgenes = setdiff(rownames(markers)[markers[,"Cell.Type"]==cell],colnames(e))
    if(length(lostgenes)>0){print(paste0(cell," genes missing:"));print(lostgenes)}
    if(length(tempgenes)>0)
    {
      cellscores = cbind(cellscores,rowMeans(e[,tempgenes,drop=F]))
      colnames(cellscores)[ncol(cellscores)]=cell
    }
  }
  return(cellscores)
}
cellscores = calc.cellscores(norm,cellmarkers)



######################################
#### further transformations of data: 
# 1. get mean and se of expression/cell score within each cohort/time
# 2. reformat expression/score matrices as deltas from baselines
# 3. get mean and se of expression/score deltas within each cohort/time
######################################

## combined genes/scores data matrix:
normpluscells = cbind(cellscores, norm)
is.cellscore = c(rep(TRUE, ncol(cellscores)), rep(FALSE, ncol(norm)))

## timepoint/cohort variables:
#timepoints = c("Pre-Tx","C1D22","C3D01","C4D01","FUW04","FUW12")
timepoints = levels(annot$time)
txtimes = paste0(rep(c("cohort1","cohort2"), each=length(timepoints))," ",rep(timepoints,2))
txtimes = setdiff(txtimes, "cohort2 C1D01")
annot$txtime2 = paste0("cohort",annot$cohort," ",annot$time)

# 1. get mean and se of expression/cell score within each time point
# means within each timepoint/cohort
means = matrix(NA,length(txtimes),ncol(normpluscells),dimnames = list(txtimes,colnames(normpluscells)))
for(txtime in txtimes)
{
  means[txtime,] = colMeans(normpluscells[annot$txtime2==txtime,])
}
# ses within each timepoint/cohort:
getse = function(x){return(sd(x,na.rm=T)/sqrt(sum(!is.na(x))))}
ses = matrix(NA, length(txtimes), ncol(normpluscells), dimnames = list(txtimes, colnames(normpluscells)))
for(txtime in txtimes)
{
  ses[txtime,] = apply(normpluscells[annot$txtime2==txtime, ], 2, getse)
}

# 2. define a matrix of each patients' gene expression deltas from baseline:
deltas = normpluscells*NA
for(tp in timepoints)
{
  # only look at the current timepoint's data:
  use = which(annot$time==tp)
  # for each observation, get the mean expression vector of the patient's time0 samples:
  for(i in use)
  {
    pt = annot$patient[i]
    baseline = which((annot$patient==pt)&(annot$time=="baseline"))
    if(length(baseline) > 0){
      deltas[i,] = normpluscells[i,] - colMeans(normpluscells[baseline,,drop=F],na.rm=T)
    }
  }
}

# 3. get the avg (and se) delta from pre-tx within each cohort*time:
meandeltas = matrix(NA, length(txtimes), ncol(normpluscells), dimnames = list(txtimes, colnames(normpluscells)))
for(txtime in txtimes)
{
  meandeltas[txtime,] = colMeans(deltas[annot$txtime2==txtime,],na.rm=T)
}

# and ses of deltas:
sedeltas = matrix(NA,length(txtimes),ncol(normpluscells),dimnames = list(txtimes,colnames(normpluscells)))
for(txtime in txtimes)
{
  sedeltas[txtime,] = apply(deltas[annot$txtime2==txtime,],2,getse)
}

# get indices of cohorts:
inc1 = grepl("cohort1", rownames(means))
inc2 = grepl("cohort2", rownames(means))


######################################
#### trend plot of cells over time:
######################################
cells = setdiff(colnames(cellscores), "CD45")
cellcols = c("darkblue","orange","chartreuse3","purple","blue","forestgreen","cornflowerblue","red","firebrick","grey50","lightblue","pink","green")[1:length(cells)]
names(cellcols)=cells


## plot cell types over time
pdf("manuscript - cells over time v2.pdf", width=9, height = 5)
par(mar = c(6,4.5,.5,1))
layout(matrix(1:3, 1, 3), width = c(4.5,4,2.7))
plot(0,0,xlim=c(1,sum(inc1)),ylim=range(meandeltas[,cells]), xaxt="n",xlab="",
     ylab="Mean log2 fold-change from baseline",main="", cex.lab = 1.8)
abline(h=0,col = "grey20", lty = 3, lwd = 2)
axis(1,1:sum(inc1),timepoints,las=2, cex.axis = 1.5)
for(cell in cells)
{
  lines(1:sum(inc1),meandeltas[inc1,cell],col=cellcols[cell],lwd=2)
  # points for changes sig at p=0.05:
  is.sig = (meandeltas[(inc1),cell] - 2*sedeltas[(inc1),cell] > 0)|(meandeltas[(inc1),cell] + 2*sedeltas[(inc1),cell] < 0)
  points((1:sum(inc1))[is.sig],meandeltas[which(inc1)[is.sig],cell],col=cellcols[cell],pch=16,cex=2)
}
legend("top", legend = "Cohort 1   ", text.col = "orange", cex = 3, bty = "n")

par(mar = c(6,0,.5,1))
plot(0,0,xlim=c(1,sum(inc2)),ylim=range(meandeltas[,cells]),xaxt="n",xlab="",
     ylab="",main="", cex.lab = 1.5, yaxt = "n")
abline(h=0,col = "grey20", lty = 3, lwd = 2)
axis(1,1:sum(inc2),setdiff(timepoints, "C1D01"),las=2, cex.axis = 1.5)
for(cell in cells)
{
  lines(1:sum(inc2),meandeltas[inc2,cell],col=cellcols[cell],lwd=2)
  # points for changes sig at p=0.05:
  is.sig = (meandeltas[(inc2),cell] - 2*sedeltas[(inc2),cell] > 0)|(meandeltas[(inc2),cell] + 2*sedeltas[(inc2),cell] < 0)
  points((1:sum(inc2))[is.sig],meandeltas[which(inc2)[is.sig],cell],col=cellcols[cell],pch=16,cex=2)
}
legend("top", legend = "Cohort 2   ", text.col = "cornflowerblue", cex = 3, bty = "n")

par(mar = c(5.6,0,0,0))
frame()
o.legend = order(meandeltas[3,names(cellcols)], decreasing = T)
legend("center", col = c("black", "black", NA, cellcols[o.legend]), 
       legend = c("p > 0.05", "p < 0.05", "", names(cellcols)[o.legend]), 
       lty = 1, pch = c(NA, 16, NA, rep(NA, length(o.legend))), lwd = 3, bty = "n", cex = 1.5)
dev.off()

# write fig 5e data csv:
write.csv(meandeltas[, cells], file = "figure data - 5e - cell scores' mean deltas from baseline.csv")
write.csv(cellscores, file = "figure data - 5e - cell scores in all samples.csv")



### table of cell type results:
#times = as.character(unique(annot$time))
celltable = cbind(cells,
                  paste0(round(meandeltas["cohort1 C1D08", cells], 2),
                 " (", 
                  round(sedeltas["cohort1 C1D08", cells], 2),
                  ")"),
                  paste0(round(meandeltas["cohort2 C1D08", cells], 2),
                  " (", 
                  round(sedeltas["cohort2 C1D08", cells], 2),
                  ")"))[o.legend, ]
colnames(celltable) = c("Cell score", "Cohort 1 C1D08 log2 fold change from baseline", "Cohort 2 C1D08 log2 fold change from baseline")      
write.csv(celltable, file = "cell scores C1D08 results.csv", row.names = F)

#### perform differential expression analyses ####



######
# model delta from baseline within cohort 1
######

# initialize results matrices: e for estimated log2 fold-change; p for p-values:
e = se = p = matrix(NA,ncol(norm),length(levels(annot$time))-1,
               dimnames=list(colnames(norm),setdiff(levels(annot$time),"baseline")))
# at each timepoint, run the model:
for(time in setdiff(levels(annot$time),"baseline"))
{
  use = (annot$cohort=="1")&(is.element(annot$time,c("baseline",time)))
  for(gene in colnames(norm))
  {
    # paired t-test:
    mod = t.test(deltas[(annot$cohort=="1")&(is.element(annot$time,time)),gene])
    e[gene,time] = mod$estimate  #mod[2,1]
    p[gene,time] = mod$p.value   #mod[2,4]
  }
}
# get fdrs:
f = p*NA
for(i in 1:ncol(p)){f[,i]=p.adjust(p[,i], method = "BH")}

# save results for C1D08:
e2 = e[,"C1D08"]
p2 = p[,"C1D08"]
f2 = f[,"C1D08"]
names(e2) = names(p2) = names(f2) = rownames(e)

tmp = cbind(e2,p2,f2)
colnames(tmp) = c("log2 fold-change from baseline","p-value","False Discovery Rate")
write.csv(tmp,file="DE results at C1D08.csv")
write.csv(tmp,file="figure data - 5b - DE results at C1D08.csv")

# summary of FDR (cohort 1):
table(rowSums(f < 0.05) > 0)
table(rowSums(f[, c("C3D01", "C4D01", "FUW04", "FUW12")] < 0.05) > 0)
table(rowSums(f[, c("C3D01", "C4D01", "FUW04", "FUW12")] < 0.50) > 0)
table(rowSums(p[, c("C3D01", "C4D01", "FUW04", "FUW12")] < 0.05) > 0)

######
# model delta from baseline within cohort 2
######

# initialize results matrices: e for estimated log2 fold-change; p for p-values:
e.c2 = se.c2 = p.c2 = matrix(NA,ncol(norm),length(levels(annot$time))-1,
                    dimnames=list(colnames(norm),setdiff(levels(annot$time),"baseline")))
# at each timepoint, run the model:
for(time in setdiff(levels(annot$time),c("baseline","C1D01")))
{
  use = (annot$cohort=="2")&(is.element(annot$time,c("baseline",time)))
  for(gene in colnames(norm))
  {
    # paired t-test:
    mod = t.test(deltas[(annot$cohort=="2")&(is.element(annot$time,time)),gene])
    e.c2[gene,time] = mod$estimate  #mod[2,1]
    p.c2[gene,time] = mod$p.value   #mod[2,4]
  }
}
# get fdrs:
f.c2 = p.c2*NA
for(i in 1:ncol(p.c2)){f.c2[,i]=p.adjust(p.c2[,i], method = "BH")}
# summary of FDR (cohort 1):
table(rowSums(f.c2 < 0.05, na.rm = T) > 0)

#### write pvalues and estimates for supplemental material:
out.c1 = cbind(e, p, f)
colnames(out.c1) = paste0(rep(colnames(e), 3), " ", rep(c("log2fc", "p.value", "BH.FDR"), each = ncol(e)))

out.c2 = cbind(e.c2[, -1], p.c2[, -1], f.c2[, -1])
colnames(out.c2) = paste0(rep(colnames(e.c2[, -1]), 3), " ", rep(c("log2fc", "p.value", "BH.FDR"), each = ncol(e.c2)-1))

write.csv(out.c1, file = "differential expression results - cohort 1.csv")
write.csv(out.c2, file = "differential expression results - cohort 2.csv")
fig5aout = cbind(out.c1, out.c2)
colnames(fig5aout) = paste0(c(rep("cohort1", ncol(out.c1)), rep("cohort2", ncol(out.c2))), " ", colnames(fig5aout))
write.csv(fig5aout, file = "figure data - 5a - differential expression results.csv")
 
######
# for cohort 2 comparison, model delta from baseline at C1D08:
######
e2.c2 = p2.c2 = c()
for(gene in colnames(norm))
{
  # paired t-test:
  mod = t.test(deltas[(annot$cohort=="2")&(is.element(annot$time,"C1D08")),gene])
  e2.c2[gene] = mod$estimate  
  p2.c2[gene] = mod$p.value   
}
f2.c2 = p.adjust(p2.c2, method = "BH")
## and exceprt top genes for main body of paper:
top.up = rownames(tmp)[order(-log(p2) * (e2 > 0), decreasing = T)[1:10]]
top.dn = rownames(tmp)[order(-log(p2) * (e2 < 0), decreasing = T)[1:10]]
# calc mean and SE of change from baseline at C1D08:
tempc1 = deltas[(annot$cohort=="1")&(is.element(annot$time,"C1D08")), c(top.up, top.dn)]
tempc2 = deltas[(annot$cohort=="2")&(is.element(annot$time,"C1D08")), c(top.up, top.dn)]

countnotna = function(x) sum(!is.na(x))
genetable = cbind(
  colnames(tempc1),
  paste0(round(colMeans(tempc1, na.rm=T),2 ), " (",
         round(apply(tempc1, 2, sd, na.rm=T) / sqrt(apply(tempc1, 2, countnotna)), 2),
         ")"),
  paste0(round(colMeans(tempc2, na.rm=T),2 ), " (",
         round(apply(tempc2, 2, sd, na.rm=T) / sqrt(apply(tempc2, 2, countnotna)), 2),
         ")")
)
# version with standard errors:
colnames(genetable) = c("Gene", "Cohort 1 C1D08 log2 fold-change from baseline", "Cohort 2 C1D08 log2 fold-change from baseline")
write.csv(genetable, file = "top genes C1D08 results.csv", row.names = F)
# version with p-values:
ttestpval = function(x) {
  mod = t.test(x)
  return(mod$p.value)
}
genetable.p = cbind(
  colnames(tempc1),
  paste0(round(colMeans(tempc1, na.rm=T),2 ), " (p = ",
         signif(apply(tempc1, 2, ttestpval), 2),
         ")"),
  paste0(round(colMeans(tempc2, na.rm=T),2 ), " (p = ",
         signif(apply(tempc2, 2, ttestpval), 2),
         ")")
)
colnames(genetable.p) = c("Gene", 
                         "Cohort 1 C1D08 log2 fold-change from baseline (p-value)", 
                         "Cohort 2 C1D08 log2 fold-change from baseline (p-value)")
write.csv(genetable.p, file = "top genes C1D08 results with p-value.csv", row.names = F)
# volc plot:
#pdf("volcanos - cohort 1 timepoints.pdf",height=4*ncol(e),width=7)
#drawvolc(estmat=e,pvalmat=p,fdrmat=f,genesets=gannot,is.score = FALSE)
#dev.off()


pdf("C1D08 volcano plot only.pdf", width = 8, height = 5.7)
#par(mfrow = c(1,2))
#layout(matrix(1:2, 1), widths = c(4.8, 6))
par(mar = c(5.5,4.5,2,1))
par(xpd = FALSE)
# genes to highlight: FDR < 0.05
top = f2<0.05
show.text = p2 < 10^-4  #f2<0.05

plot(e2, -log10(p2), pch=16, cex = 0.6, 
#     col = c(fadecols("firebrick", 0.5), fadecols("darkblue", 0.5), "white", "white")[
#       (e2 > 0) + 2 * (e2 < 0) + 2 * show.text],
     col = c("grey60", "white")[1 + show.text],
     #xlab=paste0("  \n\nlog2 fold-change\nC1D08 vs. Pre-treatment"),
     xlab = "log2 fold-change",
     ylab = "-log10(p-value)", cex.lab = 1.7)
     #main="cohort 1: C1D08 vs. baseline",ylab="-log10(p-value)",cex.lab=1.5)
# text for top 50:
text(e2[show.text], -log10(p2[show.text]), colnames(norm)[show.text], cex=.6,
     col=c("darkblue","firebrick")[1+(e2[show.text]>0)])  #aquamarine4
# draw lines for FDR cutoffs:
fdr.cutoff = 0.05 #c(0.5,0.1,0.05)
abline(h=-log10(max(p2[f2<fdr.cutoff])),lty=2)
abline(h=-log10(max(p2[f2<0.5])),lty=3)
legend("bottomright",lty=2:3,legend = paste0("FDR = ",c(fdr.cutoff, 0.5)), bg = "white")
#legend("bottomright",lty=1:length(fdr.cutoffs),legend = paste0("FDR = ",fdr.cutoffs))
dev.off()

## now run gene set analysis:
write.csv(gannot, file = "figure data - 5c - gene sets.csv")
genesets = gannot
# summary stat for each gene set: mean -log10(pval)
genesetscore = c()
for(geneset in colnames(genesets))
{
  tempgenes = intersect(rownames(genesets)[genesets[,geneset]==1],colnames(norm))
  genesetscore[geneset] = mean(-log10(p2[tempgenes]))
}
# identify the top 5 gene sets:
#top5 = names(genesetscore)[order(genesetscore,decreasing = T)[1:30]]
## selected gene set scores based on biological interest:
#top5 = c("MHC2","APM","immunoproteasome","TLR","Cytokines","Myleoid.inflam","IFN.downstream","myleoid","lymphoid","NK.Cell.Functions","T.Cell.Functions")
top5 = c("MHC2","APM","immunoproteasome","TLR","Myleoid.inflam","IFN.downstream",
         "cell.trafficking","regulation.of.inflammation","DAMP.detection","myeloid.trafficking",
         "myeloid.signaling","T.cell.activation.inhibition.1","T.cell.Signaling","cytolytic.activity","NK.Cell.Functions")
#top5 = c(top5, setdiff(colnames(genesets.small), top5))
top5.altnames = top5; names(top5.altnames) = top5
top5.altnames["MHC2"] = "MHC2" # "Antigen processing by MHC2"
top5.altnames["APM"] = "MHC1" # "Antigen processing by MHC1"
top5.altnames["immunoproteasome"] = "Proteasome" # "Immunoproteasome"
top5.altnames["myleoid"] = "Myeloid" # "Myeloid compartment"
top5.altnames["lymphoid"] = "Lymphoid" # "Lymphoid compartment"
top5.altnames["TLR"] = "TLR" # "Toll Like Receptors"
#top5.altnames["IR.Innate"] = "MHC1" # "Innate immunity"
top5.altnames["Cytokines"] = "Cytokine" # "Cytokines"
top5.altnames["Myleoid.inflam"] = "Myeloid inflam" # "Myeloid inflammatory"
top5.altnames["IFN.downstream"] = "IFN" # "Interferon downstream"
top5.altnames["NK.Cell.Functions"] = "NK cell" # "NK cell functions"
top5.altnames["T.Cell.Functions"] = "T cell" # "T cell functions"
top5.altnames["cell.trafficking"] = "Trafficking" # 
top5.altnames["regulation.of.inflammation"] = "Inflammation" # 
top5.altnames["DAMP.detection"] = "DAMP detection" # 
top5.altnames["myeloid.trafficking"] = "Myeloid trafficking" # 
top5.altnames["myeloid.signaling"] = "Myeloid signalling" # 
top5.altnames["T.cell.activation.inhibition.1"] = "T cell inhibition" # 
top5.altnames["T.cell.Signaling"] = "T cell signalling" # 
top5.altnames["cytolytic.activity"] = "Cytolytic" # 
top5.altnames["NK.Cell.Functions"] = "NK cells" # 

# now plot the gene sets' -log10 pvals:
pdf("C1D08 volcano - gene set only.pdf", width = 8)
par(mar = c(12,4.5,2,1))
plot(c(0,0),col=0,xlim=c(0.5,length(top5)+0.5),ylim=range(-log10(p2)),xaxt="n",xlab="",
     ylab = "-log10(p-value)", cex.lab = 1.7)
bp = 1:length(top5)
#bp = barplot(genesetscore[top5],col=rgb(0,0,0,.1),ylim=ylim,xaxt="n",xlab="",ylab="-log10(p-value)",border = F,main="Gene set results")
axis(1,at = bp,labels = top5.altnames[top5],las=2, cex.axis = 1.6)  
abline(h=-log10(max(p2[f2<fdr.cutoff])),lty=2)
abline(h=-log10(max(p2[f2<0.5])),lty=3)

abline(v = bp,col="grey80")
for(i in 1:length(top5))
{
  geneset = top5[i]
  rect(bp[i]-.4,0,bp[i]+.4,genesetscore[geneset],col="grey70",border = F)
  tempgenes = intersect(rownames(genesets)[genesets[,geneset]==1],names(p2))
  text(rep(bp[i],length(tempgenes)),-log10(p2[tempgenes]),tempgenes,cex=0.6,
       col = c("darkblue","firebrick")[1+(e2[tempgenes]>0)])
  #col = c(rgb(0,0,1,0.7),rgb(1,0,0,0.7))[1+(e2[tempgenes]>0)])
}
dev.off()




#### heatmap of cohort 1 DE results:

# data to show: estimated fold-changes, but only where fdr< 0.5
mat = e*(p<0.05)
mat = mat[rowSums(abs(mat))>0, ]
clust = hclust(dist(mat))
pdf("DE heatmaps sideways.pdf", height = 3, width = 6)#, onefile = F)
#pheatmap(mat[order(mat[,"C1D08"]), ], cluster_cols = F, cluster_rows = F,
pheatmap(t(mat[as.character(clust$labels[clust$order]), ]), cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(101), 
         breaks = seq(-max(abs(mat)),max(abs(mat)),length.out = 100),
         show_colnames = F) #, main = "Mean log2 fold-changes from baseline, cohort 1")
dev.off()

# cohorts 1 & 2 side-by-side:
#clustorder = as.character(clust$labels[clust$order])
#mat2 = rbind((e*(p < 0.05))[clustorder, ], (e.c2 * (p.c2 < 0.05))[clustorder, ])
#mat2 = mat2[rowSums(abs(mat2))>0, ]
mat1 = e*(p<0.05)
mat2 = e.c2*(p.c2<0.05)
ever.sig = rowSums(abs(cbind(mat1, mat2)), na.rm=T)>0
mat1 = mat1[ever.sig, ]
mat2 = mat2[ever.sig, ]
clust = hclust(dist(mat1))
clust.order = as.character(clust$labels[clust$order])
clust2 = hclust(dist(mat2))
clust.order2 = as.character(clust2$labels[clust2$order])
# reformat mat2 to account for missing C1D01:
#mat2[, "C1D01"] = 0
colnames(mat2)[1] = ""

pdf("DE heatmaps sideways cohorts12 - cohort 1.pdf", height = 3, width = 3)#, onefile = F)
pheatmap(t(mat1[clust.order, ]), cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(101), 
         breaks = seq(-max(abs(mat1), na.rm = T),max(abs(mat1), na.rm = T),length.out = 100),
         show_colnames = F, legend = FALSE)
dev.off()
pdf("DE heatmaps sideways cohorts12 - cohort 2.pdf", height = 3, width = 3.2)#, onefile = F)
pheatmap(t(mat2[clust.order2, ]), cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(101), 
         breaks = seq(-max(abs(mat2), na.rm = T),max(abs(mat2), na.rm = T),length.out = 100),
         show_colnames = F, na_col = "white")
dev.off()



#pdf("DE heatmaps.pdf", height = 8, width = 4, onefile = F)
#pheatmap(mat[as.character(clust$labels[clust$order]), ], cluster_cols = F, cluster_rows = F,
#         color = colorRampPalette(c("blue","white","red"))(101), 
#         breaks = seq(-max(abs(mat)),max(abs(mat)),length.out = 100),
#         show_rownames = F) #, main = "Mean log2 fold-changes from baseline, cohort 1")
###gannot.show = c("MHC2","myleoid","lymphoid")
##pheatmap(gannot[clust$labels[clust$order], gannot.show], cluster_cols = F, cluster_rows = F,
##         show_rownames = F, col = c("white","darkblue"))
#dev.off()


####  trendplots:
# individual genes: CD74, HLA-DRA, HLA-DPB1, FLT3LG, FLT3, CD1D
# signatures: MHC2, TLRs


### calculate scores:
#scores = data.frame(rowMeans(norm[, rownames(gannot)[gannot[, "MHC2"]==1]]))
scores = data.frame(rowMeans(norm[, setdiff(rownames(gannot)[gannot[, "MHC2"]==1], "HLA-DQB1")]))
colnames(scores)[1] = "MHC2"
scores$APM = rowMeans(norm[, rownames(gannot)[gannot[, "APM"]==1]])
scores$TLR = rowMeans(norm[, rownames(gannot)[gannot[, "TLR"]==1]])


pdf("spaghetti plots.pdf", width = 7, height = 7)
# genes:
#par(mfrow=c(3,3))
layout(mat = t(matrix(1:9, 3)), heights = rep(c(3,3,4.5), each = 3))
par(mar = c(1, 2, 0.1, 0))
name = "HLA-DMA";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
#name = "IFI27";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "IL13RA1";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "CD74";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "FLT3LG";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "FLT3";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "SIGLEC1";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
# bottom row:
par(mar = c(6, 2, 0.1, 0))
name = "CD1C";trendplot(x=norm[,name],ylab="",xlabels=NULL,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "CD1D";trendplot(x=norm[,name],ylab="",xlabels=NULL,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))

# scores:
name = "MHC2";trendplot(x=scores[,name],ylab="",xlabels=NULL,
                        patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
legend("topright", legend = "mean of HLA-DRA,\n-DPB1, -DMA, -DPA1", cex = 1.5, box.col = rgb(0,0,0,0))
#name = "APM";trendplot(x=scores[,name],ylab=top5.altnames[name],patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
#name = "TLR";trendplot(x=scores[,name],ylab=top5.altnames[name],patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
dev.off()


pdf("spaghetti plots - tier1.pdf", width = 8.5, height = 2.35)
par(mfrow=c(1,3))
par(mar = c(1, 2, 0.3, 0))
name = "HLA-DMA";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "IL13RA1";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "CD74";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
dev.off()
pdf("spaghetti plots - tier2.pdf", width = 8.5, height = 2.35)
par(mfrow=c(1,3))
par(mar = c(1, 2, 0.3, 0))
name = "FLT3LG";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "FLT3";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "SIGLEC1";trendplot(x=norm[,name],ylab="",xlabels=rep("",12),patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
dev.off()
# bottom row:
pdf("spaghetti plots - tier3.pdf", width = 8.5, height = 3)
par(mfrow=c(1,3))
par(mar = c(6, 2, 0.3, 0))
name = "CD1C";trendplot(x=norm[,name],ylab="",xlabels=NULL,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "CD1D";trendplot(x=norm[,name],ylab="",xlabels=NULL,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol); legend("topright", legend = name, cex = 1.5, box.col = rgb(0,0,0,0))
name = "MHC2";trendplot(x=scores[,name],ylab="",xlabels=NULL,
                        patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
legend("topright", legend = "mean of HLA-DRA,\n-DPB1, -DMA, -DPA1", cex = 1.5, box.col = rgb(0,0,0,0))
#name = "APM";trendplot(x=scores[,name],ylab=top5.altnames[name],patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
#name = "TLR";trendplot(x=scores[,name],ylab=top5.altnames[name],patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
dev.off()


### heatmap of all data for cartoon:
pheatmap(pmax(pmin(as.matrix(scale(norm[1:50, 1:10])),3),-3), scale = "none",
         show_colnames = FALSE, show_rownames = FALSE,
         color = colorRampPalette(c("blue","white","red"))(101), 
         breaks = seq(-3, 3,length.out = 100), legend = FALSE,
         cluster_rows = F, cluster_cols = F)
         
