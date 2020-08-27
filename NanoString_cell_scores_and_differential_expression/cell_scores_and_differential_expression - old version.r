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


### load the module-rene relationships, and format as a gene annotation matrix:
#modules = read.csv("CITN_nano_mod_Gen2_NoHK[3].csv", row.names = 1, stringsAsFactors = F)
#gannot = matrix(0, ncol(norm), length(unique(modules$Function)), dimnames = list(colnames(norm), unique(modules$Function)))
## fill it in:
#for (m in unique(modules$Function)){
#  modgenes = modules$Symbol[modules$Function == m]
#  gannot[is.element(rownames(gannot), modgenes), m] = 1
#}
## remove modules with too few genes to be useful:
#gannot = gannot[, colSums(gannot) >= 10]

##### reformat the gene annotation matrix:
gannot = gannot[match(colnames(norm), rownames(gannot)), ]
gannot = replace(gannot, is.na(gannot), 0)
rownames(gannot) = colnames(norm)

# load additional gene sets:
genesets = read.csv("gene sets.csv", stringsAsFactors = F)
gannot2 = matrix(0, ncol(norm), length(unique(genesets$sig)), dimnames = list(colnames(norm), unique(genesets$sig)))
# fill it in:
for (m in unique(genesets$sig)){
  modgenes = genesets$gene[genesets$sig == m]
  gannot2[is.element(rownames(gannot2), modgenes), m] = 1
}
# remove modules with too few genes to be useful:
gannot2 = gannot2[, colSums(gannot2) >= 3]

gannot = cbind(gannot, gannot2)


######################################
#### calculate cell type signatures: #### 
######################################
cellmarkers = read.csv("data/cell type markers.csv",row.names = 1)
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

#cells = c("B-cells","CD8 T cells","T-cells","Cytotoxic cells","DC","Macrophages","NK cells","Th1 cells","Treg")
#cells = c("B-cells","CD8 T cells","T-cells","Cytotoxic cells","Macrophages","NK cells","Th1 cells","Treg")
#cells = colnames(cellscores)
#cellcols = c("darkblue","orange","gold","chartreuse3","purple","blue","forestgreen","cornflowerblue","red","firebrick","grey50","lightblue","pink","green")[1:length(cells)]
cells = setdiff(colnames(cellscores), "CD45")
cellcols = c("darkblue","orange","chartreuse3","purple","blue","forestgreen","cornflowerblue","red","firebrick","grey50","lightblue","pink","green")[1:length(cells)]
names(cellcols)=cells

svg("manuscript - cells over time v1.svg",width=10)
nudge = 2
#layout(matrix(1:2, nrow = 1), width = c(4,4))
par(mar = c(5.5,4.5,2,1))
par(mfrow = c(1,2))
plot(0,0,xlim=c(1,sum(inc1)+1.5) + c(0,nudge),ylim=range(meandeltas[,cells]), xaxt="n",xlab="",
     ylab="Mean log2 fold-change from baseline",main="cohort1", cex.lab = 1.5)
abline(h=0,col = "grey20", lty = 3, lwd = 2)
axis(1,1:sum(inc1),timepoints,las=2, cex.axis = 1.5)
for(cell in cells)
{
  lines(1:sum(inc1),meandeltas[inc1,cell],col=cellcols[cell],lwd=2)
  text(sum(inc1)+0.5,meandeltas[max(which(inc1)),cell],cell,col=cellcols[cell],cex=1.25)
  # points for changes sig at p=0.05:
  is.sig = (meandeltas[(inc1),cell] - 2*sedeltas[(inc1),cell] > 0)|(meandeltas[(inc1),cell] + 2*sedeltas[(inc1),cell] < 0)
  points((1:sum(inc1))[is.sig],meandeltas[which(inc1)[is.sig],cell],col=cellcols[cell],pch=16,cex=1.5)
}

names(cellcols)=cells
plot(0,0,xlim=c(1,sum(inc2)+1.5) + c(0,nudge),ylim=range(meandeltas[,cells]),xaxt="n",xlab="",
     ylab="",main="cohort2", cex.lab = 1.5)
abline(h=0,col = "grey20", lty = 3, lwd = 2)
axis(1,1:sum(inc2),setdiff(timepoints, "C1D01"),las=2, cex.axis = 1.5)
for(cell in cells)
{
  lines(1:sum(inc2),meandeltas[inc2,cell],col=cellcols[cell],lwd=2)
  text(sum(inc2)+0.5,meandeltas[max(which(inc2)),cell],cell,col=cellcols[cell],cex=1.25)
  # points for changes sig at p=0.05:
  is.sig = (meandeltas[(inc2),cell] - 2*sedeltas[(inc2),cell] > 0)|(meandeltas[(inc2),cell] + 2*sedeltas[(inc2),cell] < 0)
  points((1:sum(inc2))[is.sig],meandeltas[which(inc2)[is.sig],cell],col=cellcols[cell],pch=16,cex=1.5)
}

dev.off()

svg("manuscript - cells over time v2.svg", width=9, height = 6)
par(mar = c(5.6,4.5,2,1))
layout(matrix(1:3, 1, 3), width = c(4,4,3))
plot(0,0,xlim=c(1,sum(inc1)),ylim=range(meandeltas[,cells]), xaxt="n",xlab="",
     ylab="Mean log2 fold-change from baseline",main="cohort1", cex.lab = 1.8)
abline(h=0,col = "grey20", lty = 3, lwd = 2)
axis(1,1:sum(inc1),timepoints,las=2, cex.axis = 1.5)
for(cell in cells)
{
  lines(1:sum(inc1),meandeltas[inc1,cell],col=cellcols[cell],lwd=2)
  #text(sum(inc1)+0.5,meandeltas[max(which(inc1)),cell],cell,col=cellcols[cell],cex=1.25)
  # points for changes sig at p=0.05:
  is.sig = (meandeltas[(inc1),cell] - 2*sedeltas[(inc1),cell] > 0)|(meandeltas[(inc1),cell] + 2*sedeltas[(inc1),cell] < 0)
  points((1:sum(inc1))[is.sig],meandeltas[which(inc1)[is.sig],cell],col=cellcols[cell],pch=16,cex=1.5)
}

names(cellcols)=cells
plot(0,0,xlim=c(1,sum(inc2)),ylim=range(meandeltas[,cells]),xaxt="n",xlab="",
     ylab="",main="cohort2", cex.lab = 1.5)
abline(h=0,col = "grey20", lty = 3, lwd = 2)
axis(1,1:sum(inc2),setdiff(timepoints, "C1D01"),las=2, cex.axis = 1.5)
for(cell in cells)
{
  lines(1:sum(inc2),meandeltas[inc2,cell],col=cellcols[cell],lwd=2)
  #text(sum(inc2)+0.5,meandeltas[max(which(inc2)),cell],cell,col=cellcols[cell],cex=1.25)
  # points for changes sig at p=0.05:
  is.sig = (meandeltas[(inc2),cell] - 2*sedeltas[(inc2),cell] > 0)|(meandeltas[(inc2),cell] + 2*sedeltas[(inc2),cell] < 0)
  points((1:sum(inc2))[is.sig],meandeltas[which(inc2)[is.sig],cell],col=cellcols[cell],pch=16,cex=1.5)
}
par(mar = c(5.6,0,0,0))
frame()
o.legend = order(meandeltas[3,names(cellcols)], decreasing = T)
legend("center", col = c("black", "black", NA, cellcols[o.legend]), 
       legend = c("p > 0.05", "p < 0.05", "", names(cellcols)[o.legend]), 
       lty = 1, pch = c(NA, 16, NA, rep(NA, length(o.legend))), lwd = 3, bty = "n", cex = 1.5)
dev.off()

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
#svg("volcanos - cohort 1 timepoints.svg",height=4*ncol(e),width=7)
#drawvolc(estmat=e,pvalmat=p,fdrmat=f,genesets=gannot,is.score = FALSE)
#dev.off()


svg("C1D08 volcano plot.svg", width = 10)
par(mfrow = c(1,2))
par(mar = c(14,4.5,2,0))
# genes to highlight: FDR < 0.05
top = f2<0.05
show.text = f2<0.05

plot(e2, -log10(p2),pch=16, cex = 0.5,col=c(rgb(0,0,0,.5),"white")[1 + show.text],
     xlab=paste0("  \n\nlog2 fold-change\nC1D08 vs. Pre-treatment"),
     ylab = "-log10(p-value)",cex.lab=1.2)
     #main="cohort 1: C1D08 vs. baseline",ylab="-log10(p-value)",cex.lab=1.5)
# text for top 50:
text(e2[show.text],-log10(p2[show.text]),colnames(norm)[show.text],cex=.5,col=c("darkblue","firebrick")[1+(e2[show.text]>0)])  #aquamarine4
# draw lines for FDR cutoffs:
fdr.cutoff = 0.05 #c(0.5,0.1,0.05)
abline(h=-log10(max(p2[f2<fdr.cutoff])),lty=2)
abline(h=-log10(max(p2[f2<0.5])),lty=3)
legend("bottomright",lty=2:3,legend = paste0("FDR = ",c(fdr.cutoff, 0.5)))
#legend("bottomright",lty=1:length(fdr.cutoffs),legend = paste0("FDR = ",fdr.cutoffs))


## now run gene set analysis:
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
top5 = c("MHC2","APM","immunoproteasome","TLR","Cytokines","Myleoid.inflam","IFN.downstream","myleoid","lymphoid","NK.Cell.Functions","T.Cell.Functions")
top5.altnames = top5; names(top5.altnames) = top5
top5.altnames["MHC2"] = "Antigen processing by MHC2"
top5.altnames["APM"] = "Antigen processing by MHC1"
top5.altnames["immunoproteasome"] = "Immunoproteasome"
top5.altnames["myleoid"] = "Myeloid compartment"
top5.altnames["lymphoid"] = "Lymphoid compartment"
top5.altnames["TLR"] = "Toll Like Receptors"
top5.altnames["IR.Innate"] = "Innate immunity"
top5.altnames["Cytokines"] = "Cytokines"
top5.altnames["Myleoid.inflam"] = "Myeloid inflammatory"
top5.altnames["IFN.downstream"] = "Interferon downstream"
top5.altnames["NK.Cell.Functions"] = "NK cell functions"
top5.altnames["T.Cell.Functions"] = "T cell functions"


# now plot the gene sets' -log10 pvals:
par(mar = c(14,1,2,1))
plot(c(0,0),col=0,xlim=c(0.5,length(top5)+0.5),ylim=range(-log10(p2)),xaxt="n",xlab="", yaxt="n", ylab = "")#,ylab="-log10(p-value)", cex.lab = 1.2)
bp = 1:length(top5)
#bp = barplot(genesetscore[top5],col=rgb(0,0,0,.1),ylim=ylim,xaxt="n",xlab="",ylab="-log10(p-value)",border = F,main="Gene set results")
axis(1,at = bp,labels = top5.altnames[top5],las=2)  
abline(h=-log10(max(p2[f2<fdr.cutoff])),lty=2)
abline(h=-log10(max(p2[f2<0.5])),lty=3)

abline(v = bp,col="grey80")
for(i in 1:length(top5))
{
  geneset = top5[i]
  rect(bp[i]-.4,0,bp[i]+.4,genesetscore[geneset],col=rgb(0,0,0,.5),border = F)
  tempgenes = intersect(rownames(genesets)[genesets[,geneset]==1],names(p2))
  text(rep(bp[i],length(tempgenes)),-log10(p2[tempgenes]),tempgenes,cex=0.5,
       col = c("darkblue","firebrick")[1+(e2[tempgenes]>0)])
  #col = c(rgb(0,0,1,0.7),rgb(1,0,0,0.7))[1+(e2[tempgenes]>0)])
}
dev.off()


#### heatmap of cohort 1 DE results:

# data to show: estimated fold-changes, but only where fdr< 0.5
mat = e*(p<0.05)
#mat = cbind(0, mat)
#colnames(mat)[1] = "Pre-tx"
mat = mat[rowSums(abs(mat))>0, ]
clust = hclust(dist(mat))
svg("DE heatmaps sideways.svg", height = 3, width = 6, onefile = F)
#pheatmap(mat[order(mat[,"C1D08"]), ], cluster_cols = F, cluster_rows = F,
pheatmap(t(mat[as.character(clust$labels[clust$order]), ]), cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(101), 
         breaks = seq(-max(abs(mat)),max(abs(mat)),length.out = 100),
         show_colnames = F, main = "Mean log2 fold-changes from baseline, cohort 1")
dev.off()
svg("DE heatmaps.svg", height = 8, width = 4, onefile = F)
pheatmap(mat[as.character(clust$labels[clust$order]), ], cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(101), 
         breaks = seq(-max(abs(mat)),max(abs(mat)),length.out = 100),
         show_rownames = F) #, main = "Mean log2 fold-changes from baseline, cohort 1")
#gannot.show = c("MHC2","myleoid","lymphoid")
#pheatmap(gannot[clust$labels[clust$order], gannot.show], cluster_cols = F, cluster_rows = F,
#         show_rownames = F, col = c("white","darkblue"))
dev.off()


####  trendplots:
# individual genes: CD74, HLA-DRA, HLA-DPB1, FLT3LG, FLT3, CD1D
# signatures: MHC2, TLRs


### calculate scores:
#scores = data.frame(rowMeans(norm[, rownames(gannot)[gannot[, "MHC2"]==1]]))
scores = data.frame(rowMeans(norm[, setdiff(rownames(gannot)[gannot[, "MHC2"]==1], "HLA-DQB1")]))
colnames(scores)[1] = "MHC2"
scores$APM = rowMeans(norm[, rownames(gannot)[gannot[, "APM"]==1]])
scores$TLR = rowMeans(norm[, rownames(gannot)[gannot[, "TLR"]==1]])


svg("spaghetti plots.svg", width = 10, height = 10)
# genes:
par(mfrow=c(3,3))
name = "HLA-DMA";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
#name = "IFI27";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "IL13RA1";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "CD74";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "FLT3LG";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "FLT3";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "SIGLEC1";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "CD1C";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
name = "CD1D";trendplot(x=norm[,name],ylab=name,patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)


# scores:

name = "MHC2";trendplot(x=scores[,name],ylab="mean of HLA-DRA, -DPB1, -DMA, -DPA1",
                        patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
#name = "APM";trendplot(x=scores[,name],ylab=top5.altnames[name],patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
#name = "TLR";trendplot(x=scores[,name],ylab=top5.altnames[name],patient=annot$patient,groups=annot$time,pointcols=annot$cohortcol)
dev.off()

