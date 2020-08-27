#### volcano plot source. also trendplots.

## trendplot code:
trendplot = function(x,ylab=NULL,xlabels=NULL,patient,groups,pointcols,show.ids=FALSE)
{
  pts = unique(patient)
  group = groups
  jgroup = jitter(as.numeric(groups),factor=0.5)
  xlim=NULL
  if (length(xlabels) == 0) {xlabels = levels(groups)}
  if(show.ids){xlim=range(as.numeric(as.factor(groups)))+c(-.08,.08)*diff(range(as.numeric(as.factor(groups))))}
  plot(jgroup,x,xaxt="n",xlim=xlim,ylab=ylab,xlab="",pch=16,col=pointcols,cex.lab=1.5)
  axis(1,at = as.numeric(1:length(levels(groups))),labels = xlabels,cex.axis = 1.3,las=2)
  for(pt in pts)
  {
    use = which(patient==pt)
    if(length(use)>0)
    {
      o = order(groups[use])
      lines(x=jgroup[use][o],y=x[use][o],col=pointcols[use[1]])
    }
    if(show.ids)
    {
      text(min(as.numeric(groups[use]))-0.5,x[use][o[1]],pt,cex=0.5,col="grey30")
      text(max(as.numeric(groups[use]))+0.5,x[use][o[length(o)]],pt,cex=0.5,col="grey30")
    }
  }
  
}

  
## volc plot code: estmat,pvalmat,fdrmat: matrices of log2fcs, pvals, etc. each row in a gene, each column is a contrast. (so a simple t-test model will have 1-column matrices)
# genesets is a 1-0 matrix of gene set membership. It might be required to align with the gene names.
# is.score is a logical vector saying which genes to highlight as a "score"
drawvolc = function(estmat,pvalmat,fdrmat,genesets,is.score=rep(FALSE,nrow(estmat)))
{
  if(length(colnames(estmat))==0){colnames(estmat)=colnames(pvalmat)=colnames(fdrmat)=""}
  ylim = range(-log10(pvalmat))
  xlim = range(estmat)
  par(mfrow = c(ncol(estmat),2))
  par(mar = c(12,4,2,1))
  for(j in 1:ncol(estmat))
  {
    level = colnames(estmat)[j]
    main=level
    #if(length(grouping.name.for.response)==0){main = ""}
    #if(length(grouping.name.for.response)>0){main = paste0(grouping.name.for.response," = ",level)}
    # id the top 50:
    top = unique(c(order(pvalmat[,j])[1:50],which(is.score)))
    # volc plot:
    plot(estmat[,j],-log10(pvalmat[,j]),pch=1,col=c("grey50","white")[1+is.element(1:nrow(pvalmat),top)],xlim=xlim,ylim=ylim,
         xlab=paste0("log2 fold-change"),
         main=main,ylab="-log10(p-value)",cex.lab=1.5)
    # text for top 50:
    text(estmat[top,j],-log10(pvalmat[top,j]),rownames(pvalmat)[top],cex=.7,col=c("aquamarine4","darkorange2")[1+is.element(top,which(is.score))])
    # draw lines for FDR cutoffs:
    fdr.cutoffs = c(0.5,0.1,0.05)
    for(i in 1:length(fdr.cutoffs))
    {
      suppressWarnings(abline(h=-log10(suppressWarnings((max(pvalmat[fdrmat[,j]<fdr.cutoffs[i],j])))),lty=i))
    }
    legend("bottomleft",lty=1:length(fdr.cutoffs),legend = paste0("FDR = ",fdr.cutoffs))
    
    ## now run gene set analysis:
    # summary stat for each gene set: mean -log10(pval)
    genesetscore = c()
    for(geneset in colnames(genesets))
    {
      tempgenes = intersect(rownames(genesets)[genesets[,geneset]==1],rownames(pvalmat))
      genesetscore[geneset] = mean(-log10(pvalmat[tempgenes,j]))
    }
    # identify the top 5 gene sets:
    top5 = names(genesetscore)[order(genesetscore,decreasing = T)[1:10]]
    
    # now plot the gene sets' -log10 pvals:
    plot(c(0,0),col=0,xlim=c(1,length(top5)),ylim=ylim,xaxt="n",xlab="",ylab="-log10(p-value)")
    bp = 1:length(top5)
    #bp = barplot(genesetscore[top5],col=rgb(0,0,0,.1),ylim=ylim,xaxt="n",xlab="",ylab="-log10(p-value)",border = F,main="Gene set results")
    axis(1,at = bp,labels = top5,las=2)  
    #abline(v = 1:length(top5),col="grey80")
    abline(v = bp,col="grey80")
    for(i in 1:length(top5))
    {
      geneset = top5[i]
      rect(bp[i]-.4,ylim[1],bp[i]+.4,genesetscore[geneset],col=rgb(0,0,0,.5),border = F)
      tempgenes = intersect(rownames(genesets)[genesets[,geneset]==1],rownames(pvalmat))
      text(rep(bp[i],length(tempgenes)),-log10(pvalmat[tempgenes,j]),tempgenes,cex=0.5,
           col = c(rgb(0,0,1,0.7),rgb(1,0,0,0.7))[1+(estmat[tempgenes,j]>0)])
    }
  }
}

