##### This script loads in the raw nanostring data, removes low-signal genes and samples,
#     normalizes the gene expression data, and creates well-formatted sample annotations.

source("RCC_gather_n_collect.r")

library(NormqPCR)
library(RColorBrewer)
library(colorspace)

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


#####
# load RCCs
#####

datasets = collate.lanes.2.dataset(allRCC.path = paste0("input_data/CITN_lysate RCC_258_Jan2018/"),
                                   sup.path = paste0("input_data/CITN_lysate RCC_258_Jan2018","//"))$datasets
dat = datasets[[1]]

###### look out for duplicates: #####
is.dup = names(table(dat$sample.annot$Sample.ID))[table(dat$sample.annot$Sample.ID)>1]
rownames(dat$sample.annot[is.element(dat$sample.annot$Sample.ID, is.dup), ])
# any difference in signal?
barplot(rowMeans(log2(pmax(dat$count[is.element(dat$sample.annot$Sample.ID, is.dup), ],1))), las = 2)

# there's one sample run in dup. Delete one of the duplicates:
dup.rcc = "20171228_CITN 22_107-27-013 FUW12_04"
remove = dat$sample.annot$sample.Name==dup.rcc
dat$count = dat$count[!remove,]
dat$sample.annot = dat$sample.annot[!remove,]

for(i in c("count","sample.annot")){rownames(dat[[i]])=dat$sample.annot$Sample.ID}

# extract raw data from the results object:
raw = log2(pmax(dat$count,1))

#####
# load annot
#####
annot = read.csv("input_data/CITN07 sample cohort info.csv",stringsAsFactors = F,row.names=1)
# remove the dup rcc from here as well:
head(annot)
remove = rownames(annot)==paste0(dup.rcc,".RCC")
annot = annot[!remove,]
rownames(annot) = annot$sample.ID

#####
# align annot to rccs, evaluate overlap
#####
sharedids = intersect(rownames(dat$count),rownames(annot))
setdiff(rownames(dat$count),rownames(annot))
setdiff(rownames(annot),rownames(dat$count))
annot = annot[sharedids,]
raw = log2(pmax(dat$count[sharedids,],1))
#rm(dat)
#### split out controls from gene expression data:
neg = raw[,grepl("NEG_",colnames(raw))]
pos = raw[,grepl("POS_",colnames(raw))]
raw = raw[,!grepl("NEG_",colnames(raw))&!grepl("POS_",colnames(raw))]

######
## Signal QC: genes
######

# drop genes in background: all with max below a threshold number of counts:
minrawthresh = 25  
in.bg = (apply(raw,2,max)<log2(minrawthresh))
plot(apply(raw,2,mean),apply(raw,2,sd))
plot(apply(raw,2,mean),apply(raw,2,max),col=1+in.bg)
abline(h=log2(minrawthresh))

# remove low signal genes:
raw = raw[,!in.bg]

######
# get HKs and normalize:
######
is.hk = intersect(colnames(raw),rownames(dat$gannot)[dat$gannot$CodeClass=="Housekeeping"])
# run genorm to get good subset:
temp = selectHKs(raw[,is.hk],log=TRUE,minNrHKs=10,Symbols=is.hk)
selectedhks = temp$ranking[names(temp$ranking)<=14]

# calc hk geomean:
hk.factor = apply(raw[,selectedhks],1,mean) 

# remove low-signal lanes:
hist(hk.factor,breaks=40,col="grey")
abline(v=5)
lowsignal = hk.factor<5
raw = raw[!lowsignal,]
neg = neg[!lowsignal,]
pos = pos[!lowsignal,]
annot = annot[!lowsignal,]
hk.factor = hk.factor[!lowsignal]
# normalize:
norm = sweep(raw,1,hk.factor,"-")+mean(hk.factor)
# remove HKs from analysis:
norm = norm[,setdiff(colnames(norm),is.hk)]





### assign colors:
timecols = colorspace::diverge_hcl(length(unique(annot$time)))
names(timecols)=levels(annot$time)
annot$timecol = timecols[annot$time]

cohortcols = c("cornflowerblue","orange")
names(cohortcols) = c("2","1")
annot$cohortcol = cohortcols[annot$cohort]

######
# load gene set matrix and align to norm:
######
gannot = as.matrix(read.csv("input_data/PCI gene set matrix.csv",row.names = 1,stringsAsFactors = F))
gannot = gannot[match(colnames(norm),rownames(gannot)),]
gannot = replace(gannot,is.na(gannot),0)

save(norm,annot,gannot,file="data_post_QC_and_normalization/CITN07 normalized data and parsed sample annotation.RData")
save(norm,annot,gannot,file="../NanoString_cell_scores_and_differential_expressoin/data/CITN07 normalized data and parsed sample annotation.RData")
write.csv(norm, file = "data_post_QC_and_normalization/CITN07 normalized expression data ready for analysis October2018.csv")
write.csv(annot, file = "data_post_QC_and_normalization/CITN07 sample annotation data ready for analysis October2018.csv")
write.csv(gannot, file = "data_post_QC_and_normalization/CITN07 gene annotation data ready for analysis October2018.csv")
