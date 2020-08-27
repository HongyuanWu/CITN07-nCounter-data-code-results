#### questions for CITN:
# - C1D-7 timepoint: only in cohort 1. Can I consider C1D01 the "pre-tx" timepoint and ignore C1D-7 in my models?
# - 

rm(list=ls())
source("RCC_gather_n_collect.r")

library(NormqPCR)
library(RColorBrewer)

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

datasets = collate.lanes.2.dataset(allRCC.path = paste0("input_data/CITN_lysate RCC_258_Jan2018"),
                                   sup.path = paste0("input_data/CITN_lysate RCC_258_Jan2018","//"))$datasets
dat = datasets[[1]]
# there's one sample run in dup. Delete it:
dup.rcc = "20171228_CITN 22_107-27-013 FUW12_04"
remove = dat$sample.annot$sample.Name==dup.rcc
dat$count = dat$count[!remove,]
dat$sample.annot = dat$sample.annot[!remove,]

for(i in c("count","sample.annot")){rownames(dat[[i]])=dat$sample.annot$Sample.ID}
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
raw = dat$count[sharedids,]
#raw = log2(pmax(dat$count[sharedids,],1))


# split out controls:
neg = raw[,grepl("NEG_",colnames(raw))]
pos = raw[,grepl("POS_",colnames(raw))]
raw = raw[,!grepl("NEG_",colnames(raw))&!grepl("POS_",colnames(raw))]

######
## Signal QC: genes
######

# calculate maximum count of each gene, and save it as a signal threshold:
max.counts = apply(raw, 2, max)
#hist(log2(pmax(max.counts,1)))

# drop genes in background: all with max < 32 counts
#minrawthresh = 25
#in.bg = (apply(raw,2,max)<log2(minrawthresh))
#plot(apply(raw,2,mean),apply(raw,2,sd))
#plot(apply(raw,2,mean),apply(raw,2,max),col=1+in.bg)
#abline(h=log2(minrawthresh))

# remove low signal genes:
#raw = raw[,!in.bg]

######
# get HKs and normalize:
######
is.hk = intersect(colnames(raw),rownames(dat$gannot)[dat$gannot$CodeClass=="Housekeeping"])
# run genorm to get good subset:
temp = selectHKs(log2(pmax(raw[,is.hk], 1)), log=TRUE, minNrHKs=10, Symbols=is.hk)
selectedhks = temp$ranking[names(temp$ranking)<=12]

# calc hk geomean:
hk.factor = apply(log2(pmax(raw[,selectedhks], 1)),1,mean) 

# remove low-signal lanes:
hist(hk.factor,breaks=40,col="grey")
abline(v=4)
lowsignal = hk.factor<4
raw = raw[!lowsignal,]
neg = neg[!lowsignal,]
pos = pos[!lowsignal,]
annot = annot[!lowsignal,]
hk.factor = hk.factor[!lowsignal]
# normalize:
norm = sweep(raw, 1, 2^hk.factor, "/") * 2^mean(hk.factor)
# remove HKs from analysis:
norm = norm[,setdiff(colnames(norm),is.hk)]
# log-transformed normalized data:
log2norm = log2(pmax(norm, 1))

######
## look at structure of annot:
######
# cohort as a string:
annot$cohort=factor(as.character(annot$cohort),levels=c("2","1"))
# extract patient ID:
annot$patient = substr(annot$sample.ID,1,10)
# timepoint variables:
table(annot$cohort,annot$timepoint)
table(annot$course,annot$timepoint)
# there's not an single time variable, but rather the combo of timepoint and course gives time. 
annot[annot$patient=="107-11-015",]
# and it looks like the timepoint and course columns have errors anyway - or at least, they differ from the sample.id column

# extract time from sample ID:
annot$time = substr(annot$sample.ID,12,16)
unique(annot$time)
annot$time[annot$time=="CID08"]="C1D08"
annot$time[annot$time=="CID22"]="C1D22"
annot$time[annot$time=="CID01"]="C1D01"
annot$time[annot$time=="CID-7"]="C1D-7"
# recode time variable to accommodate weird baselines:
annot$time[annot$time=="C1D-7"] = "baseline"
annot$time[(annot$time=="C1D01")&(annot$cohort==2)] = "baseline"

annot$time = factor(annot$time,
                    levels = c("baseline","C1D01","C1D08","C1D15","C1D22",
                               "C2D01","C2D08","C2D15","C3D01","C4D01",
                               "FUW04","FUW12"))
table(annot$patient,annot$time)
table(annot$time,annot$cohort)


### assign colors:
timecols = colorspace::diverge_hcl(length(unique(annot$time)))
names(timecols)=levels(annot$time)
annot$timecol = timecols[annot$time]

cohortcols = c("cornflowerblue","orange")
names(cohortcols) = c("2","1")
annot$cohortcol = cohortcols[annot$cohort]
## -> so create a time variable:
#annot$time = factor(paste0(annot$course,annot$timepoint),
#                       levels = c( "C1d1","C1d8","C1d15","C1d22",
#                                   "C2d1","C2d8","C2d15","C2d22",
#                                   "C3d1","C4d1","FUWd4","FUWd12"))


######
# load gene set matrix and align to norm:
######
gannot = as.matrix(read.csv("input_data/PCI gene set matrix.csv",row.names = 1,stringsAsFactors = F))
gannot = gannot[match(colnames(norm),rownames(gannot)),]
gannot = replace(gannot,is.na(gannot),0)

save(raw,norm,log2norm,annot,gannot,hk.factor,max.counts,file="data_post_QC_and_normalization/CITN07 data ready for analysis October2018.RData")
save(raw,norm,log2norm,annot,gannot,hk.factor,max.counts,file="../NanoString_cell_scores_and_differential_expression/data/CITN07 data ready for analysis October2018.RData")
write.csv(log2norm, file = "data_post_QC_and_normalization/CITN07 log2-transformed normalized expression data ready for analysis October2018.csv")
write.csv(norm, file = "data_post_QC_and_normalization/CITN07 normalized expression data ready for analysis October2018.csv")
write.csv(annot, file = "data_post_QC_and_normalization/CITN07 sample annotation data ready for analysis October2018.csv")
write.csv(gannot, file = "data_post_QC_and_normalization/CITN07 gene annotation data ready for analysis October2018.csv")
