#' This function grab and copy RCC files from a bunch of folders into a single folder
#' @param sup.path <- xxx # top folder containing RCC folders Example: "S:/AMashadi-Hossein/Sprint temp modeling RCC/"
#' @param allRCC.path <- yyy # Folder where the RCCs to be dumped in Example: "S:/AMashadi-Hossein/Take2/2016_10_04"
#' @return a data frame showing which files successfully got copied

Gather.RCC <- function(sup.path, allRCC.path, overwrite = T,recursive = T){
  if(!dir.exists(allRCC.path))
    dir.create(allRCC.path)
  
  getwd()->path0
  setwd(sup.path)
  folders <- dir()
  
  #This gets all the RCC files
  fls <- list.files(path=folders,full.names = TRUE,pattern = "\\.RCC$",recursive = T)
  
  #This copies all the RCC files
  all.passed<- file.copy(from = fls, to = allRCC.path, overwrite = overwrite, recursive = recursive) 
  if(!all(all.passed)){
    cat("some files failed to be copied")
  }else
    cat(sum(all.passed),"files copied")
  
  setwd(path0)
  return(data.frame(fls=fls,copied = all.passed))
  
  # Below we can get a mapping of RCC to folder name where they were stored in case sample annotations is contained there 
  #--------
  # strsplit(fls,"/") -> tmp
  # names(tmp) <- lapply(tmp, function(x)x[[2]])
  # tmp <- lapply(tmp, function(x) x[[1]])
  # tmp <- data.frame(file.name = t(data.frame(tmp)))
  # write.csv(x = tmp, file = "S:/AMashadi-Hossein/Take2/2016_10_04/RCC_to_folder_map.csv")
}


#' This funciton reads and combines RCC files into a single file
#' @param allRCC.path is the directory
#' @param fl.name is the file name

extract.lane.from.RCC <- function(allRCC.path,fl.name){
  tmp <- readLines(paste(allRCC.path,"/",fl.name,sep=""),warn = F)
  
  # smp attr
  #---------
  smp.attr.indices <- c((grep("<Sample_Attributes>",tmp) + 1) : (grep("</Sample_Attributes>",tmp)-1))
  smp.attr <- unlist(lapply(tmp[smp.attr.indices],function(x){spl<- unlist(strsplit(x,split = ","));
  y<-ifelse(length(spl)==2,spl[2],"");
  names(y)<-spl[1];
  return(y)}))
  names(smp.attr)[names(smp.attr)=="ID"] <- "Sample ID"
  # lane attr
  #---------
  lane.attr.indices <- c((grep("<Lane_Attributes>",tmp) + 1) : (grep("</Lane_Attributes>",tmp)-1))
  lane.attr <- unlist(lapply(tmp[lane.attr.indices],function(x){spl<- unlist(strsplit(x,split = ","));
  y<-ifelse(length(spl)==2,spl[2],"");
  names(y)<-spl[1];
  return(y)}))
  names(lane.attr)[names(lane.attr)=="ID"] <- "Lane Number"
  
  # Counts
  #--------
  count.attr.indices <- c((grep("<Code_Summary>",tmp) + 2) : (grep("</Code_Summary>",tmp)-1))
  count.attr <- data.frame(t(data.frame(lapply(tmp[count.attr.indices],function(x){spl<- unlist(strsplit(x,split = ","))}),stringsAsFactors = F)),stringsAsFactors = F)
  colnames(count.attr) <- unlist(strsplit(x =tmp[grep("<Code_Summary>",tmp) + 1] ,split = ","))
  rownames(count.attr) <- NULL
  count.attr$Count <- as.integer(as.character(count.attr$Count))
  
  sample.annot <-  data.frame(t(data.frame(c(sample.Name = gsub(".RCC","",fl.name),smp.attr,lane.attr))),stringsAsFactors = F)
  rownames(sample.annot) <- NULL
  
  return(list(sample.annot =sample.annot, count = count.attr))
}


collate.lanes.2.dataset <- function(allRCC.path,sup.path = NULL){
  
  
  # require(plyr)
  
  # Get all RCC files
  #------------------
  all.RCC.files <- dir(allRCC.path)
  all.RCC.files <- all.RCC.files[grep(".RCC",all.RCC.files)]
  if(length(all.RCC.files)==0){
    cat("No file is exist in all. will try to copy from sup.path\n")
    if(is.null(sup.path))
      stop("sup.path is null and all.RCC.files is empty. Provide sup.path or populate all.RCC.files")
    cat("Copying the files from sup.path into all.RCC.files")
    fls.copied <- Gather.RCC(sup.path = sup.path,allRCC.path = allRCC.path)
    all.RCC.files <- dir(allRCC.path)
    all.RCC.files <- all.RCC.files[grep(".RCC",all.RCC.files)]
  }
  
  
  # Initiate storage objects
  datasets <- list()
  current.genenames <- NULL
  fl.num <- 0
  
  for(fl in all.RCC.files){
    fl.num <- fl.num + 1
    tmp <- suppressWarnings(try(extract.lane.from.RCC(allRCC.path = allRCC.path,fl.name = fl),silent = T))

    
    if(class(tmp)=="try-error"){
      cat("File", fl, " Failed\n")
      next
    }
    tmpcounts <- tmp$count$Count
    names(tmpcounts) <- tmp$count$Name
    
    
    gene.names.match <- length(tmp$count$Name) == length(current.genenames)
    if(gene.names.match)
      gene.names.match <- all(tmp$count$Name == current.genenames)
    
    if(gene.names.match){

      datasets[[this.set]]$count <-  rbind(datasets[[this.set]]$count,tmpcounts)
      rownames(datasets[[this.set]]$count)[nrow(datasets[[this.set]]$count)] <- make.names(ifelse(length(tmp$sample.annot$sample.Name)>0,as.character(tmp$sample.annot$sample.Name),nrow(datasets[[this.set]]$count)))
      datasets[[this.set]]$sample.annot <- plyr:::rbind.fill(datasets[[this.set]]$sample.annot,tmp$sample.annot)
      rownames(datasets[[this.set]]$sample.annot) <- make.names(datasets[[this.set]]$sample.annot$sample.Name)
      
      cat("dataset= ",this.set,", file= ",fl.num,"\n")
      
      
    }else{
      
      this.set <- names(datasets)[unlist(lapply(datasets,FUN = function(x) isTRUE(all.equal(x$gannot$Name,tmp$count$Name))))]
      
      if(length(this.set)==0){
        dataset.count <- length(datasets)+1
        this.set <- paste0("ds",dataset.count)
        datasets[[this.set]]$gannot <- tmp$count[c("CodeClass", "Name","Accession")]
        rownames(datasets[[this.set]]$gannot) <- make.names(datasets[[this.set]]$gannot$Name)
      }
      current.genenames <- tmp$count$Name
      datasets[[this.set]]$count <-  rbind(datasets[[this.set]]$count,tmpcounts)
      rownames(datasets[[this.set]]$count)[nrow(datasets[[this.set]]$count)] <- make.names(ifelse(length(tmp$sample.annot$sample.Name)>0,as.character(tmp$sample.annot$sample.Name),nrow(datasets[[this.set]]$count)))
      datasets[[this.set]]$sample.annot <- plyr:::rbind.fill(datasets[[this.set]]$sample.annot ,tmp$sample.annot)
      rownames(datasets[[this.set]]$sample.annot) <- make.names(datasets[[this.set]]$sample.annot$sample.Name)
      
      cat("dataset= ",this.set,", file= ",fl.num,"\n")
      
    }
  }
  
  return(list(datasets = datasets, RCCs = all.RCC.files))
  
}



#==============================================================
# This script is an example of using the funcitonality to gather
# read and collect RCCs into a single file

# if(F){
#   sup.path <- "S:/AMashadi-Hossein/Take2/2017_02_03/data/"
#   allRCC.path <- "S:/AMashadi-Hossein/Take2/2017_02_03/data/allRCC/"
#   
#   Gather.RCC(sup.path = sup.path,allRCC.path = allRCC.path)
#   
#   require(plyr)
#   
#   all.RCC.files <- dir(allRCC.path)
#   all.RCC.files <- all.RCC.files[grep(".RCC",all.RCC.files)]
#   
#   dataset.count <- 1
#   fl <- all.RCC.files[1]
#   #current.Accession <- extract.lane.from.RCC(path = path,fl.name = fl)$count$Accession
#   current.genenames <- extract.lane.from.RCC(allRCC.path = allRCC.path,fl.name = fl)$count$Name
#   count.matrix <- smp.annot <- NULL
#   fl.num <- 0
#   for(fl in all.RCC.files){
#     fl.num <- fl.num + 1
#     tmp <- suppressWarnings(try(extract.lane.from.RCC(allRCC.path = allRCC.path,fl.name = fl),silent = T))
#     if(class(tmp)=="try-error"){
#       cat("File", fl, " Failed\n")
#       next
#     }
#     
#     
#     gene.names.match <- length(tmp$count$Name) == length(current.genenames)
#     if(gene.names.match)
#       gene.names.match <- all(tmp$count$Name == current.genenames)
#     
#     if(gene.names.match){
#       tmpcounts <- tmp$count$Count
#       names(tmpcounts) <- tmp$count$Name
#       count.matrix <-  rbind(count.matrix,tmpcounts)
#       rownames(count.matrix)[nrow(count.matrix)] <- make.names(ifelse(length(tmp$sample.annot$sample.Name)>0,as.character(tmp$sample.annot$sample.Name),nrow(count.matrix)))
#       smp.annot <- plyr:::rbind.fill(smp.annot,tmp$sample.annot)
#       gannot <- tmp$count[c("CodeClass", "Name","Accession")]
#       
#       rownames(smp.annot) <- make.names(smp.annot$sample.Name)
#       rownames(gannot) <- make.names(gannot$Name)
#       
#       cat("dataset= ",dataset.count,", file= ",fl.num,"\n")
#     }else{
#       assign(x = paste("dataset",dataset.count,sep=""),value = list(smp.annot = smp.annot, gannot = gannot,counts = count.matrix),envir = globalenv())
#       dataset.count <- dataset.count+1
#       current.genenames <- tmp$count$Name
#       count.matrix <- smp.annot <- NULL
#       
#       
#       tmpcounts <- tmp$count$Count
#       names(tmpcounts) <- tmp$count$Name
#       count.matrix <-  rbind(count.matrix,tmpcounts)
#       rownames(count.matrix)[nrow(count.matrix)] <- make.names(ifelse(length(tmp$sample.annot$sample.Name)>0,as.character(tmp$sample.annot$sample.Name),nrow(count.matrix)))
#       smp.annot <- plyr:::rbind.fill(smp.annot,tmp$sample.annot)
#       gannot <- tmp$count[c("CodeClass", "Name","Accession")]
#       
#       rownames(smp.annot) <- make.names(smp.annot$sample.Name)
#       rownames(gannot) <- make.names(gannot$Name)
#       
#     }
#   }
#   
#   assign(x = paste("dataset",dataset.count,sep=""),value = list(smp.annot = smp.annot, gannot = gannot,counts = count.matrix),envir = globalenv())
#   
# }
