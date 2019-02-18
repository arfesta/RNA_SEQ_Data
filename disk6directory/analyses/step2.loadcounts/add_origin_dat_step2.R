# Combine pedigree origin location with count data
load("~/Projects/RNA_SEQ_Data/shared/Step2_Load_Counts.RData")
library(readxl)
pedigree_2018_04_19_20_36_15 <- read_excel("disk6directory/resources/pedigree/pedigree_2018_04_19_20_36_15.xlsx")
dat.set <- pedigree_2018_04_19_20_36_15[which(pedigree_2018_04_19_20_36_15$Id %in% unique(as.character(load.counts$phenos$fam_id))),]

colnames(dat.set)[1] <- "fam_id"
test.merge <- merge(load.counts$phenos,dat.set,all = T)
no.match <- test.merge[which(is.na(test.merge$State)),]
fs.crosses <- no.match[ grep("x",as.character(no.match$fam_id)),]
fs.cross.names <- strsplit(as.character(fs.crosses$fam_id),split = "x")

fs.cross.names <- lapply(fs.cross.names,function(x){
pedigree_2018_04_19_20_36_15$State[which(pedigree_2018_04_19_20_36_15$Id %in% x)]
})
fs.cross.names <- lapply(fs.cross.names,function(x){
  if(length(x) == 2){
    if(anyNA(x)){
      x[which(is.na(x))] <- "WG"
    }
      paste0(x[1],"x",x[2])
  } else {
      paste0(x,"x","WG")
    } 
})
fs.crosses$State <- unlist(fs.cross.names)

wg.crosses <- no.match[- grep("x",as.character(no.match$fam_id)),]
wg.crosses$State <- "WG"

all.crosses <- rbind(fs.crosses,wg.crosses,test.merge[-which(is.na(test.merge$State)),])

identical(all.crosses$animal_id[match(load.counts$phenos$animal_id,all.crosses$animal_id)],load.counts$phenos$animal_id)
all.crosses <- all.crosses[match(load.counts$phenos$animal_id,all.crosses$animal_id),]

load.counts$phenos <- all.crosses

save(load.counts,file="~/Projects/RNA_SEQ_Data/shared/Step2_Load_Counts_withped.RData",compress=T)
