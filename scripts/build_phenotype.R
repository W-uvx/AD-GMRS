#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("GMRS_ROOT")
source(file.path(root,"scripts","model_utils.R"))
opt <- parse_cli(commandArgs(trailingOnly=TRUE))

required <- c("disease","baseline","fam","codes","out")
if(any(vapply(required,function(k)is.null(opt[[k]]),logical(1))))
  abort("Required: --disease --baseline --fam --codes --out")
codes <- strsplit(opt$codes,",",fixed=TRUE)[[1]]
outcome <- opt$`outcome-name` %||% "PHENO"
if(!grepl("^[A-Za-z][A-Za-z0-9_]*$",outcome)) abort("Invalid --outcome-name")

disease <- read_tsv(opt$disease)
if(!all(c("eid","Date","Category")%in%names(disease)))
  abort("Disease table requires eid, Date and Category")
baseline <- data.table::fread(opt$baseline,data.table=FALSE,check.names=FALSE)
if(!"eid"%in%names(baseline)) abort("Baseline table requires eid")
date_cols <- setdiff(names(baseline),"eid")
if(!length(date_cols)) abort("Baseline table requires at least one date column")

first_date <- function(z) {
  z <- as.Date(z[nzchar(z) & !is.na(z)])
  if(length(z)) min(z) else as.Date(NA)
}
baseline$BaselineDate <- as.Date(apply(baseline[,date_cols,drop=FALSE],1,first_date),origin="1970-01-01")
baseline <- baseline[!is.na(baseline$BaselineDate),c("eid","BaselineDate")]
baseline$eid <- as.character(baseline$eid)
disease$eid <- as.character(disease$eid)
disease$Date <- as.Date(disease$Date)

target <- disease[disease$Category%in%codes & !is.na(disease$Date),]
if(nrow(target)) first_target <- aggregate(Date~eid,target,min) else
  first_target <- data.frame(eid=character(),Date=as.Date(character()))
cohort <- merge(baseline,first_target,by="eid",all.x=TRUE)
mode <- opt$mode %||% "prevalent"
if(mode=="prevalent") {
  cohort[[outcome]] <- as.integer(!is.na(cohort$Date) & cohort$Date<=cohort$BaselineDate)
} else if(mode=="ever") {
  cohort[[outcome]] <- as.integer(!is.na(cohort$Date))
} else abort("--mode must be prevalent or ever")

fam <- data.table::fread(opt$fam,header=FALSE,data.table=FALSE,select=1:2)
names(fam) <- c("FID","IID")
fam$FID <- as.character(fam$FID); fam$IID <- as.character(fam$IID)
out <- merge(fam,cohort[,c("eid",outcome)],by.x="IID",by.y="eid",all=FALSE)
out <- out[,c("FID","IID",outcome)]
if(!nrow(out)) abort("No IID overlap between PLINK FAM and phenotype tables")
if(!setequal(unique(out[[outcome]]),c(0,1))) abort("Matched cohort requires both cases and controls")
if(anyDuplicated(out[c("FID","IID")])) abort("Duplicate FID/IID after phenotype construction")
dir.create(dirname(opt$out),recursive=TRUE,showWarnings=FALSE)
write.table(out,opt$out,sep="\t",quote=FALSE,row.names=FALSE)
message(sprintf("Saved %d subjects (%d cases, %d controls) to %s",nrow(out),sum(out[[outcome]]==1),sum(out[[outcome]]==0),opt$out))
