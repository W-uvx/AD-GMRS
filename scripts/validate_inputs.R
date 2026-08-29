#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("GMRS_ROOT"); source(file.path(root,"scripts","model_utils.R"))
for(name in c("GWAS_MANIFEST","EQTL_MANIFEST","PQTL_MANIFEST")) {
  path <- Sys.getenv(name); if(!startsWith(path,"/")) path <- file.path(root,path)
  x <- read_tsv(path); required <- c("feature_id","chromosome","summary_file","sample_size")
  if(!all(required%in%names(x))) abort("%s requires: %s",name,paste(required,collapse=", "))
  if(!nrow(x)||anyDuplicated(x$feature_id)||any(!nzchar(x$feature_id))) abort("Empty/duplicate feature IDs in %s",name)
  if(any(grepl("/",x$feature_id,fixed=TRUE) | grepl("\t",x$feature_id,fixed=TRUE))) abort("feature_id cannot contain slash or tab in %s",name)
  bad_chr <- !x$chromosome%in%c(as.character(1:22),"ALL")
  if(any(bad_chr)) abort("Invalid chromosome in %s",name)
  paths <- as.character(x$summary_file); paths[!startsWith(paths,"/")] <- file.path(root,paths[!startsWith(paths,"/")])
  if(any(!file.exists(paths))) abort("Missing summary files in %s: %s",name,paste(x$feature_id[!file.exists(paths)],collapse=", "))
}
profile <- Sys.getenv("POPULATION_PROFILE"); if(!startsWith(profile,"/")) profile <- file.path(root,profile)
invisible(read_population(profile))
for(pkg in c("data.table","glmnet")) if(!requireNamespace(pkg,quietly=TRUE)) abort("Missing R package: %s",pkg)
cat("R inputs and population profile OK\n")
