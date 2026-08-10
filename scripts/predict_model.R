#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(glmnet))
root <- Sys.getenv("GMRS_ROOT"); source(file.path(root,"scripts","model_utils.R"))
opt <- parse_cli(commandArgs(trailingOnly=TRUE)); if(is.null(opt$model)||is.null(opt$out)||is.null(opt$`score-root`)) abort("--model, --out and --score-root are required")
model <- readRDS(opt$model); if(!identical(model$schema_version,"2.0")) abort("Unsupported model schema")
if(is.null(model$weight_checksums)) abort("Model does not contain weight provenance")
current_checksums <- read_weight_checksums(opt$`score-root`)
if(!identical(model$weight_checksums,current_checksums)) abort("PRS-CS weights differ from those used to train the model")
if(length(model$covariates)) abort("This genotype-only prediction entry cannot use a model trained with covariates: %s",paste(model$covariates,collapse=","))
if(!identical(model$predictors,c("PRS","TRS","ProRS"))) abort("Model predictors are incompatible with genotype-only prediction")
feat <- combine_layers(opt$`score-root`,manifests=list(),stored=model$omics_preprocessing)$data
x <- as.matrix(feat[,c("PRS","TRS","ProRS"),drop=FALSE]); x <- x[,model$predictors,drop=FALSE]
xz <- sweep(sweep(x,2,model$predictor_center,"-"),2,model$predictor_scale,"/")
prob <- as.numeric(predict(model$fit,xz,s="lambda.1se",type="response"))
raw <- feat; names(raw)[match(c("PRS","TRS","ProRS"),names(raw))] <- c("PRS_raw","TRS_raw","ProRS_raw")
standardised <- data.frame(PRS=xz[,"PRS"],TRS=xz[,"TRS"],ProRS=xz[,"ProRS"])
out <- cbind(raw,standardised,linear_predictor=qlogis(pmin(pmax(prob,1e-12),1-1e-12)),disease_probability=prob)
dir.create(dirname(opt$out),recursive=TRUE,showWarnings=FALSE); write.table(out,opt$out,quote=FALSE,sep="\t",row.names=FALSE)
