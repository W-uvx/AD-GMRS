#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("GMRS_ROOT"); source(file.path(root,"scripts","model_utils.R"))
opt <- parse_cli(commandArgs(trailingOnly=TRUE))
if(is.null(opt$model)||is.null(opt$out)||is.null(opt$`score-root`)) abort("--model, --out and --score-root are required")
model <- readRDS(opt$model)
if(!identical(model$schema_version,"2.0")) abort("Unsupported model schema")
if(is.null(model$omics_preprocessing)) abort("Model does not contain omics preprocessing metadata")
if(is.null(model$predictors)||anyDuplicated(model$predictors)) abort("Model has invalid predictors")
if(!all(c("PRS","TRS","ProRS")%in%model$predictors)) abort("Model predictors do not contain PRS, TRS and ProRS")
if(!is.null(model$weight_checksums)) {
  current_checksums <- read_weight_checksums(opt$`score-root`)
  if(!identical(model$weight_checksums,current_checksums)) abort("PRS-CS weights differ from those used to train the model")
}

feat <- combine_layers(opt$`score-root`,manifests=list(),stored=model$omics_preprocessing)$data
covariates <- model$covariates %||% setdiff(model$predictors,c("PRS","TRS","ProRS"))
if(length(covariates)) {
  if(is.null(opt$covariates)) abort("--covariates is required; model expects: %s",paste(covariates,collapse=", "))
  covar <- key_id(read_tsv(opt$covariates),"covariates")
  missing_columns <- setdiff(covariates,names(covar))
  if(length(missing_columns)) abort("Covariate file is missing: %s",paste(missing_columns,collapse=", "))
  score_ids <- paste(feat$FID,feat$IID,sep="\r")
  covar_ids <- paste(covar$FID,covar$IID,sep="\r")
  missing_ids <- !score_ids%in%covar_ids
  extra_ids <- !covar_ids%in%score_ids
  if(any(missing_ids)||any(extra_ids))
    abort("FID/IID sets differ between genotype scores and covariates (missing=%d, extra=%d)",sum(missing_ids),sum(extra_ids))
  covar_order <- match(score_ids,covar_ids)
  dat <- cbind(feat,covar[covar_order,covariates,drop=FALSE])
} else dat <- feat

x <- dat[,model$predictors,drop=FALSE]
if(!all(vapply(x,is.numeric,logical(1)))) abort("All model predictors must be numeric")
if(anyNA(x)||any(!is.finite(as.matrix(x)))) abort("Model predictors contain missing or non-finite values")

if(identical(model$engine,"caret::train/glmnet")) {
  if(!requireNamespace("caret",quietly=TRUE)) abort("R package caret is required for this model")
  if(!inherits(model$fit,"train")) abort("Model fit is not a caret train object")
  if(is.null(model$positive_class)||!model$positive_class%in%model$fit$levels) abort("Model has an invalid positive_class")
  prob_table <- predict(model$fit,newdata=x,type="prob")
  if(!model$positive_class%in%names(prob_table)) abort("Positive class is absent from caret probability output")
  prob <- as.numeric(prob_table[[model$positive_class]])
  processed <- if(is.null(model$fit$preProcess)) x else predict(model$fit$preProcess,x)
} else if(identical(model$engine,"cv.glmnet")) {
  if(!requireNamespace("glmnet",quietly=TRUE)) abort("R package glmnet is required for this model")
  if(is.null(model$predictor_center)||is.null(model$predictor_scale))
    abort("cv.glmnet model lacks predictor standardisation parameters")
  processed <- as.data.frame(sweep(sweep(as.matrix(x),2,model$predictor_center,"-"),2,model$predictor_scale,"/"))
  prob <- as.numeric(predict(model$fit,as.matrix(processed),s="lambda.1se",type="response"))
} else abort("Unsupported model engine: %s",model$engine %||% "missing")

raw <- feat
names(raw)[match(c("PRS","TRS","ProRS"),names(raw))] <- c("PRS_raw","TRS_raw","ProRS_raw")
standardised <- processed[,c("PRS","TRS","ProRS"),drop=FALSE]
out <- cbind(raw,standardised,
  linear_predictor=qlogis(pmin(pmax(prob,1e-12),1-1e-12)),
  disease_probability=prob)
dir.create(dirname(opt$out),recursive=TRUE,showWarnings=FALSE)
write.table(out,opt$out,quote=FALSE,sep="\t",row.names=FALSE)
