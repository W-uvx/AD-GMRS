#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("GMRS_ROOT"); source(file.path(root,"scripts","model_utils.R"))
opt <- parse_cli(commandArgs(trailingOnly=TRUE))
if (is.null(opt$phenotype)||is.null(opt$outcome)) abort("--phenotype and --outcome are required")
score_root <- opt$`score-root` %||% file.path(root,"results","runs","training","scores")
model_path <- opt$model %||% file.path(root,"model","GMRS_model.rds")
if(!grepl("\\.rds$",model_path)) abort("--model must end with .rds")
population_path <- opt$`population-profile` %||% Sys.getenv("POPULATION_PROFILE")
if(!nzchar(population_path)) abort("--population-profile or POPULATION_PROFILE is required")
if (!startsWith(population_path,"/")) population_path <- file.path(root,population_path)
population <- read_population(population_path)
manifests <- list(gwas=Sys.getenv("GWAS_MANIFEST"),eqtl=Sys.getenv("EQTL_MANIFEST"),pqtl=Sys.getenv("PQTL_MANIFEST"))
manifests <- lapply(manifests,function(p) if(startsWith(p,"/")) p else file.path(root,p))
evidence <- list(eqtl=Sys.getenv("EQTL_EVIDENCE"),pqtl=Sys.getenv("PQTL_EVIDENCE"))
evidence <- lapply(evidence,function(p) if(nzchar(p)&&!startsWith(p,"/")) file.path(root,p) else p)
epsilon <- as.numeric(Sys.getenv("INTACT_EPSILON","0.1")); gamma <- as.numeric(Sys.getenv("INTACT_GAMMA","1"))
if(!is.finite(epsilon) || epsilon<0 || epsilon>1 || !is.finite(gamma) || gamma<=0) abort("Invalid INTACT_EPSILON or INTACT_GAMMA")
features <- combine_layers(score_root,manifests,population,evidence=evidence,
  epsilon=epsilon,gamma=gamma)
pheno <- key_id(read_tsv(opt$phenotype),"phenotype")
if(!is.null(opt$covariates)) abort("This pipeline trains genotype-only models; --covariates is not supported")
covariates <- character()
required <- c("FID","IID",opt$outcome,covariates)
if(!all(required%in%names(pheno))) abort("Phenotype missing: %s",paste(setdiff(required,names(pheno)),collapse=", "))
dat <- merge(pheno[,required,drop=FALSE],features$data,by=c("FID","IID"),all=FALSE)
if(!nrow(dat)) abort("No FID/IID overlap between phenotype and genetic scores")
message(sprintf("Training on %d of %d scored subjects with phenotype records",nrow(dat),nrow(features$data)))
if(anyNA(dat)) abort("Missing outcome/covariate values are not allowed")
y <- suppressWarnings(as.numeric(dat[[opt$outcome]])); if(anyNA(y) || !setequal(unique(y),c(0,1))) abort("Outcome must contain both 0 and 1")
predictors <- c("PRS","TRS","ProRS",covariates)
x <- model.matrix(~.-1,data=dat[,predictors,drop=FALSE]); center <- colMeans(x); scale <- apply(x,2,sd)
if(any(!is.finite(scale)|scale==0)) abort("Constant or invalid model predictor")
xz <- sweep(sweep(x,2,center,"-"),2,scale,"/")
seed <- suppressWarnings(as.integer(opt$seed %||% 20260810)); alpha <- suppressWarnings(as.numeric(opt$alpha %||% .5))
if(!is.finite(seed) || !is.finite(alpha) || alpha<0 || alpha>1) abort("--seed must be an integer and --alpha must be between 0 and 1")
set.seed(seed); folds <- min(10L,min(table(y)))
if(folds<3) abort("At least 3 cases and 3 controls are required")
foldid <- integer(length(y))
for(cls in c(0,1)) { idx <- sample(which(y==cls)); foldid[idx] <- rep(seq_len(folds),length.out=length(idx)) }
fit <- glmnet::cv.glmnet(xz,y,family="binomial",alpha=alpha,foldid=foldid,keep=TRUE)
oof <- plogis(as.numeric(fit$fit.preval[,which.min(abs(fit$lambda-fit$lambda.1se))]))
object <- list(schema_version="2.0",disease=opt$outcome,engine="cv.glmnet",fit=fit,
  predictors=colnames(xz),covariates=covariates,predictor_center=center,predictor_scale=scale,
  composite_distribution=data.frame(score=names(center),Mean=as.numeric(center),SD=as.numeric(scale)),
  omics_preprocessing=features$preprocessing,population_reference=normalizePath(population_path),
  weight_checksums=read_weight_checksums(score_root),
  training=list(n=nrow(dat),cases=sum(y==1),controls=sum(y==0),prevalence=mean(y)),
  created_at=format(Sys.time(),tz="UTC"),session_info=utils::sessionInfo())
dir.create(dirname(model_path),recursive=TRUE,showWarnings=FALSE); saveRDS(object,model_path)
write.table(object$composite_distribution,sub("\\.rds$",".score_distribution.tsv",model_path),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(FID=dat$FID,IID=dat$IID,observed=y,oof_probability=oof),sub("\\.rds$",".oof.tsv",model_path),sep="\t",quote=FALSE,row.names=FALSE)
message(sprintf("Saved model (%d samples) to %s",nrow(dat),model_path))

