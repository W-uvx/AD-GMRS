abort <- function(...) stop(sprintf(...), call. = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x

parse_cli <- function(args) {
  out <- list(); i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--")) abort("Unexpected argument: %s", args[[i]])
    key <- sub("^--", "", args[[i]])
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) out[[key]] <- TRUE else {
      out[[key]] <- args[[i + 1L]]; i <- i + 1L
    }
    i <- i + 1L
  }
  out
}

read_tsv <- function(path) {
  if (is.null(path) || !file.exists(path)) abort("File not found: %s", path %||% "NULL")
  data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

key_id <- function(x, context) {
  if (!all(c("FID", "IID") %in% names(x))) abort("%s requires FID and IID columns", context)
  x$FID <- as.character(x$FID); x$IID <- as.character(x$IID)
  if (anyDuplicated(x[c("FID", "IID")])) abort("Duplicate FID/IID in %s", context)
  x
}

read_profile <- function(path, feature) {
  x <- key_id(read_tsv(path), path)
  candidates <- c("SCORE", "SCORESUM", "SCORE1_SUM", "SCORE1_AVG")
  score <- intersect(candidates, names(x))[1]
  if (is.na(score)) {
    numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
    score <- tail(setdiff(numeric_cols, c("PHENO", "CNT", "CNT2")), 1)
  }
  if (length(score) != 1L || is.na(score)) abort("Cannot identify PLINK score column in %s", path)
  setNames(data.frame(FID=x$FID, IID=x$IID, value=as.numeric(x[[score]]), check.names=FALSE),
           c("FID", "IID", feature))
}

read_layer <- function(score_root, layer, expected) {
  directory <- file.path(score_root, layer)
  paths <- file.path(directory, paste0(expected, ".profile"))
  missing <- expected[!file.exists(paths)]
  if (length(missing)) abort("Missing %s score files: %s", layer, paste(missing, collapse=", "))
  extras <- setdiff(sub("\\.profile$", "", list.files(directory, pattern="\\.profile$")), expected)
  if (length(extras)) message(sprintf("Ignoring %d stale %s profiles", length(extras), layer))
  parts <- Map(read_profile, paths, expected)
  Reduce(function(a,b) merge(a,b,by=c("FID","IID"),all=TRUE,sort=FALSE), parts)
}

read_manifest_features <- function(path) {
  x <- read_tsv(path)
  if (!"feature_id" %in% names(x) || anyDuplicated(x$feature_id)) abort("Invalid/duplicate feature_id in %s", path)
  as.character(x$feature_id)
}

read_population <- function(path) {
  x <- read_tsv(path)
  if (!all(c("File","Mean","SD") %in% names(x))) abort("Population profile requires File, Mean and SD")
  if (anyDuplicated(x$File)) abort("Duplicate features in population profile")
  x$Mean <- as.numeric(x$Mean); x$SD <- as.numeric(x$SD)
  if (any(!is.finite(x$Mean)) || any(!is.finite(x$SD) | x$SD < 0)) abort("Population Mean/SD must be finite and SD >= 0")
  if (all(c("Min","Max") %in% names(x))) {
    x$Min <- as.numeric(x$Min); x$Max <- as.numeric(x$Max)
    if (any(!is.finite(x$Min)) || any(!is.finite(x$Max))) abort("Population Min/Max must be finite")
    extreme <- rep(0,nrow(x)); valid <- x$SD > 0
    extreme[valid] <- pmax(abs((x$Min[valid]-x$Mean[valid])/x$SD[valid]), abs((x$Max[valid]-x$Mean[valid])/x$SD[valid]))
    if (any(extreme > 100, na.rm=TRUE))
      abort("Population profile is inconsistent: %d features exceed 100 SD; SD values may still be SE", sum(extreme>100,na.rm=TRUE))
  }
  x
}

evidence_weights <- function(features, evidence_path="", epsilon=.1, gamma=1) {
  if (!nzchar(evidence_path)) return(setNames(rep(1/length(features),length(features)),features))
  ev <- read_tsv(evidence_path)
  if (!all(c("feature_id","PP4","Z") %in% names(ev))) abort("Evidence requires feature_id, PP4 and Z")
  if (anyDuplicated(ev$feature_id)) abort("Duplicate feature_id in evidence: %s", evidence_path)
  ev <- ev[match(features,ev$feature_id),]
  if (anyNA(ev$feature_id) || any(!is.finite(ev$PP4)) || any(!is.finite(ev$Z))) abort("Evidence does not fully cover selected features")
  if (!requireNamespace("INTACT",quietly=TRUE)) abort("INTACT is required when an evidence table is configured")
  prob <- as.numeric(INTACT::intact(GLCP_vec=ev$PP4,z_vec=ev$Z))
  if (length(prob)!=length(features) || any(!is.finite(prob))) abort("INTACT returned invalid probabilities")
  raw <- (epsilon+(1-epsilon)*prob)^gamma
  setNames(raw/sum(raw),features)
}

read_weight_checksums <- function(score_root) {
  rows <- list()
  for(layer in c("gwas","eqtl","pqtl")) {
    path <- file.path(score_root,layer,"weight_checksums.tsv")
    x <- read_tsv(path)
    if(!all(c("feature_id","sha256")%in%names(x)) || anyDuplicated(x$feature_id)) abort("Invalid weight checksum table: %s",path)
    x$key <- paste(layer,x$feature_id,sep="::")
    rows[[layer]] <- x[,c("key","sha256")]
  }
  out <- do.call(rbind,rows)
  setNames(as.character(out$sha256),out$key)
}

aggregate_layer <- function(x, layer, population, evidence="", epsilon=.1, gamma=1, stored=NULL) {
  observed <- setdiff(names(x),c("FID","IID"))
  expected <- if (is.null(stored)) observed else stored$features
  if (!setequal(observed,expected)) abort("%s features differ from fitted model; missing=[%s], extra=[%s]", layer,
    paste(setdiff(expected,observed),collapse=","),paste(setdiff(observed,expected),collapse=","))
  mat <- as.matrix(x[,expected,drop=FALSE]); storage.mode(mat) <- "double"
  if (anyNA(mat)) abort("Missing sample scores in %s layer after FID/IID merge", layer)
  if (is.null(stored)) {
    idx <- match(expected,population$File)
    if (anyNA(idx) && layer != "gwas") abort("Population profile misses %s features: %s", layer,paste(expected[is.na(idx)],collapse=","))
    if (layer == "gwas" && anyNA(idx)) {
      message("GWAS population parameters absent; estimating Mean/SD from the training cohort")
      center <- setNames(colMeans(mat),expected); scale <- setNames(apply(mat,2,sd),expected)
    } else {
      if (any(population$SD[idx] <= 0)) abort("Selected %s features have zero SD: %s",layer,paste(expected[population$SD[idx]<=0],collapse=","))
      center <- setNames(population$Mean[idx],expected); scale <- setNames(population$SD[idx],expected)
    }
    if (any(!is.finite(scale) | scale <= 0)) abort("Cannot obtain positive population SD")
    weights <- evidence_weights(expected,evidence,epsilon,gamma)
  } else { center <- stored$center; scale <- stored$scale; weights <- stored$weights }
  if (!identical(names(center),expected) || !identical(names(scale),expected) || !identical(names(weights),expected)) abort("Corrupt preprocessing metadata for %s",layer)
  z <- sweep(sweep(mat,2,center,"-"),2,scale,"/")
  name <- switch(layer,gwas="PRS",eqtl="TRS",pqtl="ProRS")
  list(data=setNames(data.frame(FID=x$FID,IID=x$IID,score=as.numeric(z%*%weights)),c("FID","IID",name)),
       preprocessing=list(features=expected,center=center,scale=scale,weights=weights))
}

combine_layers <- function(score_root, manifests, population=NULL, stored=NULL, evidence=list(), epsilon=.1, gamma=1) {
  result <- NULL; prep <- list(); reference_ids <- NULL
  for (layer in c("gwas","eqtl","pqtl")) {
    expected <- if (is.null(stored)) read_manifest_features(manifests[[layer]]) else stored[[layer]]$features
    raw <- read_layer(score_root,layer,expected)
    ids <- paste(raw$FID, raw$IID, sep="\r")
    if (is.null(reference_ids)) reference_ids <- ids
    else if (!setequal(reference_ids, ids)) abort("FID/IID sets differ across score layers; mismatch found in %s", layer)
    ag <- aggregate_layer(raw,layer,population,evidence[[layer]] %||% "",epsilon,gamma,if(is.null(stored)) NULL else stored[[layer]])
    result <- if(is.null(result)) ag$data else merge(result,ag$data,by=c("FID","IID"),all=FALSE,sort=FALSE)
    prep[[layer]] <- ag$preprocessing
  }
  list(data=result,preprocessing=prep)
}
