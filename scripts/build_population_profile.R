#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
get <- function(k) { i <- match(k, args); if (is.na(i) || i == length(args)) stop("Missing ", k); args[i + 1] }
x <- fread(get("--score-matrix"), data.table = FALSE)
sample_id <- intersect(c("IID", "eid", "ID"), names(x))
if (!length(sample_id)) stop("Score matrix requires IID, eid or ID")
id <- intersect(c("FID", "IID", "eid", "ID"), names(x))
values <- x[, setdiff(names(x), id), drop = FALSE]
if (!ncol(values)) stop("Score matrix contains no score columns")
if (!all(vapply(values, is.numeric, logical(1)))) stop("All score columns must be numeric")
out <- data.frame(File = names(values), Mean = vapply(values, mean, numeric(1), na.rm = TRUE),
                  SD = vapply(values, sd, numeric(1), na.rm = TRUE),
                  Min = vapply(values, min, numeric(1), na.rm = TRUE),
                  Max = vapply(values, max, numeric(1), na.rm = TRUE),
                  N = vapply(values, function(z) sum(is.finite(z)), integer(1)))
if (any(out$N < 2) || any(!is.finite(out$Mean)) || any(!is.finite(out$SD)) || any(out$SD <= 0))
  stop("Every score requires at least two finite observations and a positive SD")
out_path <- get("--out")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, quote = FALSE, sep = "\t", row.names = FALSE)

