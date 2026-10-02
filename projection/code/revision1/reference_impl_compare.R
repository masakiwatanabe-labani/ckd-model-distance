# AI2: 参照実装（Bioconductor）で同じスコア行列を作る。
#
# R 側は AUC を計算しない。スコア行列だけを書き出し、比較と AUC は Python が行う。
# GSVA は新 API（plageParam + gsva）を優先し、無ければ旧 API（method="plage"）に落とす。
# どちらを使ったかを必ず記録する。

suppressPackageStartupMessages({
  library(GSVA); library(singscore); library(GSEABase)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else getwd()
indir  <- file.path(root, "results", "roundR", "r_input")
outdir <- file.path(root, "results", "roundR")
stopifnot(dir.exists(indir))

read_mat <- function(p) {
  d <- read.delim(p, row.names = 1, check.names = FALSE)
  as.matrix(d)
}

gmt <- getGmt(file.path(indir, "pathways_422.gmt"))
cat("gene sets:", length(gmt), "\n")

api_used <- NA_character_

score_plage <- function(M) {
  ok <- FALSE
  res <- NULL
  if (exists("plageParam")) {
    res <- tryCatch({
      p <- plageParam(M, gmt, minSize = 1L, maxSize = Inf)
      out <- GSVA::gsva(p, verbose = FALSE)
      api_used <<- "plageParam + gsva (current API)"
      ok <- TRUE
      out
    }, error = function(e) { cat("  plageParam failed:", conditionMessage(e), "\n"); NULL })
  }
  if (!ok) {
    res <- GSVA::gsva(M, gmt, method = "plage", min.sz = 1, max.sz = Inf, verbose = FALSE)
    api_used <<- "gsva(method = \"plage\") (legacy API)"
  }
  res
}

score_singscore <- function(M) {
  rk <- singscore::rankGenes(M)
  sc <- sapply(names(gmt), function(nm) {
    gs <- GSEABase::geneIds(gmt[[nm]])
    s <- singscore::simpleScore(rk, upSet = gs, knownDirection = TRUE)
    s$TotalScore
  })
  t(sc)                                  # set x sample
}

for (tag in c("uncentered", "centered")) {
  M <- read_mat(file.path(indir, sprintf("delta_%s.tsv", tag)))
  cat("\n", tag, ": ", nrow(M), " x ", ncol(M), "\n", sep = "")

  p <- score_plage(M)
  p <- as.matrix(p)
  cat("  PLAGE:", nrow(p), "x", ncol(p), " API:", api_used, "\n")
  write.table(p, file.path(outdir, sprintf("r_plage_%s.tsv", tag)),
              sep = "\t", quote = FALSE, col.names = NA)

  s <- score_singscore(M)
  colnames(s) <- colnames(M)
  cat("  singscore:", nrow(s), "x", ncol(s), "\n")
  write.table(s, file.path(outdir, sprintf("r_singscore_%s.tsv", tag)),
              sep = "\t", quote = FALSE, col.names = NA)
}

writeLines(c(sprintf("GSVA\t%s", as.character(packageVersion("GSVA"))),
             sprintf("singscore\t%s", as.character(packageVersion("singscore"))),
             sprintf("GSEABase\t%s", as.character(packageVersion("GSEABase"))),
             sprintf("R\t%s", paste(R.version$major, R.version$minor, sep = ".")),
             sprintf("plage_api\t%s", api_used)),
           file.path(outdir, "r_versions.tsv"))

sink(file.path(root, "results", "r_sessioninfo.txt"))
print(sessionInfo())
cat("\nGSVA:", as.character(packageVersion("GSVA")), "\n")
cat("singscore:", as.character(packageVersion("singscore")), "\n")
cat("GSEABase:", as.character(packageVersion("GSEABase")), "\n")
cat("PLAGE API used:", api_used, "\n")
sink()
cat("\ndone\n")
