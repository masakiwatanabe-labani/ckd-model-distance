# AI2b: 4 手法すべてを Bioconductor の参照実装で計算する。
#
#   ssGSEA   GSVA::gsva(ssgseaParam(...))   normalize は TRUE / FALSE の両方を出す
#   GSVA     GSVA::gsva(gsvaParam(...))     kcdf は "Gaussian" と既定の "auto" の両方
#   PLAGE    GSVA::gsva(plageParam(...))
#   singscore singscore::simpleScore(rankGenes(...), upSet, knownDirection = TRUE)
#
# R 側は AUC を計算しない。スコア行列だけを書き出し、比較と AUC は Python が行う。

suppressPackageStartupMessages({
  library(GSVA); library(singscore); library(GSEABase)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else getwd()
indir  <- file.path(root, "results", "roundR", "r_input")
outdir <- file.path(root, "results", "roundR")
stopifnot(dir.exists(indir))

read_mat <- function(p) as.matrix(read.delim(p, row.names = 1, check.names = FALSE))
gmt <- getGmt(file.path(indir, "pathways_422.gmt"))
cat("gene sets:", length(gmt), "\n")

MINSZ <- 1L        # 集合は Python 側で 30 遺伝子以上に絞ってある。ここで再度削らない
notes <- character(0)

wr <- function(M, name, tag) {
  M <- as.matrix(M)
  write.table(M, file.path(outdir, sprintf("r_%s_%s.tsv", name, tag)),
              sep = "\t", quote = FALSE, col.names = NA)
  cat(sprintf("  %-22s %d x %d\n", name, nrow(M), ncol(M)))
}

for (tag in c("uncentered", "centered")) {
  M <- read_mat(file.path(indir, sprintf("delta_%s.tsv", tag)))
  cat("\n", tag, ": ", nrow(M), " x ", ncol(M), "\n", sep = "")

  wr(GSVA::gsva(ssgseaParam(M, gmt, minSize = MINSZ, maxSize = Inf,
                            alpha = 0.25, normalize = FALSE), verbose = FALSE),
     "ssgsea", tag)
  wr(GSVA::gsva(ssgseaParam(M, gmt, minSize = MINSZ, maxSize = Inf,
                            alpha = 0.25, normalize = TRUE), verbose = FALSE),
     "ssgseanorm", tag)
  wr(GSVA::gsva(gsvaParam(M, gmt, minSize = MINSZ, maxSize = Inf,
                          kcdf = "Gaussian", tau = 1, maxDiff = TRUE,
                          absRanking = FALSE), verbose = FALSE),
     "gsva", tag)
  auto <- gsvaParam(M, gmt, minSize = MINSZ, maxSize = Inf, kcdf = "auto")
  kc <- tryCatch(as.character(auto@kcdf), error = function(e) "unknown")
  notes <- c(notes, sprintf("gsva_kcdf_auto_%s\t%s", tag, paste(kc, collapse = ",")))
  wr(GSVA::gsva(auto, verbose = FALSE), "gsvaauto", tag)
  wr(GSVA::gsva(plageParam(M, gmt, minSize = MINSZ, maxSize = Inf), verbose = FALSE),
     "plage", tag)

  rk <- singscore::rankGenes(M)
  sc <- sapply(names(gmt), function(nm)
    singscore::simpleScore(rk, upSet = GSEABase::geneIds(gmt[[nm]]),
                           knownDirection = TRUE, centerScore = TRUE)$TotalScore)
  sc <- t(sc); colnames(sc) <- colnames(M)
  wr(sc, "singscore", tag)
}

vers <- c(sprintf("R\t%s", paste(R.version$major, R.version$minor, sep = ".")),
             sprintf("Bioconductor\t%s", as.character(BiocManager::version())),
             sprintf("GSVA\t%s", as.character(packageVersion("GSVA"))),
             sprintf("singscore\t%s", as.character(packageVersion("singscore"))),
             sprintf("GSEABase\t%s", as.character(packageVersion("GSEABase"))),
             "plage_api\tplageParam + gsva (current API)",
             "ssgsea_api\tssgseaParam + gsva (current API)",
             "gsva_api\tgsvaParam + gsva (current API)",
             notes)
writeLines(vers, file.path(outdir, "r_versions.tsv"))

sink(file.path(root, "results", "r_sessioninfo.txt"))
print(sessionInfo())
cat("\nGSVA:", as.character(packageVersion("GSVA")), "\n")
cat("singscore:", as.character(packageVersion("singscore")), "\n")
cat("GSEABase:", as.character(packageVersion("GSEABase")), "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
cat("APIs: ssgseaParam / gsvaParam / plageParam with GSVA::gsva()\n")
sink()
cat("\ndone\n")
