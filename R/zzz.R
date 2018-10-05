#
.onLoad <- function(libname, pkgname) {
  sys.file <- system.file('inst/extdata', 'codon.subst.mat.rds', package = 'pansimulatoR')
  .codon.subst.mat <- readRDS(sys.file)
  assign(".codon.subst.mat", .codon.subst.mat, envir = parent.env(environment()))
}
