#
.onLoad <- function(libname, pkgname) {
  sys.file <- system.file('data', 'codon.subst.mat.rds', package = 'pansimulatoR')
  .codon.subst.mat <- readRDS(sys.file)
  assign(".codon.subst.mat", .codon.subst.mat, envir = parent.env(environment()))
}
