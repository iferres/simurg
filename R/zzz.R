#
.onLoad <- function(libname, pkgname) {
  sys.file <- system.file('extdata', 'codon.subst.mat.rds', package = 'simurg')
  .codon.subst.mat <- readRDS(sys.file)
  assign(".codon.subst.mat", .codon.subst.mat, envir = parent.env(environment()))
}
