
#' @export
print.pangenomeSimulation <- function(x, ...){

  attrs <- attributes(x)
  # theta <- 2 * attrs$ne * attrs$u
  # rho <- 2 * attrs$ne * attrs$v
  cat('Object of class "pangenomeSimulation".\n')
  cat(paste(' Number of sampled organisms: ', attrs$norg, '\n'))
  cat(paste(' Effective Population Size: ', attrs$ne, '\n'))
  cat(' IMG parameters:\n')
  cat(paste('   Coregenome size: ', attrs$C, '\n'))
  cat(paste('   Probability of gene gain, per generation:  \u03C5 =', attrs$u, '\n'))
  cat(paste('   Probability of gene loss, per generation:  \u03BD =', attrs$v, '\n'))
  cat(' Substitution parameters:\n')
  cat(paste('   Probability of substitution, per generation:   \u03BC =', attrs$mu, '\n'))
  cat(paste(' Fasta reference file at', attrs$ref, '\n'))
  cat(paste(' Sequences written at', attrs$dir_out, '\n'))

}




#' @export
summary.pangenomeSimulation <- function(object, ...){

  cat('Summary of object of class pangenomeSimulation:\n')
  attrs <- attributes(object)
  ret <- list(Gene_family_frecuency = NULL, Evo_dist = NULL)

  phy <- object$coalescent
  summary(object$coalescent)
  cat('\n')

  cat('** Infinitely Many Genes Model ** \n')

  norg <- attrs$norg
  ne <- attrs$ne
  C <- attrs$C
  u <- attrs$u
  v <- attrs$v
  otheta <- attrs$otheta
  orho <- attrs$orho
  mu <- attrs$mu

  pm <- object$panmatrix
  dpm <- dim(pm)
  rs <- rowSums(pm)
  cs <- colSums(pm)


  theta <- 2 * ne * u
  rho <- 2 * ne * v
  #theoric mrca accessory size
  mrca_acc_t <- theta / rho
  #observed mrca accessory size
  mrca_acc_o <- round(otheta / orho)
  EG <- C + theta * sum(1/(rho + 0:(norg-1)))

  mssg1 <- ' Parameters:'
  mssg2 <- paste(' # Core genome size, C =', C)
  mssg3 <- paste(' # Probability of gene gain,', ' \u03C5 =', u)
  mssg4 <- paste(' # Probability of gene loss,', '\u03BD =', v)
  mssg5 <- paste(' # Number of generations, Ne = ', ne)
  mssg6 <- paste(' Derived from parameters:')
  mssg7 <- paste(' # \u03F4 = 2Ne', '\u03C5 =', theta)
  mssg8 <- paste(' # \u03C1 = 2Ne', '\u03BD =', rho)
  mssg9 <- paste(' # Theorical MRCA size: C + \u03F4 / \u03C1 = ', C+mrca_acc_t)
  mssg10 <- paste(' # Simulated MRCA size: C + Poi(\u03F4 / \u03C1) = ', C+mrca_acc_o)
  mssg11 <- paste(' # Expected pangenome size: E[G] = C + \u03F4 * sum( 1 / (\u03C1 + 0:',norg-1,') ) = ',format(EG), sep = '')
  mssg12 <- paste(' # Observed pangenome size: dim( *$panmatrix )[2] =', dpm[2])
  mssg13 <- paste(' # Expected average number of genes per genome: E[A] = C + ( \u03F4 / \u03C1 ) = ',  C+mrca_acc_t)
  mssg14 <- paste(' # Observed number of genes per genome: rowSums( *$panmatrix ) =\n')
  mssg <- c(mssg1, mssg2, mssg3, mssg4, mssg5,
            mssg6, mssg7, mssg8, mssg9, mssg10,
            mssg11, mssg12, mssg13, mssg14)
  cat(paste(mssg, collapse = '\n'))
  print(rs)
  cat('$Gene_family_frequency\n')
  fam_frec <- table(cs, dnn = NULL)
  ret$IMG <- fam_frec
  print(fam_frec)


  cat('\n')


  cat('** Evolutionary Genetic Model **\n')

  mu <- attrs$mu
  mssg1 <- ' Parameters:'
  mssg2 <- paste(' # Probability of substitution,  \u03BC =', mu)
  mssg3 <- paste(' # Number of generations, Ne = ', ne)
  mssg <- c(mssg1, mssg2, mssg4)

  cat(paste(mssg, collapse = '\n'))
  cat('\n')

  muts <- object$substitutions
  distances <- object$distances

  comb2 <- combn(rownames(object$panmatrix), 2)
  md <- apply(comb2, 2, function(x){
    x1 <- x[1]
    x2 <- x[2]
    cor <- which(colSums(pm[x, ])==2)
    mean(sapply(cor, function(y) .subset2(as.matrix(distances[[y]]),x1, x2) ))
  })

  df <- as.data.frame(t(comb2), stringsAsFactors = FALSE)
  colnames(df) <- c('G1', 'G2')
  coph <- cophenetic(phy)
  # Normalization of cophenetic distances
  coph <- coph / max(coph)
  df$norm_cophenetic <- apply(df, 1, function(x){coph[x[1], x[2]]})
  df$mean_gene_dist <- md
  ret$Substitutions <- df
  cat('$Evo_dist \n')
  str(df)



  return(invisible(ret))
}
