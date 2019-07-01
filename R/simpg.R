#' @name simpg
#' @title Simulate a bacterial pangenome
#' @description Simulate the evolution of a pangenome using a random coalescent
#' tree as guide, according to the neutral model and the Infinite Many Genes
#' (IMG) model.
#' The algorithm first simulates a random coalescent tree (\link[ape]{rcoal}).
#'
#' Then, random deviates from the expected number of gene gain and loss
#' according to branch length and parameters (\code{ggr} and \code{glr}), are
#' used to simulate gene birth and death along the tree. Since the IMG model is
#' used, genes are only transmited vertically from generation to generation.
#'
#' Point mutations are simulated following a similar process as the above:
#' random deviates from the expected number of mutations according to branch
#' length and mutation rate, are used. Mutations are then distributed along
#' sites with uniform probability. By default stop codons are avoided. (note:
#' actually, codons and not sites are conditioned to mutation. See \code{smat}
#' parameter below).
#'
#' Finally, sequences are sampled from \code{ref} and a pangenome is generated
#' following the evolutionary history simulated above.
#' @param ref A multi-fasta file from where to sample genes.
#' @param norg The number of organisms sampled.
#' @param ne Effective population number. For E coli is ~25 million
#' (Charlesworth et al., 2009). This is also the number of generations.
#' @param C The coregenome size.
#' @param u Probability of gene gain per generation.
#' @param v Probability of gene loss rate per generation.
#' @param mu The per site per generation substitution rate.
#' @param dir_out The non-existent directory where to put the simulated
#' orthologous groups.
#' @param smat A codon substitution matrix with probability weights for
#' transition. It must be a \code{data.frame} with \code{dim(smat) == c(64, 64)}, where
#' column and row names are the possible combination of codons (i.e. "aaa",
#' "taa", "caa", etc). The values in the \code{data.frame} correspond to probability
#' weights for obtaining the elements of the vector being sampled. See
#' \link[base]{sample} for more information. The default is a transition matrix
#' with 0 probability of transition from a coding codon to any stop codon and to
#' itself, and equal probability weights of transition to any other codon which
#' differ in, at most, one base (point mutations; one site per codon at most is
#' mutated each time). If you want to use a different codon matrix, I recommend
#' to use the default as reference: \code{simba:::.codon.subst.mat} .
#' See package source (R/subs_mat.R) for a guide on how to generate a similar
#' matrix.
#' @return A \code{list} of length 4. (1) The coalescent phylogenetic tree, (2)
#' a \code{data.frame} with the gene gain-loss events, (3) a \code{data.frame}
#' with the substitution events, and (4) the panmatrix. Also a series of
#' attributes are returned.
#' @importFrom seqinr read.fasta write.fasta
#' @importFrom ape rcoal coalescent.intervals branching.times
#' @importFrom phangorn Descendants
#' @author Ignacio Ferres
#' @export
simpg <- function(ref='pan_genome_reference.fa',
                  norg=10,
                  ne = 1e11,
                  C = 100,
                  u = 1e-8,
                  v = 1e-11,
                  mu = 5e-10,
                  dir_out='sim_pg',
                  smat,
                  force = FALSE,
                  verbose = TRUE){

  #################
  ## Check input ##
  #################

  if (!file.exists(ref)){
    stop('ref does not exists')
  }

  isnu <- sapply(list(norg=norg,
                      C=C,
                      ne=ne,
                      u=u,
                      v=v,
                      mu=mu), is.numeric)
  if (!all(isnu)){
    stop(paste('Class numeric required at:', names(isnu)[!isnu]))
  }


  if (dir.exists(dir_out)){
    stop('dir_out already exists. Use force=TRUE to overwrite it.')
  }

  if (dir.exists(dir_out) & force){
    unlink(dir_out, recursive = TRUE)
    dir.create(dir_out)
  }


  if (!missing(smat)){

    if (!all(dim(smat) == c(64, 64))){
      stop('smat must have dim(smat) == c(64, 64)')
    }

    if (class(smat)!='data.frame'){
      stop('smat must be a data.frame')
    }

    crnames <- c("aaa", "taa", "caa", "gaa", "ata", "tta", "cta", "gta", "aca",
                 "tca", "cca", "gca", "aga", "tga", "cga", "gga", "aat", "tat",
                 "cat", "gat", "att", "ttt", "ctt", "gtt", "act", "tct", "cct",
                 "gct", "agt", "tgt", "cgt", "ggt", "aac", "tac", "cac", "gac",
                 "atc", "ttc", "ctc", "gtc", "acc", "tcc", "ccc", "gcc", "agc",
                 "tgc", "cgc", "ggc", "aag", "tag", "cag", "gag", "atg", "ttg",
                 "ctg", "gtg", "acg", "tcg", "ccg", "gcg", "agg", "tgg", "cgg",
                 "ggg")

    if (!all(colnames(smat)%in%crnames)){
      stop('smat incorrect format: column and row names must be possible codons.
           See ?simpg , or pansimulatoR:::.codon.subst.mat as reference')
    }

    if (!all(rownames(smat)%in%crnames)){
      stop('smat incorrect format: column and row names must be possible codons.
           See ?simpg , or pansimulatoR:::.codon.subst.mat as reference')
    }

  }else{

    smat <- .codon.subst.mat

  }

  ##############################
  ## Load reference sequences ##
  ##############################

  rf <- read.fasta(ref,
                   seqtype = 'DNA',
                   as.string = FALSE,
                   forceDNAtolower = TRUE)

  #Remove genes with 'n's.
  w <- which(sapply(rf, function(x) 'n'%in%x))
  if (length(w)>0){
    if (verbose){
      message(paste0('Discarded ', length(w), ' sequences due to presence of "n".'))
    }
   rf <- rf[-w]
  }

  #Remove genes with length not multiple of 3
  w <- which(sapply(rf, length)%%3!=0)
  if (length(w)>0){
    if (verbose){
      message(paste0('Discarded ', length(w), ' sequences with length not multiple of 3.'))
    }
    rf <- rf[-w]
  }

  ############################
  ## Derived IMG parameters ##
  ############################
  #
  theta <- 2 * ne * u
  rho <- 2 * ne * v
  #theoric mrca accessory size
  mrca_acc_t <- theta / rho
  #observed mrca accessory size
  otheta <- mean(rpois(1000, theta))
  ou <- otheta / (2*ne) #Observed gene gain rate
  orho <- mean(rpois(1000, rho))
  ov <- orho / (2*ne) #Observed gene loss rate
  mrca_acc_o <- round(otheta / orho)

  if (verbose){
    mssg1 <- ' IMG parameters:'
    mssg2 <- paste(' # Core genome size, C =', C)
    mssg3 <- paste(' # Probability of gene gain,', ' \u03C5 =', u)
    mssg4 <- paste(' # Probability of gene loss,', '\u03BD =', v)
    mssg5 <- paste(' # Number of generations, Ne = ', ne)
    mssg6 <- paste(' Derived from parameters:')
    mssg7 <- paste(' # \u03F4 = 2Ne', '\u03C5 =', theta)
    mssg8 <- paste(' # \u03C1 = 2Ne', '\u03BD =', rho)
    mssg9 <- paste(' # Theorical (expected) MRCA size: C + \u03F4 / \u03C1 = ', C+mrca_acc_t)
    mssg10 <- paste(' # Simulated (observed) MRCA size: C + Poi(\u03F4 / \u03C1) = ', C+mrca_acc_o)
    mssg <- c(mssg1, mssg2, mssg3, mssg4, mssg5,
              mssg6, mssg7, mssg8, mssg9, mssg10)
    message(paste(mssg, collapse = '\n'))
  }

  ##############################
  ## Simulate coalescent tree ##
  ##############################
  if (verbose) message('Simulating coalescent tree, \u03C4')
  # cat('Simulating coalescent tree.\n')
  phy <- rcoal(norg, tip.label = paste0('genome', 1:norg))

  m <- as.data.frame(phy$edge)
  m$length <- phy$edge.length

  depth <- coalescent.intervals(phy)$total.depth
  brti <- c(structure(rep(0, norg), names=1:norg), branching.times(phy))
  if (verbose){
    mssg <- paste(capture.output(summary(phy))[-c(1:3)], collapse = '\n')
    message(mssg)
  }

  #################################
  ## Simulate gene gain and loss ##
  #################################

  # On this step, gene birth and death is simulated in order to obtain a
  # panmatrix at the end of this stage (IMG model).
  cat('Simulating gene gain and loss.\n')
  gl <- .sim_gl(phy = phy,
                m = m,
                ne = ne,
                brti = brti,
                norg = norg,
                ou = ou,
                ov = ov,
                mrca_acc_o = mrca_acc_o)

  dfgl <- gl[[1]]
  pm <- gl[[2]]
  rn <- rownames(pm)
  # Attach core genes panmatrix
  pm <- cbind(matrix(1L, nrow = norg, ncol = C), pm)
  dpm <- dim(pm)
  rownames(pm) <- rn
  colnames(pm) <- paste0('gene', seq_len(dpm[2]))

  # Sample from the reference many genes as columns in the panmatrix are.
  genes <- sample(rf, size = dpm[2])
  rm(rf)

  #######################
  ## Simulate mutation ##
  #######################

  # On this step, mutation is simulated according the neutral model. The
  # output is a set of sequences evolved from a common ancestor for ne
  # generations at a given substitution rate. Mutations generated avoid stop
  # codons.
  # mu : substitution rate per site per generation.

  cat('Simulating point mutations.\n')
  mmmut <- .sim_mut(phy = phy,
                    genes = genes,
                    m = m,
                    depth = depth,
                    brti = brti,
                    ne = ne,
                    norg = norg,
                    mu = mu,
                    smat = smat)


  ########################
  ## Generate sequences ##
  ########################

  # Takes the mut data.frame and applies substitutions cronologically to each
  # gene. Then removes sequences according the final panmatrix.
  rownames(dfgl) <- 1:(dim(dfgl)[1])
  dd <- dim(mmmut[[1]])
  rownames(mmmut[[1]]) <- 1:(dd[1])
  change.from <- vector('character', dd[1])
  change.to <- vector('character', dd[1])
  allcodons <- dimnames(smat)[[1]]
  allDes <- Descendants(phy, type = 'tips')

  cat(paste0('Writing groups of orthologous at "', dir_out, '" .\n'))
  dir.create(dir_out)
  for (i in seq_len(dpm[2])){

    ge <- genes[[i]]
    nn <- attr(ge, 'name')
    ge <- paste0(ge[c(T,F,F)],ge[c(F,T,F)],ge[c(F,F,T)])
    ch <- mmmut[[1]][which(mmmut[[2]]==nn),,drop=FALSE]
    rns <- as.integer(rownames(ch))
    gna <- colnames(pm)[i]
    og <- rep(list(ge), norg)
    names(og) <- paste0(phy$tip.label, '_', gna, ' ref:', nn)

    for (j in seq_len(dim(ch)[1])){
      # dsc <- Descendants(phy, ch[j, 3L])[[1]]
      # dsc <- phangorn:::bip(phy)[ch[j, 3L]][[1]]
      dsc <- allDes[[.subset2(ch, j, 3)]]
      rnsj <- .subset2(rns, j)
      chj5 <- .subset2(ch, j, 5L)
      change.from[rnsj] <- og[[dsc[1]]][chj5]
      change.to[rnsj] <- sample(allcodons, 1, prob = .subset2(smat, change.from[rnsj]))
      for (k in seq_along(dsc)){
        og[[dsc[k]]][chj5] <- change.to[rnsj]
      }
    }

    og <- lapply(og, paste0, collapse='')

    aaa <- pm[, i]
    class(aaa) <- 'logical'
    og <- og[paste0(names(which(aaa)), '_', gna, ' ref:', nn)]

    write.fasta(sequences = og,
                names = names(og),
                file.out = paste0(dir_out, '/', gna,'.fasta'))
  }

  dfmut <- as.data.frame(mmmut[[1]])
  dfmut$change.from <- change.from
  dfmut$change.to <- change.to
  dfmut$gene <- mmmut[[2]]

  ############
  ## Return ##
  ############

  # Return a list with a coalescent tree, a data.frame representing the
  # gain-loss events, another one representing mutation events, and the final
  # panmatrix. Also some attributes regarding input parameters are returned.

  out <- list(phy, dfgl, dfmut, pm)
  names(out) <- c('coalescent', 'gain-loss', 'substitutions', 'panmatrix')
  attr(out, 'reference') <- normalizePath(ref)
  attr(out, 'norg') <- norg
  attr(out, 'ngenes') <- ngenes
  attr(out, 'ne') <- ne
  attr(out, 'ggr') <- ggr
  attr(out, 'glr') <- glr
  attr(out, 'mu') <- mu
  attr(out, 'dir_out') <- dir_out
  attr(out, 'class') <- 'pangenomeSimulation'

  cat('DONE\n')
  return(out)

}


