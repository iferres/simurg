#' @param ref A fasta file from where to sample genes.
#' @param norg The number of organisms sampled.
#' @param ngenes The number of genes the MRCA will have.
#' @param ne Effective population number. For E coli is ~25 million
#' (Charlesworth et al., 2009). This is also the number of generations.
#' @param theta Gene gain rate per generation.
#' @param rho Gene loss rate per generation.
#' @param mu The per site per generation substitution rate.
#' @param dir_out The non-existent directory where to put the simulated
#' orthologous groups.
#' @return A \code{list} of length 4. (1) The coalescent phylogenetic tree, (2)
#' a \code{data.frame} with the gene gain-loss events, (3) a \code{data.frame}
#' with the substitution events, and (4) the panmatrix. Also a series of
#' attributes are also returned.
#' @author Ignacio Ferres

simpg <- function(ref='pan_genome_reference.fa',
                  norg=10,
                  ngenes=100,
                  ne = 5e7,
                  theta = 2e-6,
                  rho = 1e-6,
                  mu = 5e-10,
                  dir_out='sim_pg'){

  #################
  ## Check input ##
  #################

  if (!file.exists(ref)){
    stop('ref does not exists')
  }

  isnu <- sapply(list(norg=norg,
                      ngenes=ngenes,
                      ne=ne,
                      theta=theta,
                      rho=rho, mu=mu), is.numeric)
  if (!all(isnu)){
    stop(paste('Class numeric required at:', names(isnu)[!isnu]))
  }


  if (dir.exists(dir_out)){
    stop('dir_out already exists')
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
    cat(paste0('Discarded ', length(w), ' sequences due to presence of "n".\n'))
    rf <- rf[-w]
  }

  #Remove genes with length not multiple of 3
  w <- which(sapply(rf, length)%%3!=0)
  if (length(w)>0){
    cat(paste0('Discarded ', length(w), ' sequences due to presence of "n".\n'))
    rf <- rf[-w]
  }

  ##############################
  ## Simulate coalescent tree ##
  ##############################
  cat('Simulating coalescent tree.\n')
  phy <- rcoal(norg, tip.label = paste0('genome', 1:norg))

  m <- as.data.frame(phy$edge)
  m$length <- phy$edge.length

  depth <- coalescent.intervals(phy)$total.depth
  brti <- c(structure(rep(0, norg), names=1:norg), branching.times(phy))

  #################################
  ## Simulate gene gain and loss ##
  #################################

  #' On this step, gene birth and death is simulated in order to obtain a
  #' panmatrix at the end of this stage (IMG model).
  #' norg : number of organisms to sample.
  #' ngenes : number of *starting* genes at the MRCA.
  #' phy : simulated coalescent tree.
  #' theta : gene gain rate per generation.
  #' rho : gene loss rate per generation.

  cat('Simulating gene gain and loss.\n')
  gl <- .sim_gl(phy = phy,
                m = m,
                ne = ne,
                depth = depth,
                brti = brti,
                norg = norg,
                ngenes = ngenes,
                theta = theta,
                rho = rho)

  dfgl <- gl[[1]]
  pm <- gl[[2]]
  dpm <- dim(pm)

  # Sample from the reference many genes as columns in the panmatrix are.
  genes <- sample(rf, size = dpm[2])
  rm(rf)

  #######################
  ## Simulate mutation ##
  #######################

  #' On this step, mutation is simulated according the neutral model. The
  #' output is a set of sequences evolved from a common ancestor for ne
  #' generations at a given substitution rate. Mutations generated avoid stop
  #' codons.
  #' mu : substitution rate per site per generation.

  cat('Simulating point mutations.\n')
  dfmut <- .sim_mut(phy = phy,
                    genes = genes,
                    m = m,
                    depth = depth,
                    brti = brti,
                    ne = ne,
                    norg = norg,
                    ngenes = ngenes,
                    mu = mu)

  ########################
  ## Generate sequences ##
  ########################

  #' Takes the mut data.frame and applies substitutions cronologically to each
  #' gene. Then removes sequences according the final panmatrix.

  rownames(dfgl) <- 1:(dim(dfgl)[1])
  rownames(dfmut) <- 1:(dim(dfmut)[1])

  cat(paste0('Writing groups of orthologous at "', dir_out, '" .\n'))
  dir.create(dir_out)
  for (i in seq_len(dpm[2])){

    ge <- genes[[i]]
    nn <- attr(ge, 'name')
    ge <- paste0(ge[c(T,F,F)],ge[c(F,T,F)],ge[c(F,F,T)])
    ch <- dfmut[which(dfmut$gene==nn),]
    gna <- colnames(pm)[i]
    og <- rep(list(ge), norg)
    names(og) <- paste0(phy$tip.label, '_', gna, ' ref:', nn)

    for (j in seq_len(dim(ch)[1])){
      dsc <- Descendants(phy, ch$to.node[j])[[1]]
      for (k in seq_along(dsc)){
        og[[dsc[k]]][ch$gene.codon.pos[j]] <- ch$change.to[j]
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



  ############
  ## Return ##
  ############

  #' Return a list with a coalescent tree, a data.frame representing the
  #' gain-loss events, another one representing mutation events, and the final
  #' panmatrix. Also some attributes regarding input parameters are returned.

  out <- list(phy, dfgl, dfmut, pm)
  names(out) <- c('coalescent', 'gain-loss', 'substitutions', 'panmatrix')
  attr(out, 'reference') <- normalizePath(ref)
  attr(out, 'norg') <- norg
  attr(out, 'ngenes') <- ngenes
  attr(out, 'ne') <- ne
  attr(out, 'theta') <- theta
  attr(out, 'rho') <- rho
  attr(out, 'mu') <- mu
  attr(out, 'dir_out') <- dir_out
  attr(out, 'class') <- 'pangenomeSimulation'

  cat('DONE\n')
  return(out)

}


