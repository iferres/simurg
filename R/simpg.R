#' @name simpg
#' @title Simulate a bacterial pangenome
#' @description Simulate the evolution of a pangenome using a random coalescent
#' tree as guide, according to the neutral model and the Infinitely Many Genes
#' (IMG) model.
#' The algorithm first simulates a random coalescent tree (\link[ape]{rcoal})
#' describing the evolutionary history of organisms sampled at the final time.
#' Effective population size (ne) is also interpreted as the number of
#' generations from the MRCA to the sampling time.
#'
#' Then, a IMG process is simulated using this tree, and parameters C (coregenome
#' size), u (probability of gene gain, per generation), and v (probability of
#' gene loss, per generation). At the end of this process a final panmatrix is
#' obtained, describing the presence or absence of each gene for each organism,
#' at the sampling time. The final gene presence/absence pattern is in
#' accordance with the previously simulated coalescent tree, as the model
#' describes.
#'
#' As many genes as in the above panmatrix are sampled from the reference
#' file (ref). For each one the MRCA node in the tree is identified, and from
#' this node to the final time the gene is subjected to mutations at rate mu
#' (probability of substitution, per generation). Genes are first replicated
#' as many times it appears at the final time, and mutations are computed in
#' accordance to the simulated tree (each gene will follow an independent path).
#' This implementation mutates codons, rather than single nucleotides, because
#' we wanted to be able to avoid stop codons of being created at the middle of
#' a gene. The probability of substitution from one codon to another is decribed
#' by the \code{smat} parameter, which by defaults assigns equal probability of
#' changing one codon to another with a single nucleotide difference between
#' them, 0 probability of changing from one codon to another with more than one
#' difference, and 0 probability of changing from any codon to a stop codon.
#' This behavior can be changed by modifying smat parameter, which by default
#' takes the substitution matrix \code{'.codon.subst.mat'} in the namespace of
#' this package. To see a reference of how it was created, the file subs_mat.R
#' in the package source shows details, and can be taken as guide to create a
#' custom one.
#'
#' Finally, sequences are written to the output directory (dir_out) by gene
#' (group of orthologous sequences), or by genome, in fasta format. Each header
#' contains information about which gene in the reference file was used to
#' generate it. One can set the algorithm to be able to sample genes from the
#' reference file with replace (replace = TRUE), but it is not recommended
#' since the IMG model states that each gene only is gained once in the
#' evolutionary history of the pan organism (so replace = FALSE, by default,
#' and gives error if number of columns in the final panmatrix is greater than
#' available genes in ref), although if it is set to TRUE, each copy will
#' follow an independent evolutionary path.
#'
#' @param ref A multi-fasta file from where to sample genes.
#' @param norg The number of organisms sampled.
#' @param ne Effective population number. Baumdicker et al. (2012) estimates
#' the effective population sizes of Prochlorococcus and Synechococcus to be
#' around 10e11, which was taken as default. This is also the number of
#' generations from the MRCA to the sampled (final) organisms.
#' @param C The coregenome size (default is 100).
#' @param u Probability of gene gain per generation. Default is 1e-8, to be in
#' accordance to values presented in Baumdicker et al. (2012). (theta = 2Ne*u,
#' and taking Ne = 10e11, and estimates of theta ~2000 for Prochlorococcus).
#' @param v Probability of gene loss rate per generation. Default is 1e-11, to
#' be in accordance to values presented in Baumdicker et al. (2012). (rho = 2Ne*v,
#' and taking Ne = 10e11, and estimates of rho ~2 for Prochlorococcus).
#' @param mu The per site per generation substitution rate. According to Duchene
#' et al. (2016), a typical substitution rate is around 1e-5 and 1e-9 per site,
#' per annum. Taking 1e-9, and assuming 200 generation per annum (Shierup and Wuif,
#' 2010), this gives mu = 5e-12, which was taken as default.
#' @param write_by One of \code{'gene'} of \code{'genome'}. Whether to write
#' sequences in files grouped by orthology or by genome, respectively.
#' @param dir_out The non-existent directory where to put the simulated
#' orthologous groups.
#' @param replace \code{logical}. Whether to replace (see \link[base]{sample}) when sampling
#' genes from reference. Default and recommended is FALSE, since IMG model states
#' that each new gene is a different one. Setting \code{replace = TRUE} will probably
#' cause duplicated genes, although each one will follow a different evolutionary
#' history in the simulation.
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
#' @param force \code{logical}. If \code{dir_out} exists, whether to force
#' directory creation by overwriting the existing one. It will eliminate any
#' existing content.
#' @param verbose \code{logical}. Show (or not) progress messages.
#' @return A \code{list} of length 4. (1) The coalescent phylogenetic tree, (2)
#' a \code{data.frame} with the gene gain-loss events, (3) a \code{data.frame}
#' with the substitution events, and (4) the panmatrix. Also a series of
#' attributes are returned.
#' @importFrom seqinr read.fasta write.fasta
#' @importFrom stats rpois setNames
#' @importFrom ape rcoal coalescent.intervals branching.times
#' @importFrom phangorn Descendants
#' @importFrom utils txtProgressBar setTxtProgressBar capture.output str
#' @author Ignacio Ferres
#' @export
simpg <- function(ref='pan_genome_reference.fa',
                  norg=10,
                  ne = 1e11,
                  C = 100,
                  u = 1e-8,
                  v = 1e-11,
                  mu = 5e-12,
                  write_by = 'gene',
                  dir_out='sim_pg',
                  replace = FALSE,
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

  write_by <- match.arg(write_by, choices = c('gene', 'genome'), several.ok = F)


  if (dir.exists(dir_out)){
    stop('dir_out already exists. Use force=TRUE to overwrite it.')
  }else{
    dir.create(dir_out)
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
  pm <- .sim_gl(phy = phy,
                m = m,
                ne = ne,
                brti = brti,
                norg = norg,
                ou = ou,
                ov = ov,
                mrca_acc_o = mrca_acc_o)

  # dfgl <- gl[[1]]
  # pm <- gl[[2]]
  rn <- rownames(pm)
  # Attach core genes panmatrix
  pm <- cbind(matrix(1L, nrow = norg, ncol = C), pm)
  dpm <- dim(pm)
  rownames(pm) <- rn
  colnames(pm) <- paste0('gene', seq_len(dpm[2]))

  # Sample from the reference many genes as columns in the panmatrix are.
  if (dpm[2] > length(rf) & replace == FALSE){
    stop("cannot take a sample larger than the panmatrix when 'replace = FALSE'")
  }
  if (replace){
    genes <- sample(rf, size = dpm[2], replace = TRUE)
  }else{
    genes <- sample(rf, size = dpm[2], replace = FALSE)
  }

  rm(rf)
  nsites <- sapply(genes, length)

  #######################
  ## Simulate mutation ##
  #######################

  # On this step, mutation is simulated according the neutral model. The
  # output is a set of sequences evolved from a common ancestor for ne
  # generations at a given substitution rate. Mutations generated avoid stop
  # codons.
  # mu : substitution rate per site per generation.

  if (verbose){
    mssg1 <- ' Simulating point mutations:'
    mssg2 <- paste(' # Mutation rate, \u03BC =', mu)
    mssg <- c(mssg1, mssg2)
    message(paste(mssg, collapse = '\n'))
  }

  mmmut <- .sim_mut(pm = pm,
                    phy = phy,
                    nsites = nsites,
                    m = m,
                    depth = depth,
                    brti = brti,
                    ne = ne,
                    norg = norg,
                    mu = mu,
                    smat = smat)

  if (verbose){
    mssg1 <- 'Number and position of substitutions branch-wise simulated:\n'
    mssg2 <- paste(capture.output(str(mmmut, list.len = 3)), collapse = '\n')
    mssg <- c(mssg1, mssg2)
    message(mssg)
  }


  ########################
  ## Generate sequences ##
  ########################

  allcodons <- dimnames(smat)[[1]]
  allDes <- Descendants(phy, type = 'tips')
  # Paths from root to tips
  roottotip <- setNames(ape::nodepath(phy), nm = seq_len(norg))
  gnas <- colnames(pm)

  if (verbose){
    mssg <- paste('Writing sequences at', dir_out)
    message(mssg)
    pb <- txtProgressBar(min = 0, max = dpm[2], style = 3)
  }

  for (i in seq_len(dpm[2])){
    gnms <- names(which(pm[, i]==1))
    ngnm <- length(gnms)
    mrca <- getMRCA(phy, tip = gnms)
    # alds <- allDes[[mrca]]
    ge <- genes[[i]]
    nn <- attr(ge, 'name')
    muts <- mmmut[[nn]]
    # Unique, since more than 2 substitutions on the same place is
    # seen as a single change. Saturation:
    muts <- lapply(muts, unique)
    gna <- gnas[i]

    ge <- paste0(ge[c(T,F,F)],ge[c(F,T,F)],ge[c(F,F,T)])
    og <- rep(list(ge), ngnm)
    # names(og) <- paste0(gnms, '_', gna, ' ref:', nn)
    names(og) <- gnms

    if (!is.null(mrca)){

      # Segments.
      ismrca <- which(sapply(roottotip, function(x) mrca %in%x))
      rttp <- lapply(roottotip[ismrca], function(x){
        x[which(x==mrca):length(x)]
      })

      sgmt <- lapply(rttp, function(x){
        paste(x[1:(length(x)-1)], x[2:length(x)], sep = '-')
      })

    }else{
      # Segments
      sgmt <- paste(rev(rev(roottotip[[which(pm[, i]==1)]])[1:2]), collapse = '-')
    }

    # Iterate over segments
    sgmt <- unique(unlist(sgmt, use.names = FALSE))
    for (j in sgmt){
      #Get mutation position (in codon coordinates)
      mt <- as.integer(ceiling(muts[[j]]/3))
      #Get sequences to apply changes
      dss <- allDes[[as.integer(strsplit(j, '-')[[1]][2])]]
      orgs <- phy$tip.label[dss] # orgs are all organisms which descend from this node
      geno <- gnms[which(gnms%in%orgs)] # gnms are all organisms which DO HAVE this gene
      # If length(geno)==0 means that none of the descendants genomes in this
      # branch inherit that gene (gene loss on that branch). If it is the case,
      # do not compute anything.
      if (length(geno)){
        # Iterate over mutations
        for (k in mt){
          old <- og[[geno[1]]][k]
          new <- sample(allcodons, 1, prob = .subset2(smat, old))
          # Apply changes to gene of descendant organisms
          for (o in geno){
            og[[o]][k] <- new
          }
        }
      }
    }

    og <- lapply(og, paste0, collapse = '')
    gnn <- names(og)
    names(og) <- paste0(gnn, '_', gna, ' ref:', nn)

    if (write_by == 'gene'){
      write.fasta(sequences = og,
                  names = names(og),
                  file.out = paste0(dir_out, '/', gna,'.fasta'))
    }else{

      lapply(names(og), function(x){
        sspl <- strsplit(x, '_')[[1]][1]
        fn <- paste0(sspl, '.fasta')
        write.fasta(sequences = og[x],
                    names = x,
                    file.out = paste0(dir_out, '/', fn),
                    open = 'a')

      })

    }

    if (verbose) setTxtProgressBar(pb, i)

  }

  if (verbose) close(pb)




  ############
  ## Return ##
  ############









  # ########################
  # ## Generate sequences ##
  # ########################
  #
  # # Takes the mut data.frame and applies substitutions cronologically to each
  # # gene. Then removes sequences according the final panmatrix.
  # rownames(dfgl) <- 1:(dim(dfgl)[1])
  # dd <- dim(mmmut[[1]])
  # rownames(mmmut[[1]]) <- 1:(dd[1])
  # change.from <- vector('character', dd[1])
  # change.to <- vector('character', dd[1])
  # allcodons <- dimnames(smat)[[1]]
  # allDes <- Descendants(phy, type = 'tips')
  #
  # cat(paste0('Writing groups of orthologous at "', dir_out, '" .\n'))
  # dir.create(dir_out)
  # for (i in seq_len(dpm[2])){
  #
  #   ge <- genes[[i]]
  #   nn <- attr(ge, 'name')
  #   ge <- paste0(ge[c(T,F,F)],ge[c(F,T,F)],ge[c(F,F,T)])
  #   ch <- mmmut[[1]][which(mmmut[[2]]==nn),,drop=FALSE]
  #   rns <- as.integer(rownames(ch))
  #   gna <- colnames(pm)[i]
  #   og <- rep(list(ge), norg)
  #   names(og) <- paste0(phy$tip.label, '_', gna, ' ref:', nn)
  #
  #   for (j in seq_len(dim(ch)[1])){
  #     # dsc <- Descendants(phy, ch[j, 3L])[[1]]
  #     # dsc <- phangorn:::bip(phy)[ch[j, 3L]][[1]]
  #     dsc <- allDes[[.subset2(ch, j, 3)]]
  #     rnsj <- .subset2(rns, j)
  #     chj5 <- .subset2(ch, j, 5L)
  #     change.from[rnsj] <- og[[dsc[1]]][chj5]
  #     change.to[rnsj] <- sample(allcodons, 1, prob = .subset2(smat, change.from[rnsj]))
  #     for (k in seq_along(dsc)){
  #       og[[dsc[k]]][chj5] <- change.to[rnsj]
  #     }
  #   }
  #
  #   og <- lapply(og, paste0, collapse='')
  #
  #   aaa <- pm[, i]
  #   class(aaa) <- 'logical'
  #   og <- og[paste0(names(which(aaa)), '_', gna, ' ref:', nn)]
  #
  #   write.fasta(sequences = og,
  #               names = names(og),
  #               file.out = paste0(dir_out, '/', gna,'.fasta'))
  # }
  #
  # dfmut <- as.data.frame(mmmut[[1]])
  # dfmut$change.from <- change.from
  # dfmut$change.to <- change.to
  # dfmut$gene <- mmmut[[2]]
  #
  # ############
  # ## Return ##
  # ############
  #
  # # Return a list with a coalescent tree, a data.frame representing the
  # # gain-loss events, another one representing mutation events, and the final
  # # panmatrix. Also some attributes regarding input parameters are returned.
  #
  # out <- list(phy, dfgl, dfmut, pm)
  # names(out) <- c('coalescent', 'gain-loss', 'substitutions', 'panmatrix')
  # attr(out, 'reference') <- normalizePath(ref)
  # attr(out, 'norg') <- norg
  # attr(out, 'ngenes') <- ngenes
  # attr(out, 'ne') <- ne
  # attr(out, 'ggr') <- ggr
  # attr(out, 'glr') <- glr
  # attr(out, 'mu') <- mu
  # attr(out, 'dir_out') <- dir_out
  # attr(out, 'class') <- 'pangenomeSimulation'
  #
  # cat('DONE\n')
  # return(out)

}


