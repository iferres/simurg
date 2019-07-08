
###############################################################################
#                             HELPER FUNCTIONS                                #
###############################################################################

## Simulate gene gain and loss
#' @importFrom phangorn Descendants
#' @importFrom stats runif rbinom setNames
#' @importFrom reshape2 melt
#' @author Ignacio Ferres
.sim_gl <- function(phy, m, brti, ne, norg, ou, ov, mrca_acc_o){

  # Generate random deviate integers respect the expected number of gains using
  # a poisson distribution, given a gain rate and a branch length, and
  # distribute them uniformly along branches.
  revbrti <-  max(brti) - brti
  dpt <- max(revbrti)
  gain <- apply(m, 1, function(x){

    ti <- (x[3] * ne) / dpt

    min <- revbrti[as.character(x[1])]
    max <- revbrti[as.character(x[2])]

    min.ti <- round((min * ne) / dpt)
    max.ti <- round((max * ne) / dpt)

    # Note for me:
    # unique() call should eventually be removed and below code
    # should cope with these marginal cases better.
    unique(round(runif(ti * ou, min=min.ti, max=max.ti)))

  })

  if (!length(gain)) gain <- rep(list(numeric(0)), dim(m)[1])

  mgain <- melt(setNames(gain, paste(m$V1, m$V2, sep = '-')))
  mgain <- mgain[order(mgain$value, decreasing = FALSE), ]

  # Segments. Initialize list where to save gene loss events.
  sgmt <- paste(m$V1, m$V2, sep  =  '-')
  loss <- setNames(vector('list', length(sgmt)), nm = sgmt)

  # Paths from root to tips
  roottotip <- setNames(ape::nodepath(phy), nm = seq_len(norg))

  # Empthy vector to fill acc genome size on each node
  node_ag_size <- vector('integer', length(brti))
  names(node_ag_size) <- names(brti)
  # MRCA node is the first one after the tips:
  node_ag_size[as.character(norg + 1)] <- mrca_acc_o


  for (i in seq_len(norg)){

    nds <- roottotip[[i]]
    # Iterate over branches from root to tip
    for (j in seq_len(length(nds) - 1L)){
      sg <- c(nds[j], nds[j+1])
      # Generations at which tree bifurcates (nodes)
      tnode <- round(revbrti[sg] * ne / dpt)
      psg <- paste(sg, collapse = '-')
      # Do not repeat already computed branches
      if (psg %in% sgmt){
        # Initialize variable to save generations at which gl event occur
        gloss <- c()
        # Get gene gains for this segment (branch)
        ggains <- mgain[which(mgain$L1==psg), 1, drop=TRUE]
        ti <- tnode[1]
        tfnode <- tnode[2]
        ag_size <- node_ag_size[as.character(nds[j])]
        for (tf in c(ggains, tfnode)){
          # Length of sub-segment
          tl <- tf - ti
          # P loosing existing genes
          pl <- ov * ag_size ##
          # Simulate gene loss events. Number of successes in 'tl' generations.
          # Here the algorithm don't mimic exactly what IMG model says:
          x <- rbinom(1, size = tl, prob = pl)
          # Next routine corrects behavoir to mimic IMG
          while(x>0){
            pos <- round(runif(x, min = ti, max = tf))
            pos <- min(pos)
            gloss <- c(gloss, pos)
            ag_size <- ag_size - 1L
            pl <- ov * ag_size
            ti <- pos
            x <- rbinom(1, size = tf - ti, prob = pl)
          }
          # Update initial time of sub-segment
          ti <- tf + 1L
          # Add a gene since every iteration is over computed gene gains.
          ag_size <- ag_size + 1L
        }
        loss[[psg]] <- if (length(gloss)>0) gloss else integer()
        # "Remove" segment (mark as already computed)
        sgmt <- sgmt[-which(sgmt==psg)]
        # Save accessory genome size of end node
        node_ag_size[as.character(nds[j+1])] <-  ag_size
      }
    }

  }

  if (!length(loss)) loss <- rep(list(numeric(0)), dim(m)[1])

  dfg <- lapply(1:(dim(m)[1]), function(x){

    if (length(gain[[x]])>0){
      cbind(gain[[x]], m[x, 1], m[x, 2], 'G')
    }else{
      NULL
    }

  })
  dfg <- do.call(rbind, dfg)

  dfl <- lapply(1:(dim(m)[1]), function(x){

    if (length(loss[[x]])>0){
      cbind(loss[[x]], m[x, 1], m[x, 2], 'L')
    }else{
      NULL
    }

  })
  dfl <- do.call(rbind, dfl)


  df <- do.call(rbind,  list(dfg, dfl))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  colnames(df) <- c('generation', 'from.node', 'to.node', 'type')
  df$generation <- as.numeric(df$generation)
  df$from.node <- as.integer(df$from.node)
  df$to.node <- as.integer(df$to.node)
  # df$type <- as.factor(df$type)
  df <- df[order(df$generation, decreasing = FALSE), ]

  pm <- matrix(data = 1L,
               nrow = norg,
               ncol = mrca_acc_o,
               dimnames = list(paste0('genome', 1:norg),
                               NULL))

  allDes <- Descendants(phy, type = 'tips')
  allorgs <- seq_len(norg)

  for (i in 1:dim(df)[1]){
    # If it is a gain, then add a new column with ones and zeros according the
    # descendants of the node where the gene was gained.
    if (.subset2(df, c(4, i)) == 'G'){

      # desc <- phy$tip.label[allDes[[.subset2(df, c(3, i))]]]
      desc <- .subset2(phy, 3)[allDes[[.subset2(df, c(3, i))]]]
      nwc <- rep(0, norg)
      nwc[as.integer(sub('genome','',desc))] <- 1L
      pm <- cbind(pm, nwc) # Future versions: modify in place with data.table
      # ndm <- dim(pm)
      # colnames(pm)[ndm[2]] <- paste0('gene', ndm[2])

      # If it is a loss, then put zeros at the corresponding column and rows.
    } else {

      # desc <- phy$tip.label[allDes[[.subset2(df, c(3, i))]]]
      desc <- .subset2(phy, 3)[allDes[[.subset2(df, c(3, i))]]]
      # Loose existing genes from node descendants
      ex <- which(colSums(pm[desc, , drop=FALSE])!=0) # Expensive
      # sample
      lo <- round(runif(1, min = 1, max = length(ex)))
      exlo <- ex[lo]
      pm[desc, exlo] <- 0L
      if (sum(.subset(pm, allorgs, exlo))==0)
        pm <- .subset(pm, allorgs, seq_len(dim(pm)[2])[-exlo])  #pm <- pm[, -exlo]

    }

  }

  # colnames(pm) <- paste0('gene', seq_len(dim(pm)[2]))



  return(pm)
}



#' @importFrom ape getMRCA nodepath
#' @importFrom stats setNames rpois
#' @importFrom phangorn Descendants
.sim_mut <- function(pm, phy, nsites, m, depth, brti, ne, norg, mu, smat){

  revbrti <-  max(brti) - brti
  dpt <- max(revbrti)

  if (missing(smat)){

    smat <- .codon.subst.mat

  }

  # Paths from root to tips
  roottotip <- setNames(nodepath(phy), nm = seq_len(norg))
  # Create structure to store result
  res_mut <- setNames(vector('list', dim(pm)[2]), names(nsites))


  # Iterate over genes
  for (g in seq_len(dim(pm)[2])){

    # ge <- genes[[i]]
    # nsites <- length(ge)
    mrca <- getMRCA(phy, names(which(pm[, g]==1L)))
    if (!is.null(mrca)){
      tips <- Descendants(phy, mrca, 'tips')[[1]]

      # Segments.
      ismrca <- which(sapply(roottotip, function(x) mrca %in%x))
      rttp <- lapply(roottotip[ismrca], function(x){
        x[which(x==mrca):length(x)]
      })
      sgmt <- m[which(m$V1 %in% unlist(rttp)), c(1,2)]
      sgmt <- paste(sgmt$V1, sgmt$V2, sep  =  '-')
    }else{
      tips <- unname(which(pm[, g]==1))
      # Segments
      sgmt <- paste(rev(rev(roottotip[[which(pm[, g]==1L)]])[1:2]), collapse = '-')
    }


    # Create structure to save simulated mutations
    muts <- setNames(vector('list', length = length(sgmt)), nm = sgmt)

    #Iterate over descendant organisms
    for (i in tips){
      nds <- roottotip[[i]]
      # Iterate over branches
      for (j in seq_len(length(nds) - 1L)){
        sg <- c(nds[j], nds[j+1])
        psg <- paste(sg, collapse = '-')
        # Do not repeat if already computed branch
        if (psg %in% sgmt){
          # Simulate number of substitutions for this segment (binomial)
          slth <- revbrti[sg[2]] - revbrti[sg[1]]
          ngen <- unname( slth * ne / dpt)
          trials <- ngen * nsites[g]

          # # Too big, coerces it into numeric for double presision :,(
          # trials <- ngen * nsites
          # # This produces NAs because expects class(trials) == 'integer' :
          # rbinom(1, size = trials, prob = mu)
          # # Seen on stackoverflow but also produces NaNs:
          # qbinom(runif(1), trials, mu, FALSE)

          # Approximation with poisson  ¯\_(ツ)_/¯ :
          nmut <- rpois(1, trials * mu)

          if (nmut){
            # Distribute substitutions with uniform probability
            smut <- sample(seq_len(nsites[g]), nmut, replace = TRUE)
            muts[[psg]] <- smut
          }else{
            muts[[psg]] <- NA_integer_
          }

          # "Remove" segment (mark as already computed)
          sgmt <- sgmt[-which(sgmt==psg)]
        }

      }

    }

    #print(g)
    res_mut[[g]] <- muts

    #Finish gene[g]

  }


  res_mut
}
