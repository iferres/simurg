
###############################################################################
#                             HELPER FUNCTIONS                                #
###############################################################################

## Simulate gene gain and loss
#' @importFrom phangorn Descendants
#' @importFrom stats runif rpois
#' @author Ignacio Ferres
.sim_gl <- function(phy, m, depth, brti, ne, norg, ngenes, theta, rho){


  # theta <- 2e-6
  # Generate random deviate integers respect the expected number of gains using
  # a poisson distribution, given a gain rate and a branch length, and
  # distribute them uniformly along branches.
  gain <- apply(m, 1, function(x){

    ti <- (x[3] * ne) / depth

    min <- brti[as.character(x[2])]
    max <- brti[as.character(x[1])]

    min.ti <- round((min * ne) / depth)
    max.ti <- round((max * ne) / depth)

    round(runif(rpois(1, ti * theta), min=min.ti, max=max.ti))

  })

  # rho <- 1e-6
  # Anologous as above, but with loss rate.
  loss <- apply(m, 1, function(x){

    ti <- (x[3] * ne) / depth

    min <- brti[as.character(x[2])]
    max <- brti[as.character(x[1])]

    min.ti <- round((min * ne) / depth)
    max.ti <- round((max * ne) / depth)

    round(runif(rpois(1, ti * rho), min=min.ti, max=max.ti))

  })

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
  df$generation <- as.integer(df$generation)
  df$from.node <- as.integer(df$from.node)
  df$to.node <- as.integer(df$to.node)
  df$type <- as.factor(df$type)
  df <- df[order(df$generation, decreasing = TRUE), ]

  pm <- matrix(data = 1L,
               nrow = norg,
               ncol = ngenes,
               dimnames = list(paste0('genome', 1:norg),
                               paste0('gene', 1:ngenes)))

  for (i in 1:dim(df)[1]){
    # If it is a gain, then add a new column with ones and zeros according the
    # descendants of the node where the gene was gained.
    if (df$type[i] == 'G'){

      desc <- phy$tip.label[Descendants(phy, df$to.node[i])[[1]]]
      nwc <- rep(0, norg)
      nwc[as.integer(sub('genome','',desc))] <- 1L
      pm <- cbind(pm, nwc)
      ndm <- dim(pm)
      colnames(pm)[ndm[2]] <- paste0('gene', ndm[2])

      # If it is a loss, then put zeros at the corresponding column and rows.
    } else {

      lo <- sample(colnames(pm), 1)
      desc <- phy$tip.label[Descendants(phy, df$to.node[i])[[1]]]
      pm[as.integer(sub('genome','',desc)), lo] <- 0L

    }

  }

  # Remove empthy columns and re-write gene names
  if (any(colSums(pm)==0L)){
    pm <- pm[, -which(colSums(pm)==0L)]
    colnames(pm) <- paste0('gene', seq_len(dim(pm)[2]))
  }

  return(list(df, pm))
}








## Simulate mutations
#' @importFrom phangorn Ancestors
#' @importFrom stats runif rpois
#' @author Ignacio Ferres
.sim_mut <- function(phy, genes, m, depth, brti, ne, norg, ngenes, mu, smat){

  if (missing(smat)){

    smat <- .codon.subst.mat

  }

  # Count total number of sites.
  nsites <- sum(sapply(genes, length))

  # At each branch, simulate generations at which substitutions happens.
  # Generate random deviate integers respect the expected number of
  # substitutions using a poisson distribution, given a substitution rate, a
  # number of sites and a branch length,and distribute them uniformly along
  # branches.
  mutations <- apply(m, 1, function(x){

    ti <- (x[3] * ne) / depth

    min <- brti[as.character(x[2])]
    max <- brti[as.character(x[1])]

    min.ti <- round((min * ne) / depth)
    max.ti <- round((max * ne) / depth)

    round(runif(rpois(1, ti * mu * nsites), min=min.ti, max=max.ti))

  })

  dfmut <- lapply(1:(dim(m)[1]), function(x){

    if (length(mutations[[x]])>0){
      cbind(mutations[[x]], m[x, 1], m[x, 2])
    }else{
      NULL
    }

  })
  dfmut <- as.data.frame(do.call(rbind, dfmut), stringsAsFactors = FALSE)
  colnames(dfmut) <- c('generation', 'from.node', 'to.node')
  dfmut$generation <- as.integer(dfmut$generation)
  dfmut$from.node <- as.integer(dfmut$from.node)
  dfmut$to.node <- as.integer(dfmut$to.node)
  dfmut <- dfmut[order(dfmut$generation, decreasing = TRUE),]
  dfmut$position <- as.integer(round(runif(dim(dfmut)[1], min = 1, max = nsites)))
  dfmut$codon.pos <- ceiling(dfmut$position/3)

  # codones <- strsplit(paste0(genes, collapse = ''), "(?<=.{3})", perl = TRUE)[[1]]
  codones <- lapply(genes, function(x){
    paste0(x[c(T,F,F)], x[c(F,T,F)], x[c(F,F,T)])
  })
  codones <- unlist(codones)
  dfmut$change.from <- codones[dfmut$codon.pos]


  # Mutate positions according Jukes-Cantor model (??? sort of...).
  ncodon <- colnames(smat)
  dfmut$change.to <- sapply(dfmut$change.from,
                            function(x, ncodon, smat) {
                              sample(ncodon, 1, prob = smat[,x])
                            },
                            ncodon = ncodon,
                            smat = smat)


  # Correct multiples changes on same position if first change affect the branch
  # of the second.
  dups <- duplicated(dfmut$codon.pos) | duplicated(dfmut$codon.pos, fromLast = TRUE)
  un <- unique(dfmut$codon.pos[dups])
  for (i in seq_along(un)){
    aa <- which(dfmut$codon.pos==un[i])
    for (j in 2:length(aa)){
      K <- j - 1
      while(K>0){
        ab <- dfmut[aa[K], ]
        a1 <- c(Ancestors(phy, ab$to.node), ab$to.node)
        ac <- dfmut[aa[j], ]
        a2 <- c(Ancestors(phy, ac$to.node), ac$to.node)
        #branch and position are ancestral?
        if (all(a1 %in% a2)){
          dfmut$change.from[aa[j]] <- dfmut$change.to[aa[K]]
          dfmut$change.to[aa[j]] <- sample(ncodon, 1, prob = smat[,dfmut$change.from[aa[j]]])
          K <- 0
        }else{
          K <- ifelse(K==1, 0, K-1)
        }
      }
    }
  }

  #Identify gene from absolute position. Identify relative position (in gene).
  gle <- sapply(genes, length)
  csu <- ceiling(cumsum(gle)/3)
  dfmut$gene <- sapply(dfmut$codon.pos, function(x){names(which(csu>=x))[1]})
  dfmut$gene.codon.pos <- sapply(dfmut$codon.pos, function(x){
    ifelse(length(which(csu<=x)), x-csu[rev(which(csu<x))[1]], x)
  })

  dfmut$position <- NULL

  return(dfmut)
}






