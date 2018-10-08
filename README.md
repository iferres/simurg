pansimulatoR
================
I Ferrés
October 8, 2018

Simulate a bacterial pangenome in R
-----------------------------------

This R package is intended to produce simulated pangenomes using reference sequences as starting point (MRCA), and both Neutral \[1\] and Infinitely Many Genes (IMG) \[2\] models to produce changes along branches of a simulated coalescent tree.

Implementation
--------------

The algorithm first simulates a random coalescent tree.

Then, random deviates from the expected number of gene gain and loss according to branch length and other parameters (*theta*: gene gain rate per generation; and *rho*; gene loss rate per generation), are used to simulate gene birth and death along the tree. Since the IMG model is used, genes are only transmited vertically from generation to generation.

Point mutations are simulated following a similar process as the above: random deviates from the expected number of mutations according to branch length and mutation rate (*mu*), are used. Mutations are then distributed along sites with uniform probability. By default stop codons are avoided. (**note:** actually, codons and not sites are conditioned to mutation).

Finally, sequences are sampled from a reference multi-fasta file, and a pangenome is generated following the evolutionary history simulated previously.

Example
-------

First load package and decompress the attached reference sequences (for tutorial purposes only).

``` r
setwd(tempdir())

#Load library
library(pansimulatoR)
```

    ## Loading required package: ape

``` r
#List attached file
tgz <- system.file('extdata', 'ref_tutorial.tar.gz', package = 'pansimulatoR')
tgz
```

    ## [1] "/home/ignacio/R/x86_64-redhat-linux-gnu-library/3.3/pansimulatoR/extdata/ref_tutorial.tar.gz"

``` r
ref <- untar(tarfile = tgz, list = T)

#Decompress
untar(tarfile = tgz,exdir = tempdir())
ref <- list.files(path = tempdir(), 
                  pattern = ref, 
                  full.names = TRUE)
ref
```

    ## [1] "/tmp/RtmpXz0hvd/ref_tutorial.fasta"

Then run the simulation using the extracted reference sequences, and default parameters.

``` r
set.seed(123)
p <- simpg(ref = ref)
```

    ## Discarded 13 sequences due to presence of "n".
    ## Simulating coalescent tree.
    ## Simulating gene gain and loss.
    ## Simulating point mutations.
    ## Writing groups of orthologous at "sim_pg" .
    ## DONE

``` r
#length
length(p)
```

    ## [1] 4

``` r
#names
names(p)
```

    ## [1] "coalescent"    "gain-loss"     "substitutions" "panmatrix"

``` r
#class
lapply(p, class)
```

    ## $coalescent
    ## [1] "phylo"
    ## 
    ## $`gain-loss`
    ## [1] "data.frame"
    ## 
    ## $substitutions
    ## [1] "data.frame"
    ## 
    ## $panmatrix
    ## [1] "matrix"

``` r
#attributes
attributes(p)
```

    ## $names
    ## [1] "coalescent"    "gain-loss"     "substitutions" "panmatrix"    
    ## 
    ## $reference
    ## [1] "/tmp/RtmpXz0hvd/ref_tutorial.fasta"
    ## 
    ## $norg
    ## [1] 10
    ## 
    ## $ngenes
    ## [1] 100
    ## 
    ## $ne
    ## [1] 5e+07
    ## 
    ## $theta
    ## [1] 2e-06
    ## 
    ## $rho
    ## [1] 1e-06
    ## 
    ## $mu
    ## [1] 5e-10
    ## 
    ## $dir_out
    ## [1] "sim_pg"
    ## 
    ## $class
    ## [1] "pangenomeSimulation"

As seen, a `list` of 4 elements is returned. The first element is the simulated coalescent tree (class `phylo`, from `ape` package).

``` r
p$coalescent
```

    ## 
    ## Phylogenetic tree with 10 tips and 9 internal nodes.
    ## 
    ## Tip labels:
    ##  genome8, genome1, genome4, genome6, genome2, genome7, ...
    ## 
    ## Rooted; includes branch lengths.

``` r
plot(p$coalescent)
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

The second element is a `data.frame` with the simulated history of gene gain (G) and loss (L) events along the tree.

``` r
head(p$`gain-loss`)
```

    ##   generation from.node to.node type
    ## 1   49972449        11      12    L
    ## 2   49323837        11      13    G
    ## 3   49303995        11      12    G
    ## 4   49269860        11      12    G
    ## 5   49173677        11      12    L
    ## 6   48921604        11      13    L

The third element is a `data.frame` with the simulated history of mutation events along the tree.

``` r
head(p$substitutions)
```

    ##   generation from.node to.node codon.pos change.from change.to
    ## 1   49997287        11      12      1827         gga       gca
    ## 2   49995229        11      13     37915         caa       gaa
    ## 3   49990772        11      13     17194         aac       tac
    ## 4   49987055        11      13     62287         aat       act
    ## 5   49986604        11      13     10304         gac       gcc
    ## 6   49983853        11      12     45083         tgc       ggc
    ##              gene gene.codon.pos
    ## 1  32019.11_01806            271
    ## 2 1093099.4_01806             40
    ## 3  593452.3_00206            113
    ## 4 1244529.3_00647              5
    ## 5 1093099.4_00853            392
    ## 6   32020.6_00398             80

And the forth element is the final panmatrix, which represent the gene presence/absence of a gene on each sampled organism.

``` r
#Dimensions
dim(p$panmatrix)
```

    ## [1]  10 303

``` r
#First 5 rows and columns
p$panmatrix[1:5, 1:5]
```

    ##         gene1 gene2 gene3 gene4 gene5
    ## genome1     0     1     1     0     1
    ## genome2     0     1     1     0     1
    ## genome3     1     1     1     1     1
    ## genome4     0     1     1     0     1
    ## genome5     1     1     1     1     1

The last relevant result is a directory where simulated final sequences are saved as fasta files.

``` r
dir_out <- attributes(p)$dir_out
normalizePath(dir_out)
```

    ## [1] "/tmp/RtmpOk92QC/sim_pg"

``` r
list.files(path = normalizePath(dir_out), full.names = TRUE)[1:5]
```

    ## [1] "/tmp/RtmpOk92QC/sim_pg/gene1.fasta"  
    ## [2] "/tmp/RtmpOk92QC/sim_pg/gene10.fasta" 
    ## [3] "/tmp/RtmpOk92QC/sim_pg/gene100.fasta"
    ## [4] "/tmp/RtmpOk92QC/sim_pg/gene101.fasta"
    ## [5] "/tmp/RtmpOk92QC/sim_pg/gene102.fasta"

Some thoughts on parameter selection
------------------------------------

-   **norg** is the number of sampled organisms. In other words, is the number of final organisms your pangenome will have.

-   **ngenes** is the number of genes your MRCA will have. It's the starting point of the simulation. The final number is undetermined because it depends on *theta* and *rho* ratio, and on each simulation *per se*.

-   **ne** is the effective population size of the species, although it is also interpreted as the number of generations from root to tips, see \[2\]. According to \[3\], for *E. coli* is ~ 25 million. The default is two times this number as some studies have pointed out that \[3\] may be too conservative.

-   **theta** and **rho** are the gene gain/loss rate per generation, respectively. I'm considering changing these parameter's names because they are the equivalent to *u* and *v* in \[2\]; `theta` = 2 x `Ne` *u*, and `rho` = 2 x `Ne` *v* in the IMG model. So, for now beware of that technical detail. In this package you are setting the probability of gene gain and loss **per generation**, independently from the `Ne`. The values used by default were arbitrarily set, so the probability of gene gain is two times the probability of gene loss.

-   **mu** is the per site per generation substitution rate. According to \[4\], a typical bacterial substitution rate is around 1e-5 and 1e-9 per site per annum. Taking, e.g., 1e-7 substitutions per site per annum, and assuming 200 generations per annum \[5\]: 1e-7 / 200 = 5e-10 substitutions per site per generation, which was taken as default.

References
----------

\[1\] Kimura, M. (1983). **The neutral theory of molecular evolution.** *Cambridge.*

\[2\] Baumdicker, F., Hess, W. R., & Pfaffelhuber, P. (2012). **The Infinitely Many Genes Model for the Distributed Genome of Bacteria.** *Genome Biology and Evolution.*

\[3\] Charlesworth J, Eyre-Walker A. (2006). **The rate of adaptive evolution in enteric bacteria.** *Mol Biol Evol.*

\[4\] Duchêne, S., et al. (2016). **Genome-scale rates of evolutionary change in bacteria.** *Microbial Genomics.*

\[5\] Ashley Robinson, D., Falush, D., Fiel, E. J. (2010). **Bacterial Population Genetics in Infectious Disease.** *Wiley-Blackwell*, Chapter 1: The Coalescent of Bacterial Populations; Section 1.5: From Coalescent Time To Real Time. (Schierup, M. H., and Wiuf, C.)
