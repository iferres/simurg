# simba: Simulate a Bacterial Pangenome in R

This R package is intended to produce simulated pangenomes using reference sequences as starting point (MRCA), and both Neutral and Infinitely Many Genes (IMG) models to produce changes along branches of a simulated coalescent tree.

## Installation

The simplest way of installing simba is using `devtools` package:

```r
if (!require(devtools)) {
    install.packages('devtools')
    library(devtools)
}
 
install_github('iferres/simba')
```

## Help

The main function of `simba` is well documented. Once the package is loaded, please refer to:
```r
help('simpg')
```

## Vignettes

A vignette is also provided to learn more about this package. If you want to use it you have to build it first when installing `simba`:


```r
library(devtools)
install_github('iferres/simba', build_vignettes = TRUE)
vignette('Simulate_a_bacterial_pangenome')
```


