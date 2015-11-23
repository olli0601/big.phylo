big.phylo
=========

An R package to construct phylogenies for very large numbers of sequences

# Installation

The easiest way to install `big.phylo` is via the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("olli0601/big.phylo")
```

# Help

Help files and documentation are available once the package is loaded in R with 'require(big.phylo)'. 
Documentation is available via 'library(help="big.phylo")'. 
The package contains:

* functions to create ExaML bootstrap phylogenies. Try ?pipeline.ExaML.bootstrap.per.proc

* an Rscript and functions to remove drug resistance mutations. Try ?seq.rm.drugresistance. This function is also accessible from the command line via
```
Rscript rm.drm.Rscript -indir=XXX -infile=XXY  -outdir=XYY -outfile=YYY
```
where `rm.drm.Rscript` is found in the R package installation directory.
