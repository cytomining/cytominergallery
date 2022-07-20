[![Travis-CI Build Status](https://travis-ci.org/cytomining/cytominergallery.svg?branch=master)](https://travis-ci.org/cytomining/cytominergallery)

# cytominer gallery

## Installation

Install [R](https://www.r-project.org) and [RStudio](https://www.rstudio.com/).

Download this [sqlite](https://cellpainting-gallery.s3.amazonaws.com/cpg0010-caie-drugresponse/broad-az/workspace/backend/ljosa_2013/analysis_batch_stable/ljosa_jbiomolscreen_2013.sqlite) file into `~/Downloads` (required by the vignette `single_cell_analysis`).

Install the `cytominergallery` package from GitHub:

```R
# install.packages("devtools")
devtools::install_github("cytomining/cytominergallery", dependencies = TRUE, build_vignettes = TRUE)
```

You may need to do run that again in order to build the vignettes correctly (seems like a bug in `install_github`):
```R
devtools::install_github("cytomining/cytominergallery", dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
```

Occasionally, the `Suggests` dependencies [may not get installed](https://github.com/hadley/devtools/issues/1370), depending on your system, so you'd need to install those explicitly.

Browse vignettes (launches in default browser):
```R
browseVignettes()
```

Search for "Vignettes in package cytominergallery" and click on the link `HTML` to view the vignettes.


