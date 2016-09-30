[![Travis-CI Build Status](https://travis-ci.org/shntnu/cytominerworkshop.svg?branch=master)](https://travis-ci.org/shntnu/cytominerworkshop)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/shntnu/cytominerworkshop?branch=master&svg=true)](https://ci.appveyor.com/project/shntnu/cytominerworkshop)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/cytominerworkshop)](https://cran.r-project.org/package=cytominerworkshop)

# cytominer workshop

## Installation

Install [R](https://www.r-project.org) and [RStudio](https://www.rstudio.com/).

Download this [sqlite](https://s3.amazonaws.com/imaging-platform-collaborator/2016_09_09_cytominer_workshop/ljosa_jbiomolscreen_2013.sqlite) file into `~/Downloads` (required by the vignette `single_cell_analysis`).

Install the `cytominerworkshop` package from GitHub:

```R
# install.packages("devtools")
devtools::install_github("shntnu/cytominerworkshop", dependencies = TRUE, build_vignettes = TRUE)
```

You may need to do run that again in order to build the vignettes correctly (seems like a bug in `install_github`):
```R
devtools::install_github("shntnu/cytominerworkshop", dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
```

Browse vignettes (launches in default browser):
```R
browseVignettes()
```

Search for "Vignettes in package cytominerworkshop" and click on the link `HTML` to view the vignettes.


