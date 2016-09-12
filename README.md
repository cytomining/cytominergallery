[![Travis-CI Build Status](https://travis-ci.org/shntnu/cytominrworkshop.svg?branch=master)](https://travis-ci.org/shntnu/cytominrworkshop)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/shntnu/cytominrworkshop?branch=master&svg=true)](https://ci.appveyor.com/project/shntnu/cytominrworkshop)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/cytominrworkshop)](https://cran.r-project.org/package=cytominrworkshop)

# cytominr workshop

## Installation

- Install [R](https://www.r-project.org)
- Install [RStudio](https://www.rstudio.com/)
- The vignette `single_cell_analysis` requires  <https://s3.amazonaws.com/imaging-platform-collaborator/2016_09_09_cytominr_workshop/ljosa_jbiomolscreen_2013.sqlite> to be downloaded into `~/Downloads`.

Now install the `cytominrworkshop` package from GitHub:

```R
# install.packages("devtools")
devtools::install_github("shntnu/cytominrworkshop", dependencies = TRUE, build_vignettes = TRUE)
```

You may need to do run that again in order to build the vignettes correctly (seems like a bug in `install_github`):
```R
devtools::install_github("shntnu/cytominrworkshop", dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
```

Browse vignettes (launches in default browser)
```R
browseVignettes()
```

Search for "Vignettes in package cytominrworkshop" and click on the link `HTML` to view the vignettes.


