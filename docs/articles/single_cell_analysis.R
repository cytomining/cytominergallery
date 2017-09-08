## ----libraries, message=FALSE--------------------------------------------
library(dplyr)
library(ggplot2)
library(magrittr)


## ----load, message=FALSE-------------------------------------------------

backend <- file.path(Sys.getenv("HOME"), "Downloads", "ljosa_jbiomolscreen_2013.sqlite")

db <- src_sqlite(path = backend)

image <- tbl(src = db, "Image")

object <-
  tbl(src = db, "Cells") %>%
  inner_join(tbl(src = db, "Cytoplasm"),
                    by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
  inner_join(tbl(src = db, "Nuclei"),
                    by = c("TableNumber", "ImageNumber", "ObjectNumber"))

object %<>% inner_join(image, by = c("TableNumber", "ImageNumber"))


## ------------------------------------------------------------------------
variables <-
  colnames(object) %>%
  stringr::str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

## ------------------------------------------------------------------------
print(length(variables))

## ------------------------------------------------------------------------
object %>%
  count() %>%
  knitr::kable(caption = "No. of cells")

## ----fig.width=5, fig.height=5-------------------------------------------
object %>%
  filter(Image_Metadata_Plate == "Week1_22123" && Image_Metadata_Well == "E04") %>%
  select(Nuclei_Intensity_IntegratedIntensity_CorrDAPI) %>%
  collect() %>% {
    ggplot(., aes(Nuclei_Intensity_IntegratedIntensity_CorrDAPI)) +
      geom_histogram(binwidth = 5)
    }


## ----message=FALSE-------------------------------------------------------
profiles <-
  cytominer::aggregate(
    population = object,
    variables = variables,
    strata = c("Image_Metadata_Plate", "Image_Metadata_Well"),
    operation = "mean"
  )

profiles %<>%
  collect()


## ------------------------------------------------------------------------
profiles %>%
  count() %>%
  knitr::kable(caption = "No. of wells")

## ----fig.width=6, fig.height=5, message=FALSE----------------------------
p <-
  ggplot(profiles, aes(Cells_AreaShape_Area, Nuclei_AreaShape_Area)) +
  geom_hex()

print(p)

