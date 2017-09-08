## ----libraries, message=FALSE--------------------------------------------
library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(cytominergallery)

## ----message=FALSE-------------------------------------------------------
profiles <-
  readr::read_csv(system.file("extdata", "ljosa_jbiomolscreen_2013_per_well_mean.csv",
                package = "cytominergallery"))

moa <-
  readr::read_csv(system.file("extdata", "BBBC021_v1_moa.csv",
                              package = "cytominergallery")) %>%
  rename(Image_Metadata_Compound = compound,
                Image_Metadata_Concentration = concentration,
                Image_Metadata_MoA = moa
  )

metadata <-
  readr::read_csv(system.file("extdata", "BBBC021_v1_image.csv",
                              package = "cytominergallery")) %>%
  rename(Image_Metadata_Plate = Image_Metadata_Plate_DAPI,
                Image_Metadata_Well = Image_Metadata_Well_DAPI
  ) %>%
  select(matches("^Image_Metadata")) %>%
  inner_join(moa) %>%
  distinct()

profiles %<>%
  inner_join(metadata)

variables <-
  colnames(profiles) %>%
  str_subset("^Nuclei_|^Cells_|^Cytoplasm_")


## ------------------------------------------------------------------------
profiles %>%
  filter(Image_Metadata_Compound != "DMSO") %>%
  distinct(Image_Metadata_Compound) %>%
  tally() %>%
  rename(`Number of compounds` = n) %>%
  knitr::kable()

## ------------------------------------------------------------------------
profiles %>%
  filter(Image_Metadata_Compound != "DMSO") %>%
  distinct(Image_Metadata_Compound, Image_Metadata_Concentration) %>%
  tally() %>%
  rename(`Number of unique treatments` = n) %>%
  knitr::kable()

## ------------------------------------------------------------------------
profiles %>%
  filter(Image_Metadata_Compound != "DMSO") %>%
  count(Image_Metadata_Compound, Image_Metadata_Concentration) %>%
  rename(`Number of replicates` = n) %>%
  knitr::kable()

## ------------------------------------------------------------------------
profiles %>%
  filter(Image_Metadata_Compound == "DMSO") %>%
  count(Image_Metadata_Plate) %>%
  rename(`Number of DMSO wells` = n) %>%
  knitr::kable()

## ----message=FALSE-------------------------------------------------------
profiles <-
  cytominer::select(
    population = profiles,
    variables = variables,
    sample = profiles,
    operation = "variance_threshold"
  ) %>%
  collect()

variables <-
  colnames(profiles) %>%
  str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

## ----message=FALSE-------------------------------------------------------
doParallel::registerDoParallel(cores = 2)

feature_replicate_correlations <-
  profiles %>%
  cytominer::variable_importance(
    variables = variables,
    strata = c("Image_Metadata_Compound", "Image_Metadata_Concentration"),
    replicates = 3,
    cores = 2)

## ----fig.width=4, fig.height=4, message=FALSE----------------------------
ggplot(feature_replicate_correlations, aes(median))  +
  stat_ecdf() +
  geom_vline(xintercept = 0.5, color = "red") +
  xlab("median replicate correlation (Pearson)") +
  ylab("F(x)")

## ----message=FALSE-------------------------------------------------------

profiles %<>%
  select_(.dots = setdiff(x = colnames(profiles),
                          y = feature_replicate_correlations %>%
                            filter(median < 0.5) %>%
                            magrittr::extract2("variable"))
          )

variables <-
  colnames(profiles) %>%
  str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

## ----message=FALSE-------------------------------------------------------
profiles <-
  cytominer::select(
    population = profiles,
    variables = variables,
    sample = profiles,
    operation = "correlation_threshold",
    cutoff = 0.95) %>%
  collect()

variables <-
  colnames(profiles) %>%
  str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

## ----message=FALSE-------------------------------------------------------

profiles <-
  cytominer::normalize(
    population = profiles,
    variables = variables,
    strata =  c("Image_Metadata_Plate"),
    sample = profiles %>% filter(Image_Metadata_Compound == "DMSO")
  )

profiles <-
  cytominer::select(
      population = profiles,
      variables = variables,
      operation = "drop_na_columns"
  )

variables <-
  colnames(profiles) %>%
  str_subset("^Nuclei_|^Cells_|^Cytoplasm_")


## ----message=FALSE-------------------------------------------------------

profiles <-
  cytominer::aggregate(
    population = profiles,
    variables = variables,
    strata = c("Image_Metadata_Compound",
               "Image_Metadata_Concentration",
               "Image_Metadata_MoA"),
    operation = "mean"
  )

## ----fig.width=8, fig.height=6, message=FALSE----------------------------

profiles %<>%
  filter(Image_Metadata_Compound != "DMSO")

correlation <-
  profiles %>%
  select(one_of(variables)) %>%
  as.matrix() %>%
  t() %>%
  cor()

mechanism <- as.character(profiles$Image_Metadata_MoA)

set.seed(123)

df <-
  tibble::as_data_frame(
    tsne::tsne(as.dist(1-correlation))
    ) %>%
  mutate(mechanism = mechanism)

p <-
  ggplot(df, aes(V1, V2, color=mechanism)) +
  geom_point() +
  ggtitle("t-SNE visualization of compound profiles")

print(p)


## ----message=FALSE-------------------------------------------------------

compound <- profiles$Image_Metadata_Compound

mask <- as.integer(outer(compound, compound, FUN="!="))

mask[mask == 0] <- -Inf

correlation_masked <- correlation * mask

prediction <- sapply(1:nrow(correlation_masked),
               function(i) mechanism[order(correlation_masked[i,],
                                           decreasing = TRUE)[1]])

confusion_matrix <- caret::confusionMatrix(prediction, mechanism)

## ------------------------------------------------------------------------
tibble::frame_data(
  ~metric, ~value,
  "Accuracy", sprintf("%.2f", confusion_matrix$overall["Accuracy"]),
  "95% CI", sprintf("(%.2f, %.2f)", confusion_matrix$overall[["AccuracyLower"]],
                    confusion_matrix$overall[["AccuracyUpper"]])
  ) %>%
  knitr::kable(digits = 2)

## ------------------------------------------------------------------------
confusion_matrix$table %>%
  knitr::kable()


