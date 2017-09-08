## ----libraries, message=FALSE--------------------------------------------
library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)
library(cytominergallery)
library(SNFtool)

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


## ----message=FALSE-------------------------------------------------------

doParallel::registerDoParallel(cores = 4)

feature_replicate_correlations <-
  profiles %>%
  cytominer::variable_importance(
    variables = variables,
    strata = c("Image_Metadata_Compound", "Image_Metadata_Concentration"),
    replicates = 3,
    cores = 2)

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


plot_tsne <- function(correlation, mechanism){
  set.seed(42)
  df <-
    tibble::as_data_frame(
      tsne::tsne(as.dist(1 - correlation))
      ) %>%
    mutate(mechanism = mechanism)

  p <-
    ggplot(df, aes(V1, V2, color = mechanism)) +
    geom_point() +
    ggtitle("t-SNE visualization of compound profiles clean using SNF")

  print(p)
}

plot_tsne(correlation = correlation, mechanism = mechanism)

## ----message=FALSE-------------------------------------------------------

predict_moa <- function(correlation, compound, mechanism){
  mask <- as.integer(outer(compound, compound, FUN = "!="))
  mask[mask == 0] <- -Inf

  correlation_masked <- correlation * mask

  return(sapply(1:nrow(correlation_masked),
                 function(i) mechanism[order(correlation_masked[i,],
                                             decreasing = TRUE)[1]])
    )
}

compound <- profiles$Image_Metadata_Compound
prediction <- predict_moa(correlation, compound, mechanism)
confusion_matrix <- caret::confusionMatrix(prediction, mechanism)

## ------------------------------------------------------------------------

evaluate_prediction <- function(confusion_matrix ){
  tibble::frame_data(
    ~metric, ~value,
    "Accuracy", sprintf("%.2f", confusion_matrix$overall["Accuracy"]),
    "95% CI", sprintf("(%.2f, %.2f)", confusion_matrix$overall[["AccuracyLower"]],
                      confusion_matrix$overall[["AccuracyUpper"]])
    )
}

evaluate_prediction(confusion_matrix) %>%
  knitr::kable(digits = 2)

## ------------------------------------------------------------------------
confusion_matrix$table %>%
  knitr::kable()


## ------------------------------------------------------------------------
distance <- 1-correlation
affinity <- affinityMatrix(distance,  K = 20, sigma = 0.5)
snf_distance <- SNF(list(affinity, affinity), K = 20, t = 20 )

## ------------------------------------------------------------------------
prediction <- predict_moa(snf_distance, compound, mechanism)
confusion_matrix <- caret::confusionMatrix(prediction, mechanism) 

evaluate_prediction(confusion_matrix) %>%
  knitr::kable(digits = 2)

## ----fig.width=8, fig.height=6, message=FALSE----------------------------

plot_tsne(correlation = snf_distance, mechanism = mechanism)


## ------------------------------------------------------------------------

feature_label = list("_AreaShape_", "_Intensity_","_Texture_","_Neighbors_")

feature_names_list <- lapply(feature_label, function(x) (
  profiles %>%
    select(matches(x)) %>% 
    colnames())
  )

correlation_matrices <- lapply(feature_names_list, function(x) (
  profiles %>%
    select(one_of(x)) %>%
    as.matrix() %>%
    t() %>%
    cor())
)

affinity_list <- lapply(correlation_matrices, function(x) (
  affinityMatrix(1 - x, K = 20, sigma = 0.5 ))
)

snf_distance <- SNF(affinity_list, K = 20, t = 20 )

## ------------------------------------------------------------------------

prediction <- predict_moa(snf_distance, compound, mechanism)

confusion_matrix <- caret::confusionMatrix(prediction, mechanism) 

evaluate_prediction(confusion_matrix) %>%
  knitr::kable(digits = 2)

## ----fig.width=8, fig.height=6, message=FALSE----------------------------

plot_tsne(correlation = snf_distance, mechanism = mechanism)


## ------------------------------------------------------------------------
#profiles %<>%
#  filter(Image_Metadata_Compound != "DMSO")

constituents = list("Cells", "Nuclei","Cytoplasm")

feature_list <- lapply(constituents, function(x)
  (profiles %>%
    select(matches(x)) %>% 
    colnames())
  )

correlation_matrices <- lapply(feature_list, function(x) (
  profiles %>%
  select(one_of(x)) %>%
  as.matrix() %>%
  t() %>%
  cor())
)

affinity_list <- lapply(correlation_matrices, function(x) (
  affinityMatrix(1 - x, K = 20, sigma = 0.5 )
))

snf_distance <- SNF(affinity_list, K = 20, t = 20 )

## ----fig.width=8, fig.height=6, message=FALSE----------------------------
prediction <- predict_moa(snf_distance, compound, mechanism)

confusion_matrix <- caret::confusionMatrix(prediction, mechanism) 

evaluate_prediction(confusion_matrix) %>%
  knitr::kable(digits = 2)

plot_tsne(correlation = snf_distance, mechanism = mechanism)


