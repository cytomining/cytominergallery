#' Measure feature quality based on replicate correlation
#'
#' @param sample ...
#' @param variables ...
#' @param strata ...
#' @param replicates ...
#' @param key ...
#' @param batch ...
#'
#' @return data.frame of feature quality measurements
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
feature_quality <-
  function(sample, variables, strata, replicates, key = NULL, batch = NULL) {

    replicate_correlation <- function(data, variable, strata, key) {

      correlation_matrix <-
        data %>%
        dplyr::arrange_(.dots = strata) %>%
        dplyr::select_(.dots = c(strata, variable)) %>%
        tidyr::spread_(key, variable) %>%
        dplyr::select_(~-dplyr::one_of(setdiff(strata, key))) %>%
        cor()

      median(correlation_matrix[upper.tri(correlation_matrix)])
    }

    if (is.null(batch)) {
      sample %<>% dplyr::mutate(batch = 0)

      batch <- "batch"
    }

    if (is.null(key)) {
      key <- "replicate_id"

      sample %<>%
        dplyr::count_(vars = strata) %>%
        dplyr::filter(n == replicates) %>%
        dplyr::select(-n) %>%
        dplyr::inner_join(sample) %>%
        dplyr::group_by_(.dots = strata) %>%
        dplyr::mutate_(.dots = setNames(list(~dplyr::row_number(batch)), key)) %>%
        dplyr::ungroup()

      strata <- c(strata, key)
    }

    foreach::foreach(feature_col = variables, .combine = rbind) %dopar%
    {

      data %>%
        split(.[batch]) %>%
        purrr::map_df(replicate_correlation,
                      feature_col = feature_col,
                      grouping_cols = strata,
                      key_col = key) %>%
        dplyr::mutate(feature_col = feature_col)
    } %>%
      tidyr::gather_(key, "pearson", setdiff(names(.), "feature_col")) %>%
      dplyr::group_by_(.dots = c("feature_col")) %>%
      dplyr::summarize_at("pearson", dplyr::funs(median, min, max))
  }
