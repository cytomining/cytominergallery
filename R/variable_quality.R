#' Measure variable quality based on replicate correlation
#'
#' @param sample ...
#' @param variables ...
#' @param strata ...
#' @param replicates ...
#' @param key ...
#' @param batch ...
#'
#' @return data.frame of variable quality measurements
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom foreach %dopar%
#' @importFrom stats median
#'
#' @export
#'
variable_quality <-
  function(sample, variables, strata, replicates, key = NULL, batch = NULL) {

    replicate_correlation <- function(sample, variable, strata, key) {

      correlation_matrix <-
        sample %>%
        dplyr::arrange_(.dots = strata) %>%
        dplyr::select_(.dots = c(strata, variable)) %>%
        tidyr::spread_(key, variable) %>%
        dplyr::select_(~-dplyr::one_of(setdiff(strata, key))) %>%
        stats::cor()

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
        dplyr::mutate_(.dots = stats::setNames(list(~dplyr::row_number(batch)), key)) %>%
        dplyr::ungroup()

      strata <- c(strata, key)
    }


    foreach::foreach(variable = variables, .combine = rbind) %dopar%
    {

      sample %>%
        split(.[batch]) %>%
        purrr::map_df(replicate_correlation,
                      variable = variable,
                      strata = strata,
                      key = key) %>%
        dplyr::mutate(variable = variable)
    } %>%
      tidyr::gather_(key, "pearson", setdiff(names(.), "variable")) %>%
      dplyr::group_by_(.dots = c("variable")) %>%
      dplyr::summarize_at("pearson", dplyr::funs(median, min, max))
  }
