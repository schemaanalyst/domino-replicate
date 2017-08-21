#' FUNCTION: analysis_means
#'
#' Averages of test generation timing and coverage
#' @param data A data frame of read_analysis
#' @return A Data Frame of averages
#' @importFrom magrittr %>%
#' @export
analysis_means <- function(data) {
  data %>% dplyr::select(dbms, casestudy, criterion, datagenerator, coverage, testgenerationtime) %>%
    dplyr::group_by(dbms, casestudy, datagenerator) %>%
    dplyr::summarise(testgenerationtiming = mean(testgenerationtime),
                     coverage = mean(coverage))

  return(data)
}

#' FUNCTION: analysis_medians
#'
#' Medians of test generation timing and coverage
#' @param data A data frame of read_analysis
#' @return A Data Frame of medians
#' @importFrom magrittr %>%
#' @export
analysis_medians <- function(data) {
  data %>% dplyr::select(dbms, casestudy, criterion, datagenerator, coverage, testgenerationtime) %>%
    dplyr::group_by(dbms, casestudy, datagenerator) %>%
    dplyr::summarise(testgenerationtiming = median(testgenerationtime),
                     coverage = median(coverage))

  return(data)
}
