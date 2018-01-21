#' FUNCTION: get_check_schemas_analysis
#'
#' Filtering analysis data to get ony check constraints schemas.
#' @param d Data frame of analysis
#' @return A data frame of check constraints only schemas
#' @importFrom magrittr %>%
#' @export
get_check_schemas_analysis <- function(d) {
  check_constraints_schemas <- c("BookTown", "BrowserCookies", "CustomerOrder",
                                 "Employee", "Examination", "Flights", "iTrust",
                                 "NistWeather", "NistXTS748", "NistXTS749",
                                 "Person", "Products", "StudentResidence")

  newData <- d %>% dplyr::filter(casestudy %in% check_constraints_schemas)

  return(newData)
}

#' FUNCTION: get_check_schemas_mutants
#'
#' Filtering analysis data to get ony check constraints schemas.
#' @param d Data frame of analysis
#' @return A data frame of check constraints only schemas
#' @importFrom magrittr %>%
#' @export
get_check_schemas_mutants <- function(d) {
  check_constraints_schemas <- c("BookTown", "BrowserCookies", "CustomerOrder",
                                 "Employee", "Examination", "Flights", "iTrust",
                                 "NistWeather", "NistXTS748", "NistXTS749",
                                 "Person", "Products", "StudentResidence")

  newData <- d %>% dplyr::filter(schema %in% check_constraints_schemas)

  return(newData)
}

#' FUNCTION: domino_table_combined
#'
#' Timing, Coverage and mutation scores combained.
#' @param ana Data frame of analysis
#' @param mut Data frame of mutants
#' @return A data frame of check constraints only schemas or table LaTeX
#' @importFrom magrittr %>%
#' @export
domino_table_combined <- function(ana, mut, rtrn = "tex", mm = "median") {
  newAna <- dominoR::get_check_schemas_analysis(ana)
  newMut <- dominoR::get_check_schemas_mutants(mut)

  #coverage_domino <- dominoR::table_generator_coverage_domino(newAna, m = mm, rtrn = "data")
  #timing_domino <- dominoR::table_generator_timing_domino(newAna, m = mm, rtrn = "data")
  #mutants_domino <- dominoR::table_generator_mutation_score_domino(newMut, m = mm, rtrn = "data")

  coverage_domino <- dominoR::table_generator_coverage_domino_flipped(newAna, m = mm, rtrn = "data")
  timing_domino <- dominoR::table_generator_timing_domino_flipped(newAna, m = mm, rtrn = "data")
  mutants_domino <- dominoR::table_generator_mutation_score_domino_flipped(newMut, m = mm, rtrn = "data")

  colnames(coverage_domino) <- paste(colnames(coverage_domino), "cov", sep = "_")
  colnames(timing_domino) <- paste(colnames(timing_domino), "time", sep = "_")
  colnames(mutants_domino) <- paste(colnames(mutants_domino), "mutant", sep = "_")

  all_domino <- cbind(coverage_domino, timing_domino[2:7], mutants_domino[2:7])

  if (rtrn == "tex") {
    return(print(xtable::xtable(all_domino), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(all_domino)
  }
}

#' FUNCTION: table_generator_coverage_old
#'
#' Generates a latex table or a data frame for coverage table with effect size and U test.
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of coverages compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_coverage_old <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust")
  # Store the dataframe into another var
  d1 <- d
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((mean(coverage)), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((median(coverage)), 1), nsmall = 1))
  }
  #browser()
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("coverage"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:16]
  d <- d2[ , order(names(d2))]
  c <- d[1:5]
  c <- c[c(3,4,1,2,5)]
  a <- d[6:10]
  a <- a[c(3,4,1,2,5)]
  b <- d[11:15]
  b <- b[c(3,4,1,2,5)]
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get each generators
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # Effect size for PSQL
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                         (avm %>% dplyr::filter(dbms == "Postgres"))$coverage)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                          (avmd %>% dplyr::filter(dbms == "Postgres"))$coverage)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                          (rand %>% dplyr::filter(dbms == "Postgres"))$coverage)$size
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                           (dravm %>% dplyr::filter(dbms == "Postgres"))$coverage)$size

    dr_coverage <- (dr %>% dplyr::filter(dbms == "Postgres"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = dravm_coverage,
                           effect = postgres_dravm,
                           result = a[i,2])

    # get coverage
    avmr_coverage <- (avm %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,3] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmr_coverage,
                           effect = postgres_avm,
                           result = a[i,3])


    # get coverage
    avmd_coverage <- (avmd %>% dplyr::filter(dbms == "Postgres"))$coverage
    a[i,4] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmd_coverage,
                           effect = postgres_avmd,
                           result = a[i,4])


    # U-test Random vs DR
    rand_coverage <- (rand %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,5] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = rand_coverage,
                           effect = postgres_rand,
                           result = a[i,5])


    # get SQLite effect size
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                                       (avm %>% dplyr::filter(dbms == "SQLite"))$coverage)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                                        (avmd %>% dplyr::filter(dbms == "SQLite"))$coverage)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                                        (rand %>% dplyr::filter(dbms == "SQLite"))$coverage)$size
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                                         (dravm %>% dplyr::filter(dbms == "SQLite"))$coverage)$size

    # get coverage
    dr_coverage <- (dr %>% dplyr::filter(dbms == "SQLite"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "SQLite"))$coverage

    b[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = dravm_coverage,
                           effect = sqlite_dravm,
                           result = b[i,2])



    # U-test AVMR vs DR
    avmr_coverage <- (avm %>% dplyr::filter(dbms == "SQLite"))$coverage
    b[i,3] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmr_coverage,
                           effect = sqlite_avm,
                           result = b[i,3])



    # U-test AVMD vs DR
    avmd_coverage <- (avmd %>% dplyr::filter(dbms == "SQLite"))$coverage
    b[i,4] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmd_coverage,
                           effect = sqlite_avmd,
                           result = b[i,4])


    # U-test Random vs DR
    rand_coverage <- (rand %>% dplyr::filter(dbms == "SQLite"))$coverage
    b[i,5] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = rand_coverage,
                           effect = sqlite_rand,
                           result = b[i,5])

    # calculate effect size for coverage for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                                     (avm %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                                      (avmd %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                                      (rand %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                                       (dravm %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size

    # U-test DRAVM vs DR
    dr_coverage <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = dravm_coverage,
                           effect = hsql_dravm,
                           result = c[i,2])


    # U-test AVMR vs DR
    avmr_coverage <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,3] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmr_coverage,
                           effect = hsql_avm,
                           result = c[i,3])



    avmd_coverage <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,4] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmd_coverage,
                           effect = hsql_avmd,
                           result = c[i,4])


    rand_coverage <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,5] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = rand_coverage,
                           effect = hsql_rand,
                           result = c[i,5])


    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_coverage
#'
#' Generates a latex table or a data frame for coverage table with effect size and U test.
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of coverages compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_coverage <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust", datagenerator != "dravm")
  # Store the dataframe into another var
  d1 <- d
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((mean(coverage)), 0), nsmall = 0))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((median(coverage)), 0), nsmall = 0))
  }
  #browser()
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("coverage"))
  #browser()
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get each generators
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "random")

    # Effect size for PSQL
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$coverage)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$coverage)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                 (rand %>% dplyr::filter(dbms == "Postgres"))$coverage)$size

    dr_coverage <- (dr %>% dplyr::filter(dbms == "Postgres"))$coverage

    # get coverage
    avmr_coverage <- (avm %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmr_coverage,
                           effect = postgres_avm,
                           result = a[i,2])


    # get coverage
    avmd_coverage <- (avmd %>% dplyr::filter(dbms == "Postgres"))$coverage
    a[i,3] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmd_coverage,
                           effect = postgres_avmd,
                           result = a[i,3])


    # U-test Random vs DR
    rand_coverage <- (rand %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,4] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = rand_coverage,
                           effect = postgres_rand,
                           result = a[i,4])


    # get SQLite effect size
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$coverage)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$coverage)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                               (rand %>% dplyr::filter(dbms == "SQLite"))$coverage)$size

    # get coverage
    dr_coverage <- (dr %>% dplyr::filter(dbms == "SQLite"))$coverage

    # U-test AVMR vs DR
    avmr_coverage <- (avm %>% dplyr::filter(dbms == "SQLite"))$coverage
    b[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmr_coverage,
                           effect = sqlite_avm,
                           result = b[i,2])



    # U-test AVMD vs DR
    avmd_coverage <- (avmd %>% dplyr::filter(dbms == "SQLite"))$coverage
    b[i,3] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmd_coverage,
                           effect = sqlite_avmd,
                           result = b[i,3])


    # U-test Random vs DR
    rand_coverage <- (rand %>% dplyr::filter(dbms == "SQLite"))$coverage
    b[i,4] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = rand_coverage,
                           effect = sqlite_rand,
                           result = b[i,4])

    # calculate effect size for coverage for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                             (rand %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size

    # U-test DRAVM vs DR
    dr_coverage <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage

    # U-test AVMR vs DR
    avmr_coverage <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmr_coverage,
                           effect = hsql_avm,
                           result = c[i,2])



    avmd_coverage <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,3] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = avmd_coverage,
                           effect = hsql_avmd,
                           result = c[i,3])


    rand_coverage <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    c[i,4] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = rand_coverage,
                           effect = hsql_rand,
                           result = c[i,4])


    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_coverage_domino
#'
#' Generates a latex table or a data frame for coverage table with effect size and U test.
#' Only for domino AVM and Random
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of coverages compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_coverage_domino <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  d <- d %>% dplyr::filter(datagenerator %in% c("directedRandom", "dravm"))
  #d <- d %>% dplyr::filter(casestudy != "iTrust")
  #browser()
  # Store the dataframe into another var
  d1 <- d
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((mean(coverage)), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((median(coverage)), 1), nsmall = 1))
  }
  #browser()
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("coverage"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:7]
  d <- d2[ , order(names(d2))]
  c <- d[1:2]
  #c <- c[c(3,4,1,2,5)]
  a <- d[3:4]
  #a <- a[c(3,4,1,2,5)]
  b <- d[5:6]
  #b <- b[c(3,4,1,2,5)]
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get each generators
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # Effect size for PSQL
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                  (dravm %>% dplyr::filter(dbms == "Postgres"))$coverage)$size

    dr_coverage <- (dr %>% dplyr::filter(dbms == "Postgres"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = dravm_coverage,
                           effect = postgres_dravm,
                           result = a[i,2])

    # get SQLite effect size
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                                (dravm %>% dplyr::filter(dbms == "SQLite"))$coverage)$size

    # get coverage
    dr_coverage <- (dr %>% dplyr::filter(dbms == "SQLite"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "SQLite"))$coverage

    b[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = dravm_coverage,
                           effect = sqlite_dravm,
                           result = b[i,2])

    # calculate effect size for coverage for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                              (dravm %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size

    # U-test DRAVM vs DR
    dr_coverage <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$coverage

    c[i,2] = dominoR::comparing_sig(sample1 = dr_coverage,
                           sample2 = dravm_coverage,
                           effect = hsql_dravm,
                           result = c[i,2])


    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}


#' FUNCTION: table_generator_coverage_domino_flipped
#'
#' Generates a latex table or a data frame for coverage table with effect size and U test.
#' Only for domino AVM and Random
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of coverages compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_coverage_domino_flipped <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  d <- d %>% dplyr::filter(datagenerator %in% c("directedRandom", "dravm"))
  #d <- d %>% dplyr::filter(casestudy != "iTrust")
  #browser()
  # Store the dataframe into another var
  d1 <- d
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((mean(coverage)), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, coverage, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(coverage = format(round((median(coverage)), 1), nsmall = 1))
  }
  #browser()
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("coverage"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:7]
  d <- d2[ , order(names(d2))]
  c <- d[1:2]
  c <- c[c(2,1)]
  #c <- c[c(3,4,1,2,5)]
  a <- d[3:4]
  a <- a[c(2,1)]
  #a <- a[c(3,4,1,2,5)]
  b <- d[5:6]
  b <- b[c(2,1)]
  #b <- b[c(3,4,1,2,5)]
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get each generators
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # Effect size for PSQL
    postgres_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "Postgres"))$coverage,
                                                  (dr %>% dplyr::filter(dbms == "Postgres"))$coverage)$size

    dr_coverage <- (dr %>% dplyr::filter(dbms == "Postgres"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "Postgres"))$coverage

    a[i,2] = dominoR::comparing_sig(sample1 = dravm_coverage,
                                   sample2 = dr_coverage,
                                   effect = postgres_dravm,
                                   result = a[i,2])

    # get SQLite effect size
    sqlite_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "SQLite"))$coverage,
                                                (dr %>% dplyr::filter(dbms == "SQLite"))$coverage)$size

    # get coverage
    dr_coverage <- (dr %>% dplyr::filter(dbms == "SQLite"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "SQLite"))$coverage

    b[i,2] = dominoR::comparing_sig(sample1 = dravm_coverage,
                                   sample2 = dr_coverage,
                                   effect = sqlite_dravm,
                                   result = b[i,2])

    # calculate effect size for coverage for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "HyperSQL"))$coverage,
                                              (dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage)$size

    # U-test DRAVM vs DR
    dr_coverage <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$coverage
    dravm_coverage <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$coverage

    c[i,2] = dominoR::comparing_sig(sample1 = dravm_coverage,
                                   sample2 = dr_coverage,
                                   effect = hsql_dravm,
                                   result = c[i,2])


    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}


#' FUNCTION: table_generator_timing
#'
#' Generates a latex table or data frame for test generation timing table with effect size and U test.
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of test generation timing compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_timing <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust", datagenerator != "random")

  # copy values for Sig without transforming
  d3 <- d
  # Transform data with rounding down
  d1 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:16]
  d <- d2[ , order(names(d2))]
  c <- d[1:5]
  c <- c[c(3,4,1,2,5)]
  a <- d[6:10]
  a <- a[c(3,4,1,2,5)]
  b <- d[11:15]
  b <- b[c(3,4,1,2,5)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get each generators for transformed data
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")


    # Effect size for PSQL
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                         (avm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                          (avmd %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                          (rand %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                           (dravm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size

    # get generators for non-transformed
    drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmdp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    randp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    dr_time <- (drp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime
    avmr_time <- (dravmp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                           sample2 = avmr_time,
                           effect = postgres_dravm,
                           result = a[i,2])


    avmr_time <- (avmp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = postgres_avm,
                                  result = a[i,3])


    # AVM-D
    avmd_time <- (avmdp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmd_time,
                                  effect = postgres_avmd,
                                  result = a[i,4])


    # RANDOM
    rand_time <- (randp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,5] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = rand_time,
                                  effect = postgres_rand,
                                  result = a[i,5])

    # Effect size for SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                       (avm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                        (avmd %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                        (rand %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                         (dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (drp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime
    dravmp_time <- (dravmp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = dravmp_time,
                                  effect = sqlite_dravm,
                                  result = b[i,2])

    # DR vs AVM-R U-test
    avmr_time <- (avmp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = sqlite_avm,
                                  result = b[i,3])

    # AVM-D u-test
    avmd_time <- (avmdp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmd_time,
                                  effect = sqlite_avmd,
                                  result = b[i,4])

    # Random U-Test
    rand_time <- (randp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,5] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = rand_time,
                                  effect = sqlite_rand,
                                  result = b[i,5])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                                     (avm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                                      (avmd %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                                      (rand %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                                       (dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (drp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime
    dravm_time <- (dravmp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = dravm_time,
                                  effect = hsql_dravm,
                                  result = c[i,2])

    # U-Test avm-r vs dr
    avmr_time <- (avmp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = hsql_avm,
                                  result = c[i,3])

    # U-Test avm-d
    avmd_time <- (avmdp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmd_time,
                                  effect = hsql_avmd,
                                  result = c[i,4])



    # Rnadom u-test
    rand_time <- (randp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,5] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = rand_time,
                                  effect = hsql_rand,
                                  result = c[i,5])

    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_timing_nonRand_old
#'
#' Generates a latex table or data frame for test generation timing table with effect size and U test.
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of test generation timing compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_timing_nonRand_old <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust", datagenerator != "random")

  # copy values for Sig without transforming
  d3 <- d
  # Transform data with rounding down
  d1 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,4,1,2)]
  a <- d[5:8]
  a <- a[c(3,4,1,2)]
  b <- d[9:12]
  b <- b[c(3,4,1,2)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get each generators for transformed data
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")


    # Effect size for PSQL
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                  (dravm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size

    # get generators for non-transformed
    drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmdp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    randp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    dr_time <- (drp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime
    avmr_time <- (dravmp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = postgres_dravm,
                                          result = a[i,2])


    avmr_time <- (avmp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = postgres_avm,
                                          result = a[i,3])


    # AVM-D
    avmd_time <- (avmdp %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmd_time,
                                          effect = postgres_avmd,
                                          result = a[i,4])


    # Effect size for SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                (dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (drp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime
    dravmp_time <- (dravmp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = dravmp_time,
                                          effect = sqlite_dravm,
                                          result = b[i,2])

    # DR vs AVM-R U-test
    avmr_time <- (avmp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = sqlite_avm,
                                          result = b[i,3])

    # AVM-D u-test
    avmd_time <- (avmdp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmd_time,
                                          effect = sqlite_avmd,
                                          result = b[i,4])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                              (dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (drp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime
    dravm_time <- (dravmp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = dravm_time,
                                          effect = hsql_dravm,
                                          result = c[i,2])

    # U-Test avm-r vs dr
    avmr_time <- (avmp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = hsql_avm,
                                          result = c[i,3])

    # U-Test avm-d
    avmd_time <- (avmdp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmd_time,
                                          effect = hsql_avmd,
                                          result = c[i,4])

    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_timing_others
#'
#' Generates a latex table or data frame for test generation timing table with effect size and U test.
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of test generation timing compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_timing_others <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust", datagenerator != "dravm")

  # copy values for Sig without transforming
  d1 <- d
  # Transform data with rounding down
  # d3 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get generators for non-transformed
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # get each generators for transformed data
    # drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    # avmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    # avmdp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    # randp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    # dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")


    # Effect size for PSQL
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                 (rand %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size



    dr_time <- (dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    avmr_time <- (avm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = postgres_avm,
                                  result = a[i,2])


    # AVM-D
    avmd_time <- (avmd %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmd_time,
                                  effect = postgres_avmd,
                                  result = a[i,3])


    # RANDOM
    rand_time <- (rand %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = rand_time,
                                  effect = postgres_rand,
                                  result = a[i,4])

    # Effect size for SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                               (rand %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    # DR vs AVM-R U-test
    avmr_time <- (avm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = sqlite_avm,
                                  result = b[i,2])

    # AVM-D u-test
    avmd_time <- (avmd %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmd_time,
                                  effect = sqlite_avmd,
                                  result = b[i,3])

    # Random U-Test
    rand_time <- (rand %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = rand_time,
                                  effect = sqlite_rand,
                                  result = b[i,4])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                             (rand %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    # U-Test avm-r vs dr
    avmr_time <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = hsql_avm,
                                  result = c[i,2])

    # U-Test avm-d
    avmd_time <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmd_time,
                                  effect = hsql_avmd,
                                  result = c[i,3])



    # Rnadom u-test
    rand_time <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = rand_time,
                                  effect = hsql_rand,
                                  result = c[i,4])

    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_timing_nonRandom
#'
#' Generates a latex table or data frame for test generation timing table with effect size and U test.
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of test generation timing compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_timing_nonRandom <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust", datagenerator != "dravm")

  # copy values for Sig without transforming
  d1 <- d
  # Transform data with rounding down
  # d3 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get generators for non-transformed
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # get each generators for transformed data
    # drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    # avmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avs")
    # avmdp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "avsDefaults")
    # randp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "random")
    # dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")


    # Effect size for PSQL
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                 (rand %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size



    dr_time <- (dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    avmr_time <- (avm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = postgres_avm,
                                          result = a[i,2])


    # AVM-D
    avmd_time <- (avmd %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmd_time,
                                          effect = postgres_avmd,
                                          result = a[i,3])


    # RANDOM
    rand_time <- (rand %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = rand_time,
                                          effect = postgres_rand,
                                          result = a[i,4])

    # Effect size for SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                               (rand %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    # DR vs AVM-R U-test
    avmr_time <- (avm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = sqlite_avm,
                                          result = b[i,2])

    # AVM-D u-test
    avmd_time <- (avmd %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmd_time,
                                          effect = sqlite_avmd,
                                          result = b[i,3])

    # Random U-Test
    rand_time <- (rand %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = rand_time,
                                          effect = sqlite_rand,
                                          result = b[i,4])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                             (rand %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    # U-Test avm-r vs dr
    avmr_time <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmr_time,
                                          effect = hsql_avm,
                                          result = c[i,2])

    # U-Test avm-d
    avmd_time <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,3] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = avmd_time,
                                          effect = hsql_avmd,
                                          result = c[i,3])



    # Rnadom u-test
    rand_time <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,4] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                          sample2 = rand_time,
                                          effect = hsql_rand,
                                          result = c[i,4])

    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # Remove RANDPM techniques
  a$Postgres_random <- NULL
  b$SQLite_random <- NULL
  c$HyperSQL_random <- NULL

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_timing_domino
#'
#' Generates a latex table or data frame for test generation timing table with effect size and U test.
#' Only for domino AVM and Random
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of test generation timing compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_timing_domino <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust")
  d <- d %>% dplyr::filter(datagenerator %in% c("directedRandom", "dravm"))
  # Transform data with rounding down
  # d3 <- d
  # copy values for Sig without transforming
  d1 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:7]
  d <- d2[ , order(names(d2))]
  c <- d[1:2]
  # c <- c[c(3,4,1,2,5)]
  a <- d[3:4]
  #a <- a[c(3,4,1,2,5)]
  b <- d[5:6]
  #b <- b[c(3,4,1,2,5)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get generators for non-transformed
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # get each generators for transformed data
    #drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    #dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # Effect size for PSQL
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                  (dravm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size

    dr_time <- (dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime
    avmr_time <- (dravm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = avmr_time,
                                  effect = postgres_dravm,
                                  result = a[i,2])


    # Effect size for SQLite
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                (dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime
    dravmp_time <- (dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = dravmp_time,
                                  effect = sqlite_dravm,
                                  result = b[i,2])

    # Effect size for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                              (dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime
    dravm_time <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = dravm_time,
                                  effect = hsql_dravm,
                                  result = c[i,2])


    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_timing_domino_flipped
#'
#' Generates a latex table or data frame for test generation timing table with effect size and U test.
#' Only for domino AVM and Random
#' @param d Data frame of analysis
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of test generation timing compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_timing_domino_flipped <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust")
  d <- d %>% dplyr::filter(datagenerator %in% c("directedRandom", "dravm"))
  # Transform data with rounding down
  # d3 <- d
  # copy values for Sig without transforming
  d1 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>%
      dplyr::group_by(dbms, casestudy, datagenerator) %>%
      dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:7]
  d <- d2[ , order(names(d2))]
  c <- d[1:2]
  c <- c[c(2,1)]
  # c <- c[c(3,4,1,2,5)]
  a <- d[3:4]
  a <- a[c(2,1)]
  #a <- a[c(3,4,1,2,5)]
  b <- d[5:6]
  b <- b[c(2,1)]
  #b <- b[c(3,4,1,2,5)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    # get generators for non-transformed
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # get each generators for transformed data
    #drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "directedRandom")
    #dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "dravm")

    # Effect size for PSQL
    postgres_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime,
                                                  (dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime)$size

    dr_time <- (dr %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime
    avmr_time <- (dravm %>% dplyr::filter(dbms == "Postgres"))$testgenerationtime

    a[i,2] = dominoR::comparing_sig_timing(sample1 = avmr_time,
                                          sample2 = dr_time,
                                          effect = postgres_dravm,
                                          result = a[i,2])


    # Effect size for SQLite
    sqlite_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                (dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime
    dravmp_time <- (dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dravmp_time,
                                          sample2 = dr_time,
                                          effect = sqlite_dravm,
                                          result = b[i,2])

    # Effect size for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                              (dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime
    dravm_time <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dravm_time,
                                          sample2 = dr_time,
                                          effect = hsql_dravm,
                                          result = c[i,2])


    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

table_generator_timing_hsql <- function(d, rtrn = "tex", m = "median") {
  # Arrange dataframe by case study
  d <- d %>% dplyr::arrange(casestudy)
  #d <- d %>% dplyr::filter(casestudy != "iTrust")

  # copy values for Sig without transforming
  d3 <- d
  # Transform data with rounding down
  d1 <- d
  # d1 <- dominoR::transform_execution_times_for_threshold(d, 1000)
  # generate a DF for mean or median
  if (m == "mean") {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((mean(testgenerationtime) / 1000), 2), nsmall = 2))
  } else {
    d <- d %>% dplyr::select(dbms, casestudy, datagenerator, testgenerationtime, randomseed) %>% dplyr::group_by(dbms, casestudy, datagenerator) %>% dplyr::summarise(testgenerationtime = format(round((median(testgenerationtime) / 1000), 2), nsmall = 2))
  }
  # filp the data frame
  d <- reshape2::dcast(d, casestudy ~ dbms + datagenerator, value.var=c("testgenerationtime"))
  # get header
  a1 <- d[1]
  # Split by DBMS
  d2 <- d[2:5]
  d <- d2[ , order(names(d2))]
  # c <- d[1:5]
  # c <- c[c(3,4,1,2,5)]
  #browser()
  c <- d[1:2]
  c <- c[c(2,1)]
  b <- d[3:4]
  b <- b[c(2,1)]
  # get nunber of rows and itrate through them
  numberOfRows <- nrow(d)
  # change the schemas from fectors to char
  a1$casestudy <- as.character(a1$casestudy)
  for (i in 1:numberOfRows) {
    schema <- a1[i,]
    #browser()
    # get each generators for transformed data
    dr <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "concentroRandom")
    dravm <- d1 %>% dplyr::filter(casestudy == schema, datagenerator == "concentroAVS")

    # get generators for non-transformed
    drp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "concentroRandom")
    dravmp <- d3 %>% dplyr::filter(casestudy == schema, datagenerator == "concentroAVS")

    # Effect size for SQLite
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime,
                                                         (dravm %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime)$size

    # DR vs AVM-R U-test
    dr_time <- (drp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime
    dravmp_time <- (dravmp %>% dplyr::filter(dbms == "SQLite"))$testgenerationtime

    b[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = dravmp_time,
                                  effect = sqlite_dravm,
                                  result = b[i,2])



    # Effect size for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime,
                                                       (dravm %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime)$size

    # U-Test avm-r vs dr
    dr_time <- (drp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime
    dravm_time <- (dravmp %>% dplyr::filter(dbms == "HyperSQL"))$testgenerationtime

    c[i,2] = dominoR::comparing_sig_timing(sample1 = dr_time,
                                  sample2 = dravm_time,
                                  effect = hsql_dravm,
                                  result = c[i,2])
    # for latex purposes
    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  # Combain data
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  #return(d)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_mutation_score
#'
#' Generates a latex table or data frame for mutation score per schema table with effect size and U test.
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per schema compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutation_score <- function(d, rtrn = "tex", m = "median") {
  # ordering mutants per run
  #d <- d %>% dplyr::filter(schema != "iTrust")
  d <- ordering_mutants_per_schema(d)
  # copying data frame so it can be compared for A12 and U-test
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(mean(mutationScore), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(median(mutationScore), 1), nsmall = 1))
  }
  # Reshaping data frame
  d <- reshape2::dcast(d, schema ~ dbms + datagenerator, value.var=c("mutationScore"))
  a1 <- d[1]
  # Splitting data frame per DBMS
  d2 <- d[2:16]
  d <- d2[ , order(names(d2))]
  c <- d[1:5]
  c <- c[c(3,4,1,2,5)]
  a <- d[6:10]
  a <- a[c(3,4,1,2,5)]
  b <- d[11:15]
  b <- b[c(3,4,1,2,5)]
  # Schemas changed to
  a1$schema <- as.character(a1$schema)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema1 <- a1[i,]
    dr <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "random")
    dravm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "dravm")


    # PSQL A12
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                         (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                          (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                          (rand %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                           (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size

    # DR vs dravm
    dr_mutation <- (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                          sample2 = dravm_mutation,
                          effect = postgres_dravm,
                          result = a[i,2])

    # DR vs AVM-r
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmr_mutation,
                           effect = postgres_avm,
                           result = a[i,3])

    # AVM-D vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmd_mutation,
                           effect = postgres_avmd,
                           result = a[i,4])


    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,5] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = rand_mutation,
                           effect = postgres_rand,
                           result = a[i,5])

    # A12 SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                       (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                        (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                        (rand %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                         (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size

    # Dr vs DRAVM
    dr_mutation <- (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = dravm_mutation,
                           effect = sqlite_dravm,
                           result = b[i,2])

    # Dr vs AVM-R
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmr_mutation,
                           effect = sqlite_avm,
                           result = b[i,3])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmd_mutation,
                           effect = sqlite_avmd,
                           result = b[i,4])

    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,5] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = rand_mutation,
                           effect = sqlite_rand,
                           result = b[i,5])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                                     (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                                      (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                                      (rand %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                                       (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size

    # DR vs AVMR
    dr_mutation <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = dravm_mutation,
                           effect = hsql_dravm,
                           result = c[i,2])
    # DR vs AVMR
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmr_mutation,
                           effect = hsql_avm,
                           result = c[i,3])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmd_mutation,
                           effect = hsql_avmd,
                           result = c[i,4])

    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,5] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = rand_mutation,
                           effect = hsql_rand,
                           result = c[i,5])

    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}


#' FUNCTION: table_generator_mutation_score_nonrnd
#'
#' Generates a latex table or data frame for mutation score per schema table with effect size and U test.
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per schema compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutation_score_nonrnd <- function(d, rtrn = "tex", m = "median") {
  # ordering mutants per run
  #d <- d %>% dplyr::filter(schema != "iTrust")
  d <- ordering_mutants_per_schema(d)
  d <- d %>% dplyr::filter(datagenerator != "random")
  # copying data frame so it can be compared for A12 and U-test
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(mean(mutationScore), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(median(mutationScore), 1), nsmall = 1))
  }
  # Reshaping data frame
  d <- reshape2::dcast(d, schema ~ dbms + datagenerator, value.var=c("mutationScore"))
  a1 <- d[1]
  # Splitting data frame per DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,4,1,2)]
  a <- d[5:8]
  a <- a[c(3,4,1,2)]
  b <- d[9:12]
  b <- b[c(3,4,1,2)]
  # Schemas changed to
  a1$schema <- as.character(a1$schema)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema1 <- a1[i,]
    dr <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avsDefaults")
    dravm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "dravm")


    # PSQL A12
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                  (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size

    # DR vs dravm
    dr_mutation <- (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = dravm_mutation,
                                   effect = postgres_dravm,
                                   result = a[i,2])

    # DR vs AVM-r
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmr_mutation,
                                   effect = postgres_avm,
                                   result = a[i,3])

    # AVM-D vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmd_mutation,
                                   effect = postgres_avmd,
                                   result = a[i,4])

    # A12 SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size

    # Dr vs DRAVM
    dr_mutation <- (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = dravm_mutation,
                                   effect = sqlite_dravm,
                                   result = b[i,2])

    # Dr vs AVM-R
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmr_mutation,
                                   effect = sqlite_avm,
                                   result = b[i,3])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmd_mutation,
                                   effect = sqlite_avmd,
                                   result = b[i,4])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                              (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size

    # DR vs AVMR
    dr_mutation <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = dravm_mutation,
                                   effect = hsql_dravm,
                                   result = c[i,2])
    # DR vs AVMR
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmr_mutation,
                                   effect = hsql_avm,
                                   result = c[i,3])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmd_mutation,
                                   effect = hsql_avmd,
                                   result = c[i,4])

    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_mutation_score_others
#'
#' Generates a latex table or data frame for mutation score per schema table with effect size and U test.
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per schema compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutation_score_others <- function(d, rtrn = "tex", m = "median") {
  # ordering mutants per run
  #d <- d %>% dplyr::filter(schema != "iTrust", datagenerator != "dravm")
  d <- ordering_mutants_per_schema_others(d)
  # copying data frame so it can be compared for A12 and U-test
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(mean(mutationScore), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(median(mutationScore), 1), nsmall = 1))
  }
  # Reshaping data frame
  d <- reshape2::dcast(d, schema ~ dbms + datagenerator, value.var=c("mutationScore"))
  a1 <- d[1]
  # Splitting data frame per DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # Schemas changed to
  a1$schema <- as.character(a1$schema)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema1 <- a1[i,]
    dr <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "random")


    # PSQL A12
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                 (rand %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size

    # DR vs dravm
    dr_mutation <- (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    # DR vs AVM-r
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmr_mutation,
                           effect = postgres_avm,
                           result = a[i,2])

    # AVM-D vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmd_mutation,
                           effect = postgres_avmd,
                           result = a[i,3])


    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = rand_mutation,
                           effect = postgres_rand,
                           result = a[i,4])

    # A12 SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                               (rand %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size

    # Dr vs DRAVM
    dr_mutation <- (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    # Dr vs AVM-R
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmr_mutation,
                           effect = sqlite_avm,
                           result = b[i,2])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmd_mutation,
                           effect = sqlite_avmd,
                           result = b[i,3])

    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = rand_mutation,
                           effect = sqlite_rand,
                           result = b[i,4])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                             (rand %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size

    # DR vs AVMR
    dr_mutation <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore
    # DR vs AVMR
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmr_mutation,
                           effect = hsql_avm,
                           result = c[i,2])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = avmd_mutation,
                           effect = hsql_avmd,
                           result = c[i,3])

    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = rand_mutation,
                           effect = hsql_rand,
                           result = c[i,4])

    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_mutation_score_nonRandom
#'
#' Generates a latex table or data frame for mutation score per schema table with effect size and U test.
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per schema compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutation_score_nonRandom <- function(d, rtrn = "tex", m = "median") {
  # ordering mutants per run
  #d <- d %>% dplyr::filter(schema != "iTrust", datagenerator != "dravm")
  d <- ordering_mutants_per_schema_others(d)
  # copying data frame so it can be compared for A12 and U-test
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(mean(mutationScore), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>% dplyr::summarise(mutationScore = format(round(median(mutationScore), 1), nsmall = 1))
  }
  # Reshaping data frame
  d <- reshape2::dcast(d, schema ~ dbms + datagenerator, value.var=c("mutationScore"))
  a1 <- d[1]
  # Splitting data frame per DBMS
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # Schemas changed to
  a1$schema <- as.character(a1$schema)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema1 <- a1[i,]
    dr <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "directedRandom")
    avm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avs")
    avmd <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "avsDefaults")
    rand <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "random")


    # PSQL A12
    postgres_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                 (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size
    postgres_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                 (rand %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size

    # DR vs dravm
    dr_mutation <- (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    # DR vs AVM-r
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmr_mutation,
                                   effect = postgres_avm,
                                   result = a[i,2])

    # AVM-D vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmd_mutation,
                                   effect = postgres_avmd,
                                   result = a[i,3])


    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = rand_mutation,
                                   effect = postgres_rand,
                                   result = a[i,4])

    # A12 SQLite
    sqlite_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                              (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                               (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size
    sqlite_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                               (rand %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size

    # Dr vs DRAVM
    dr_mutation <- (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    # Dr vs AVM-R
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmr_mutation,
                                   effect = sqlite_avm,
                                   result = b[i,2])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmd_mutation,
                                   effect = sqlite_avmd,
                                   result = b[i,3])

    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = rand_mutation,
                                   effect = sqlite_rand,
                                   result = b[i,4])

    # Effect size for HSQL
    hsql_avm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                            (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_avmd <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                             (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size
    hsql_rand <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                             (rand %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size

    # DR vs AVMR
    dr_mutation <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore
    # DR vs AVMR
    avmr_mutation <- (avm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmr_mutation,
                                   effect = hsql_avm,
                                   result = c[i,2])

    # AVMD vs DR
    avmd_mutation <- (avmd %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,3] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = avmd_mutation,
                                   effect = hsql_avmd,
                                   result = c[i,3])

    # Random vs DR
    rand_mutation <- (rand %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,4] = dominoR::comparing_sig(sample1 = dr_mutation,
                                   sample2 = rand_mutation,
                                   effect = hsql_rand,
                                   result = c[i,4])

    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]

  # Remove RANDPM techniques
  a$Postgres_random <- NULL
  b$SQLite_random <- NULL
  c$HyperSQL_random <- NULL

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_mutation_score_domino
#'
#' Generates a latex table or data frame for mutation score per schema table with effect size and U test.
#' Only for domino AVM and Random
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per schema compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutation_score_domino <- function(d, rtrn = "tex", m = "median") {
  # ordering mutants per run
  #d <- d %>% dplyr::filter(schema != "iTrust")
  d <- d %>% dplyr::filter(datagenerator %in% c("directedRandom", "dravm"))
  d <- ordering_mutants_per_schema_domino(d)
  # copying data frame so it can be compared for A12 and U-test
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>%
      dplyr::summarise(mutationScore = format(round(mean(mutationScore), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>%
      dplyr::summarise(mutationScore = format(round(median(mutationScore), 1), nsmall = 1))
  }
  # Reshaping data frame
  d <- reshape2::dcast(d, schema ~ dbms + datagenerator, value.var=c("mutationScore"))
  a1 <- d[1]
  # Splitting data frame per DBMS
  d2 <- d[2:7]
  d <- d2[ , order(names(d2))]
  c <- d[1:2]
  #c <- c[c(3,4,1,2,5)]
  a <- d[3:4]
  #a <- a[c(3,4,1,2,5)]
  b <- d[5:6]
  #b <- b[c(3,4,1,2,5)]
  # Schemas changed to
  a1$schema <- as.character(a1$schema)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema1 <- a1[i,]
    dr <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "directedRandom")
    dravm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "dravm")


    # PSQL A12
    postgres_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                  (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size

    # DR vs dravm
    dr_mutation <- (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = dravm_mutation,
                           effect = postgres_dravm,
                           result = a[i,2])


    # A12 SQLite
    sqlite_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size

    # Dr vs DRAVM
    dr_mutation <- (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = dravm_mutation,
                           effect = sqlite_dravm,
                           result = b[i,2])


    # Effect size for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                              (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size

    # DR vs AVMR
    dr_mutation <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,2] = dominoR::comparing_sig(sample1 = dr_mutation,
                           sample2 = dravm_mutation,
                           effect = hsql_dravm,
                           result = c[i,2])

    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}


#' FUNCTION: table_generator_mutation_score_domino_flipped
#'
#' Generates a latex table or data frame for mutation score per schema table with effect size and U test.
#' Only for domino AVM and Random
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per schema compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutation_score_domino_flipped <- function(d, rtrn = "tex", m = "median") {
  # ordering mutants per run
  #d <- d %>% dplyr::filter(schema != "iTrust")
  d <- d %>% dplyr::filter(datagenerator %in% c("directedRandom", "dravm"))
  d <- ordering_mutants_per_schema_domino(d)
  # copying data frame so it can be compared for A12 and U-test
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>%
      dplyr::summarise(mutationScore = format(round(mean(mutationScore), 1), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(schema, datagenerator, dbms)  %>%
      dplyr::summarise(mutationScore = format(round(median(mutationScore), 1), nsmall = 1))
  }
  # Reshaping data frame
  d <- reshape2::dcast(d, schema ~ dbms + datagenerator, value.var=c("mutationScore"))
  a1 <- d[1]
  # Splitting data frame per DBMS
  d2 <- d[2:7]
  d <- d2[ , order(names(d2))]
  c <- d[1:2]
  c <- c[c(2,1)]
  #c <- c[c(3,4,1,2,5)]
  a <- d[3:4]
  a <- a[c(2,1)]
  #a <- a[c(3,4,1,2,5)]
  b <- d[5:6]
  b <- b[c(2,1)]
  #b <- b[c(3,4,1,2,5)]
  # Schemas changed to
  a1$schema <- as.character(a1$schema)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    schema1 <- a1[i,]
    dr <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "directedRandom")
    dravm <- d1 %>% dplyr::filter(schema == schema1, datagenerator == "dravm")


    # PSQL A12
    postgres_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore,
                                                  (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore)$size

    # DR vs dravm
    dr_mutation <- (dr %>% dplyr::filter(dbms == "Postgres"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "Postgres"))$mutationScore

    a[i,2] = dominoR::comparing_sig(sample1 = dravm_mutation,
                                   sample2 = dr_mutation,
                                   effect = postgres_dravm,
                                   result = a[i,2])


    # A12 SQLite
    sqlite_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore,
                                                (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore)$size

    # Dr vs DRAVM
    dr_mutation <- (dr %>% dplyr::filter(dbms == "SQLite"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "SQLite"))$mutationScore

    b[i,2] = dominoR::comparing_sig(sample1 = dravm_mutation,
                                   sample2 = dr_mutation,
                                   effect = sqlite_dravm,
                                   result = b[i,2])


    # Effect size for HSQL
    hsql_dravm <- dominoR::effectsize_accurate((dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore,
                                              (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore)$size

    # DR vs AVMR
    dr_mutation <- (dr %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore
    dravm_mutation <- (dravm %>% dplyr::filter(dbms == "HyperSQL"))$mutationScore

    c[i,2] = dominoR::comparing_sig(sample1 = dravm_mutation,
                                   sample2 = dr_mutation,
                                   effect = hsql_dravm,
                                   result = c[i,2])

    if (a1[i,] == "NistXTS749") {
      a1[i,] <- "NistXTSNine"
    }
    if (a1[i,] == "Iso3166") {
      a1[i,] <- "Isoiii"
    }
    if (a1[i,] == "IsoFlav_R2") {
      a1[i,] <- "IsoFlav"
    }
    if (a1[i,] == "NistDML181") {
      a1[i,] <- "NistDMLi"
    }
    if (a1[i,] == "NistDML182") {
      a1[i,] <- "NistDMLii"
    }
    if (a1[i,] == "NistDML183") {
      a1[i,] <- "NistDMLiii"
    }
    if (a1[i,] == "NistXTS748") {
      a1[i,] <- "NistXTSEight"
    }
    a1[i,] <- paste("\\", a1[i,], "ForTable", sep = "")
  }
  #a <- a[c(1,4,3,2)]
  #b <- b[c(1,4,3,2)]
  #c <- c[c(1,4,3,2)]
  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)
  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: table_generator_mutant_operators
#'
#' Generates a latex table or data frame for mutation operators table with effect size and U test.
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per operator compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutant_operators <- function(d, rtrn = "tex", m = "median") {
  # Order mutants per run
  d <- ordering_mutants_per_operator(d)
  # copying data before reshaping
  d1 <- d
  if (m == "mean") {
    a <- d %>% dplyr::group_by(dbms, datagenerator, operator) %>% dplyr::summarise(value = format(mean(mutationScore), nsmall = 1))
  } else {
    a <- d %>% dplyr::group_by(dbms, datagenerator, operator) %>% dplyr::summarise(value = format(median(mutationScore), nsmall = 1))
  }
  d <- reshape2::dcast(a,  operator ~ dbms + datagenerator)
  # Spliting df per DBMS
  a1 <- d[1]
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # Transoferming operater columns to strings
  a1$operator <- as.character(a1$operator)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    selected_operator <- a1[i,]
    # Selecting data per DBMS and datagenerator
    postgres_dr <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "directedRandom")
    postgres_avm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "avs")
    postgres_avmd <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "avsDefaults")
    postgres_rand <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "random")
    sqlite_dr <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "directedRandom")
    sqlite_avm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "avs")
    sqlite_avmd <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "avsDefaults")
    sqlite_rand <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "random")
    hsql_dr <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "directedRandom")
    hsql_avm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "avs")
    hsql_avmd <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "avsDefaults")
    hsql_rand <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "random")

    # A12 for PSQL
    postgres_avm_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_avm$mutationScore)$size
    postgres_avmd_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_avmd$mutationScore)$size
    postgres_rand_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_rand$mutationScore)$size

    # Sig for DR vs AVMR
    dr_mutation <- postgres_dr$mutationScore
    avmr_mutation <- postgres_avm$mutationScore

    p1 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      a[i,2] = paste("\\textbf{",a[i,2],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      a[i,2] = paste("\\textit{",a[i,2],"}", sep = "")
    } else {
    }

    if (postgres_avm_effectsize == "large") {
      a[i,2] = paste("$^{\\ast\\ast\\ast}$",a[i,2], sep = "")
    } else if (postgres_avm_effectsize == "medium") {
      a[i,2] = paste("$^{\\ast\\ast}$",a[i,2], sep = "")
    } else if (postgres_avm_effectsize == "small") {
      a[i,2] = paste("$^{\\ast}$",a[i,2], sep = "")
    } else {
    }

    # U-test AVMD vs DR
    avmd_mutation <- postgres_avmd$mutationScore

    p1 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      a[i,3] = paste("\\textbf{",a[i,3],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      a[i,3] = paste("\\textit{",a[i,3],"}", sep = "")
    } else {
    }

    if (postgres_avmd_effectsize == "large") {
      a[i,3] = paste("$^{\\ast\\ast\\ast}$",a[i,3], sep = "")
    } else if (postgres_avmd_effectsize == "medium") {
      a[i,3] = paste("$^{\\ast\\ast}$",a[i,3], sep = "")
    } else if (postgres_avmd_effectsize == "small") {
      a[i,3] = paste("$^{\\ast}$",a[i,3], sep = "")
    } else {

    }

    # U-test for Random vs DR
    rand_mutation <- postgres_rand$mutationScore

    p1 <- wilcox.test(dr_mutation, rand_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, rand_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      a[i,4] = paste("\\textbf{",a[i,4],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      a[i,4] = paste("\\textit{",a[i,4],"}", sep = "")
    } else {
    }

    if (postgres_rand_effectsize == "large") {
      a[i,4] = paste("$^{\\ast\\ast\\ast}$",a[i,4], sep = "")
    } else if (postgres_rand_effectsize == "medium") {
      a[i,4] = paste("$^{\\ast\\ast}$",a[i,4], sep = "")
    } else if (postgres_rand_effectsize == "small") {
      a[i,4] = paste("$^{\\ast}$",a[i,4], sep = "")
    } else {

    }

    # A12 for SQLite
    sqlite_avm_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_avm$mutationScore)$size
    sqlite_avmd_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_avmd$mutationScore)$size
    sqlite_rand_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_rand$mutationScore)$size

    # U-test dr vs avm-r
    dr_mutation <- sqlite_dr$mutationScore
    avmr_mutation <- sqlite_avm$mutationScore

    p1 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      b[i,2] = paste("\\textbf{",b[i,2],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      b[i,2] = paste("\\textit{",b[i,2],"}", sep = "")
    } else {
    }

    if (sqlite_avm_effectsize == "large") {
      b[i,2] = paste("$^{\\ast\\ast\\ast}$",b[i,2], sep = "")
    } else if (sqlite_avm_effectsize == "medium") {
      b[i,2] = paste("$^{\\ast\\ast}$",b[i,2], sep = "")
    } else if (sqlite_avm_effectsize == "small") {
      b[i,2] = paste("$^{\\ast}$",b[i,2], sep = "")
    } else {

    }

    # U-test AVMD vs dr
    avmd_mutation <- sqlite_avmd$mutationScore

    p1 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      b[i,3] = paste("\\textbf{",b[i,3],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      b[i,3] = paste("\\textit{",b[i,3],"}", sep = "")
    } else {
    }

    if (sqlite_avmd_effectsize == "large") {
      b[i,3] = paste("$^{\\ast\\ast\\ast}$",b[i,3], sep = "")
    } else if (sqlite_avmd_effectsize == "medium") {
      b[i,3] = paste("$^{\\ast\\ast}$",b[i,3], sep = "")
    } else if (sqlite_avmd_effectsize == "small") {
      b[i,3] = paste("$^{\\ast}$",b[i,3], sep = "")
    } else {

    }

    # U-test for Random vs DR
    rand_mutation <- sqlite_rand$mutationScore

    p1 <- wilcox.test(dr_mutation, rand_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, rand_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      b[i,4] = paste("\\textbf{",b[i,4],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      b[i,4] = paste("\\textit{",b[i,4],"}", sep = "")
    } else {
    }

    if (sqlite_rand_effectsize == "large") {
      b[i,4] = paste("$^{\\ast\\ast\\ast}$",b[i,4], sep = "")
    } else if (sqlite_rand_effectsize == "medium") {
      b[i,4] = paste("$^{\\ast\\ast}$",b[i,4], sep = "")
    } else if (sqlite_rand_effectsize == "small") {
      b[i,4] = paste("$^{\\ast}$",b[i,4], sep = "")
    } else {

    }

    # A12 for HSQL
    hsql_avm_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_avm$mutationScore)$size
    hsql_avmd_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_avmd$mutationScore)$size
    hsql_rand_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_rand$mutationScore)$size

    # U-test for DR vs AVM-R
    dr_mutation <- hsql_dr$mutationScore
    avmr_mutation <- hsql_avm$mutationScore

    p1 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      c[i,2] = paste("\\textbf{",c[i,2],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      c[i,2] = paste("\\textit{",c[i,2],"}", sep = "")
    } else {
    }

    if (hsql_avm_effectsize == "large") {
      c[i,2] = paste("$^{\\ast\\ast\\ast}$",c[i,2], sep = "")
    } else if (hsql_avm_effectsize == "medium") {
      c[i,2] = paste("$^{\\ast\\ast}$",c[i,2], sep = "")
    } else if (hsql_avm_effectsize == "small") {
      c[i,2] = paste("$^{\\ast}$",c[i,2], sep = "")
    } else {

    }

    # AVM-D vs DR U-test
    avmd_mutation <- hsql_avmd$mutationScore
    p1 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      c[i,3] = paste("\\textbf{",c[i,3],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      c[i,3] = paste("\\textit{",c[i,3],"}", sep = "")
    } else {
    }

    if (hsql_avmd_effectsize == "large") {
      c[i,3] = paste("$^{\\ast\\ast\\ast}$",c[i,3], sep = "")
    } else if (hsql_avmd_effectsize == "medium") {
      c[i,3] = paste("$^{\\ast\\ast}$",c[i,3], sep = "")
    } else if (hsql_avmd_effectsize == "small") {
      c[i,3] = paste("$^{\\ast}$",c[i,3], sep = "")
    } else {

    }

    # Random vs Dr U-test
    rand_mutation <- hsql_rand$mutationScore

    p1 <- wilcox.test(dr_mutation, rand_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, rand_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      c[i,4] = paste("\\textbf{",c[i,4],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      c[i,4] = paste("\\textit{",c[i,4],"}", sep = "")
    } else {
    }

    if (hsql_rand_effectsize == "large") {
      c[i,4] = paste("$^{\\ast\\ast\\ast}$",c[i,4], sep = "")
    } else if (hsql_rand_effectsize == "medium") {
      c[i,4] = paste("$^{\\ast\\ast}$",c[i,4], sep = "")
    } else if (hsql_rand_effectsize == "small") {
      c[i,4] = paste("$^{\\ast}$",c[i,4], sep = "")
    } else {

    }

    # For latex table
    a1[i,] <- paste("\\", a1[i,], sep = "")


  }
  a <- a[c(1,4,3,2)]
  b <- b[c(1,4,3,2)]
  c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)

  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}


#' FUNCTION: table_generator_mutant_operators_domino
#'
#' Generates a latex table or data frame for mutation operators table with effect size and U test.
#' @param d Data frame of mutants
#' @param rtrn Latex (tex) or a data frame (data)
#' @param m Results shown as median or mean
#' @return A A12 effect size and U-test of mutation score per operator compared pair wise
#' @importFrom magrittr %>%
#' @export
table_generator_mutant_operators_domino <- function(d, rtrn = "tex", m = "median") {
  # Order mutants per run
  d <- ordering_mutants_per_operator_domino(d)
  #browser()

  # copying data before reshaping
  d1 <- d
  if (m == "mean") {
    d <- d %>% dplyr::group_by(dbms, datagenerator, operator) %>% dplyr::summarise(value = format(mean(mutationScore), nsmall = 1))
  } else {
    d <- d %>% dplyr::group_by(dbms, datagenerator, operator) %>% dplyr::summarise(value = format(median(mutationScore), nsmall = 1))
  }
  d <- reshape2::dcast(d,  operator ~ dbms + datagenerator)
  # Spliting df per DBMS
  a1 <- d[1]
  d2 <- d[2:13]
  d <- d2[ , order(names(d2))]
  c <- d[1:4]
  c <- c[c(3,1,2,4)]
  a <- d[5:8]
  a <- a[c(3,1,2,4)]
  b <- d[9:12]
  b <- b[c(3,1,2,4)]
  # Transoferming operater columns to strings
  a1$operator <- as.character(a1$operator)
  numberOfRows <- nrow(d)
  for (i in 1:numberOfRows) {
    selected_operator <- a1[i,]
    # Selecting data per DBMS and datagenerator
    postgres_dr <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "directedRandom")
    postgres_avm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "avs")
    postgres_avmd <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "avsDefaults")
    postgres_rand <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "random")
    postgres_dravm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "Postgres", datagenerator == "dravm")

    sqlite_dr <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "directedRandom")
    sqlite_avm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "avs")
    sqlite_avmd <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "avsDefaults")
    sqlite_rand <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "random")
    sqlite_dravm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "SQLite", datagenerator == "dravm")

    hsql_dr <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "directedRandom")
    hsql_avm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "avs")
    hsql_avmd <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "avsDefaults")
    hsql_rand <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "random")
    hsql_dravm <- d1 %>% dplyr::filter(operator == selected_operator, dbms == "HyperSQL", datagenerator == "dravm")

    # A12 for PSQL
    postgres_avm_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_avm$mutationScore)$size
    postgres_avmd_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_avmd$mutationScore)$size
    postgres_rand_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_rand$mutationScore)$size
    postgres_dravm_effectsize <- dominoR::effectsize_accurate(postgres_dr$mutationScore, postgres_dravm$mutationScore)$size

    # Sig for DR vs AVMR
    dr_mutation <- postgres_dr$mutationScore
    avmr_mutation <- postgres_avm$mutationScore

    p1 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      a[i,2] = paste("\\textbf{",a[i,2],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      a[i,2] = paste("\\textit{",a[i,2],"}", sep = "")
    } else {
    }

    if (postgres_avm_effectsize == "large") {
      a[i,2] = paste("$^{\\ast\\ast\\ast}$",a[i,2], sep = "")
    } else if (postgres_avm_effectsize == "medium") {
      a[i,2] = paste("$^{\\ast\\ast}$",a[i,2], sep = "")
    } else if (postgres_avm_effectsize == "small") {
      a[i,2] = paste("$^{\\ast}$",a[i,2], sep = "")
    } else {
    }

    # U-test AVMD vs DR
    avmd_mutation <- postgres_avmd$mutationScore

    p1 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      a[i,3] = paste("\\textbf{",a[i,3],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      a[i,3] = paste("\\textit{",a[i,3],"}", sep = "")
    } else {
    }

    if (postgres_avmd_effectsize == "large") {
      a[i,3] = paste("$^{\\ast\\ast\\ast}$",a[i,3], sep = "")
    } else if (postgres_avmd_effectsize == "medium") {
      a[i,3] = paste("$^{\\ast\\ast}$",a[i,3], sep = "")
    } else if (postgres_avmd_effectsize == "small") {
      a[i,3] = paste("$^{\\ast}$",a[i,3], sep = "")
    } else {

    }

    # U-test for Random vs DR
    rand_mutation <- postgres_rand$mutationScore

    p1 <- wilcox.test(dr_mutation, rand_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, rand_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      a[i,4] = paste("\\textbf{",a[i,4],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      a[i,4] = paste("\\textit{",a[i,4],"}", sep = "")
    } else {
    }

    if (postgres_rand_effectsize == "large") {
      a[i,4] = paste("$^{\\ast\\ast\\ast}$",a[i,4], sep = "")
    } else if (postgres_rand_effectsize == "medium") {
      a[i,4] = paste("$^{\\ast\\ast}$",a[i,4], sep = "")
    } else if (postgres_rand_effectsize == "small") {
      a[i,4] = paste("$^{\\ast}$",a[i,4], sep = "")
    } else {

    }

    # A12 for SQLite
    sqlite_avm_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_avm$mutationScore)$size
    sqlite_avmd_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_avmd$mutationScore)$size
    sqlite_rand_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_rand$mutationScore)$size
    sqlite_dravm_effectsize <- dominoR::effectsize_accurate(sqlite_dr$mutationScore, sqlite_dravm$mutationScore)$size

    # U-test dr vs avm-r
    dr_mutation <- sqlite_dr$mutationScore
    avmr_mutation <- sqlite_avm$mutationScore

    p1 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      b[i,2] = paste("\\textbf{",b[i,2],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      b[i,2] = paste("\\textit{",b[i,2],"}", sep = "")
    } else {
    }

    if (sqlite_avm_effectsize == "large") {
      b[i,2] = paste("$^{\\ast\\ast\\ast}$",b[i,2], sep = "")
    } else if (sqlite_avm_effectsize == "medium") {
      b[i,2] = paste("$^{\\ast\\ast}$",b[i,2], sep = "")
    } else if (sqlite_avm_effectsize == "small") {
      b[i,2] = paste("$^{\\ast}$",b[i,2], sep = "")
    } else {

    }

    # U-test AVMD vs dr
    avmd_mutation <- sqlite_avmd$mutationScore

    p1 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      b[i,3] = paste("\\textbf{",b[i,3],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      b[i,3] = paste("\\textit{",b[i,3],"}", sep = "")
    } else {
    }

    if (sqlite_avmd_effectsize == "large") {
      b[i,3] = paste("$^{\\ast\\ast\\ast}$",b[i,3], sep = "")
    } else if (sqlite_avmd_effectsize == "medium") {
      b[i,3] = paste("$^{\\ast\\ast}$",b[i,3], sep = "")
    } else if (sqlite_avmd_effectsize == "small") {
      b[i,3] = paste("$^{\\ast}$",b[i,3], sep = "")
    } else {

    }

    # U-test for Random vs DR
    rand_mutation <- sqlite_rand$mutationScore

    p1 <- wilcox.test(dr_mutation, rand_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, rand_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      b[i,4] = paste("\\textbf{",b[i,4],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      b[i,4] = paste("\\textit{",b[i,4],"}", sep = "")
    } else {
    }

    if (sqlite_rand_effectsize == "large") {
      b[i,4] = paste("$^{\\ast\\ast\\ast}$",b[i,4], sep = "")
    } else if (sqlite_rand_effectsize == "medium") {
      b[i,4] = paste("$^{\\ast\\ast}$",b[i,4], sep = "")
    } else if (sqlite_rand_effectsize == "small") {
      b[i,4] = paste("$^{\\ast}$",b[i,4], sep = "")
    } else {

    }

    # A12 for HSQL
    hsql_avm_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_avm$mutationScore)$size
    hsql_avmd_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_avmd$mutationScore)$size
    hsql_rand_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_rand$mutationScore)$size
    hsql_dravm_effectsize <- dominoR::effectsize_accurate(hsql_dr$mutationScore, hsql_dravm$mutationScore)$size


    # U-test for DR vs AVM-R
    dr_mutation <- hsql_dr$mutationScore
    avmr_mutation <- hsql_avm$mutationScore

    p1 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmr_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      c[i,2] = paste("\\textbf{",c[i,2],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      c[i,2] = paste("\\textit{",c[i,2],"}", sep = "")
    } else {
    }

    if (hsql_avm_effectsize == "large") {
      c[i,2] = paste("$^{\\ast\\ast\\ast}$",c[i,2], sep = "")
    } else if (hsql_avm_effectsize == "medium") {
      c[i,2] = paste("$^{\\ast\\ast}$",c[i,2], sep = "")
    } else if (hsql_avm_effectsize == "small") {
      c[i,2] = paste("$^{\\ast}$",c[i,2], sep = "")
    } else {

    }

    # AVM-D vs DR U-test
    avmd_mutation <- hsql_avmd$mutationScore
    p1 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, avmd_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      c[i,3] = paste("\\textbf{",c[i,3],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      c[i,3] = paste("\\textit{",c[i,3],"}", sep = "")
    } else {
    }

    if (hsql_avmd_effectsize == "large") {
      c[i,3] = paste("$^{\\ast\\ast\\ast}$",c[i,3], sep = "")
    } else if (hsql_avmd_effectsize == "medium") {
      c[i,3] = paste("$^{\\ast\\ast}$",c[i,3], sep = "")
    } else if (hsql_avmd_effectsize == "small") {
      c[i,3] = paste("$^{\\ast}$",c[i,3], sep = "")
    } else {

    }

    # Random vs Dr U-test
    rand_mutation <- hsql_rand$mutationScore

    p1 <- wilcox.test(dr_mutation, rand_mutation, alternative = "greater", exact = FALSE)$p.value <= 0.01
    p2 <- wilcox.test(dr_mutation, rand_mutation, alternative = "less", exact = FALSE)$p.value <= 0.01

    if (p1 == TRUE & p2 == FALSE) {
      c[i,4] = paste("\\textbf{",c[i,4],"}", sep = "")
    } else if (p1 == FALSE & p2 == TRUE) {
      c[i,4] = paste("\\textit{",c[i,4],"}", sep = "")
    } else {
    }

    if (hsql_rand_effectsize == "large") {
      c[i,4] = paste("$^{\\ast\\ast\\ast}$",c[i,4], sep = "")
    } else if (hsql_rand_effectsize == "medium") {
      c[i,4] = paste("$^{\\ast\\ast}$",c[i,4], sep = "")
    } else if (hsql_rand_effectsize == "small") {
      c[i,4] = paste("$^{\\ast}$",c[i,4], sep = "")
    } else {

    }

    # For latex table
    a1[i,] <- paste("\\", a1[i,], sep = "")


  }
  a <- a[c(1,4,3,2)]
  b <- b[c(1,4,3,2)]
  c <- c[c(1,4,3,2)]

  # With HSQL
  d <- cbind(a1,c,a,b)
  # Without HSQL
  #d <- cbind(a1,a,b)

  if (rtrn == "tex") {
    return(print(xtable::xtable(d), include.rownames=FALSE ,sanitize.text.function = function(x){x}))
  } else {
    return(d)
  }
}

#' FUNCTION: ordering_mutants_per_schema
#'
#' It generates an ordered data frame of mutants (normal type) grouped by each run per schema and its mutation score.
#' @param d Data frame of mutants
#' @return A data frame of ordred mutants and grouped by runs and mutation score per run per schema
#' @importFrom magrittr %>%
#' @export
ordering_mutants_per_schema <- function(d) {

  # Only selecting normal mutants types
  d1 <- d %>% dplyr::filter(type == "NORMAL")
  dt <- NULL

  # Get each case study
  casestudy <- as.vector(dplyr::distinct(d1, schema))[[1]]
  # Get each DBMS
  dbs <- as.vector(dplyr::distinct(d1, dbms))[[1]]
  for (case in casestudy) {
    schema1 <- case
    for (db in dbs) {
      # Filter data
      filtered_data <- d1 %>% dplyr::filter(schema == schema1, dbms == db) %>% dplyr::group_by(identifier, dbms)
      # Select first schema to be grouped
      first_schema <- filtered_data[1,3]
      test <- NULL

      # Get each run for DRavm
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "dravm") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema")) %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dravm_minsitrust <- test %>% dplyr::filter(datagenerator == "dravm") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      dravm <- dravm_minsitrust
      dravm <- dravm %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))


      # Get each run for DR
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "directedRandom") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema")) %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dr_minsitrust <- test %>% dplyr::filter(datagenerator == "directedRandom") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      dr <- dr_minsitrust
      dr <- dr %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-R
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avs") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avs_minsitrust <- test %>% dplyr::filter(datagenerator == "avs") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      avm <- avs_minsitrust
      avm <- avm %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-D
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avsDefaults") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avsd_minsitrust <- test %>% dplyr::filter(datagenerator == "avsDefaults") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      avmd <- avsd_minsitrust
      avmd <- avmd %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for Random
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "random") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      ran_minsitrust <- test %>% dplyr::filter(datagenerator == "random") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      rand <- ran_minsitrust
      rand <- rand %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # save each run per DBMS
      if (db == "Postgres") {
        postgres_dr <- dr
        postgres_avm <- avm
        postgres_avmd <- avmd
        postgres_rand <- rand
        postgres_dravm <- dravm

      } else if (db == "SQLite") {
        sqlite_dr <- dr
        sqlite_avm <- avm
        sqlite_avmd <- avmd
        sqlite_rand <- rand
        sqlite_dravm <- dravm

      } else if (db == "HyperSQL") {
        hsql_dr <- dr
        hsql_avm <- avm
        hsql_avmd <- avmd
        hsql_rand <- rand
        hsql_dravm <- dravm

      }

    }

    # Arrange by Runs (numbers)
    postgres_dr <- dplyr::arrange(postgres_dr, number)
    postgres_avm <- dplyr::arrange(postgres_avm, number)
    postgres_avmd <- dplyr::arrange(postgres_avmd, number)
    postgres_rand <- dplyr::arrange(postgres_rand, number)
    postgres_dravm <- dplyr::arrange(postgres_dravm, number)

    sqlite_dr <- dplyr::arrange(sqlite_dr, number)
    sqlite_avm <- dplyr::arrange(sqlite_avm, number)
    sqlite_avmd <- dplyr::arrange(sqlite_avmd, number)
    sqlite_rand <- dplyr::arrange(sqlite_rand, number)
    sqlite_dravm <- dplyr::arrange(sqlite_dravm, number)

    hsql_dr <- dplyr::arrange(hsql_dr, number)
    hsql_avm <- dplyr::arrange(hsql_avm, number)
    hsql_avmd <- dplyr::arrange(hsql_avmd, number)
    hsql_rand <- dplyr::arrange(hsql_rand, number)
    hsql_dravm <- dplyr::arrange(hsql_dravm, number)


    # bind them together
    postgres <- rbind(postgres_dr, postgres_avm, postgres_avmd, postgres_rand,postgres_dravm)
    sqlite <- rbind(sqlite_dr, sqlite_avm, sqlite_avmd, sqlite_rand,sqlite_dravm)
    hsql <- rbind(hsql_dr, hsql_avm, hsql_avmd, hsql_rand,hsql_dravm)

    dt <- rbind(dt, postgres, sqlite, hsql)
  }
  return(dt)
}

#' FUNCTION: ordering_mutants_per_schema_others
#'
#' It generates an ordered data frame of mutants (normal type) grouped by each run per schema and its mutation score.
#' @param d Data frame of mutants
#' @return A data frame of ordred mutants and grouped by runs and mutation score per run per schema
#' @importFrom magrittr %>%
#' @export
ordering_mutants_per_schema_others <- function(d) {

  # Only selecting normal mutants types
  d1 <- d %>% dplyr::filter(type == "NORMAL")
  dt <- NULL

  # Get each case study
  casestudy <- as.vector(dplyr::distinct(d1, schema))[[1]]
  # Get each DBMS
  dbs <- as.vector(dplyr::distinct(d1, dbms))[[1]]
  for (case in casestudy) {
    schema1 <- case
    for (db in dbs) {
      # Filter data
      filtered_data <- d1 %>% dplyr::filter(schema == schema1, dbms == db) %>% dplyr::group_by(identifier, dbms)
      # Select first schema to be grouped
      first_schema <- filtered_data[1,3]
      test <- NULL

      # Get each run for DR
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "directedRandom") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema")) %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dr_minsitrust <- test %>% dplyr::filter(datagenerator == "directedRandom") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      dr <- dr_minsitrust
      dr <- dr %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-R
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avs") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avs_minsitrust <- test %>% dplyr::filter(datagenerator == "avs") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      avm <- avs_minsitrust
      avm <- avm %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-D
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avsDefaults") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avsd_minsitrust <- test %>% dplyr::filter(datagenerator == "avsDefaults") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      avmd <- avsd_minsitrust
      avmd <- avmd %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for Random
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "random") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      ran_minsitrust <- test %>% dplyr::filter(datagenerator == "random") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      rand <- ran_minsitrust
      rand <- rand %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # save each run per DBMS
      if (db == "Postgres") {
        postgres_dr <- dr
        postgres_avm <- avm
        postgres_avmd <- avmd
        postgres_rand <- rand
      } else if (db == "SQLite") {
        sqlite_dr <- dr
        sqlite_avm <- avm
        sqlite_avmd <- avmd
        sqlite_rand <- rand
      } else if (db == "HyperSQL") {
        hsql_dr <- dr
        hsql_avm <- avm
        hsql_avmd <- avmd
        hsql_rand <- rand
      }

    }

    # Arrange by Runs (numbers)
    postgres_dr <- dplyr::arrange(postgres_dr, number)
    postgres_avm <- dplyr::arrange(postgres_avm, number)
    postgres_avmd <- dplyr::arrange(postgres_avmd, number)
    postgres_rand <- dplyr::arrange(postgres_rand, number)

    sqlite_dr <- dplyr::arrange(sqlite_dr, number)
    sqlite_avm <- dplyr::arrange(sqlite_avm, number)
    sqlite_avmd <- dplyr::arrange(sqlite_avmd, number)
    sqlite_rand <- dplyr::arrange(sqlite_rand, number)

    hsql_dr <- dplyr::arrange(hsql_dr, number)
    hsql_avm <- dplyr::arrange(hsql_avm, number)
    hsql_avmd <- dplyr::arrange(hsql_avmd, number)
    hsql_rand <- dplyr::arrange(hsql_rand, number)


    # bind them together
    postgres <- rbind(postgres_dr, postgres_avm, postgres_avmd, postgres_rand)
    sqlite <- rbind(sqlite_dr, sqlite_avm, sqlite_avmd, sqlite_rand)
    hsql <- rbind(hsql_dr, hsql_avm, hsql_avmd, hsql_rand)

    dt <- rbind(dt, postgres, sqlite, hsql)
  }
  return(dt)
}

#' FUNCTION: ordering_mutants_per_schema_domino
#'
#' It generates an ordered data frame of mutants (normal type) grouped by each run per schema and its mutation score.
#' Only for domino AVM and Random
#' @param d Data frame of mutants
#' @return A data frame of ordred mutants and grouped by runs and mutation score per run per schema
#' @importFrom magrittr %>%
#' @export
ordering_mutants_per_schema_domino <- function(d) {

  # Only selecting normal mutants types
  d1 <- d %>% dplyr::filter(type == "NORMAL")
  dt <- NULL

  # Get each case study
  casestudy <- as.vector(dplyr::distinct(d1, schema))[[1]]
  # Get each DBMS
  dbs <- as.vector(dplyr::distinct(d1, dbms))[[1]]
  for (case in casestudy) {
    schema1 <- case
    for (db in dbs) {
      # Filter data
      filtered_data <- d1 %>% dplyr::filter(schema == schema1, dbms == db) %>% dplyr::group_by(identifier, dbms)
      # Select first schema to be grouped
      first_schema <- filtered_data[1,3]
      test <- NULL

      # Get each run for DRavm
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "dravm") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema")) %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dravm_minsitrust <- test %>% dplyr::filter(datagenerator == "dravm") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      dravm <- dravm_minsitrust
      dravm <- dravm %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))


      # Get each run for DR
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "directedRandom") %>% dplyr::select(identifier,dbms,schema) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema")) %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dr_minsitrust <- test %>% dplyr::filter(datagenerator == "directedRandom") %>% dplyr::group_by(identifier, dbms, datagenerator, number, schema) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))
      dr <- dr_minsitrust
      dr <- dr %>% dplyr::group_by(number, datagenerator, dbms, schema) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))



      # save each run per DBMS
      if (db == "Postgres") {
        postgres_dr <- dr
        postgres_dravm <- dravm

      } else if (db == "SQLite") {
        sqlite_dr <- dr
        sqlite_dravm <- dravm

      } else if (db == "HyperSQL") {
        hsql_dr <- dr
        hsql_dravm <- dravm

      }

    }

    # Arrange by Runs (numbers)
    postgres_dr <- dplyr::arrange(postgres_dr, number)
    postgres_dravm <- dplyr::arrange(postgres_dravm, number)

    sqlite_dr <- dplyr::arrange(sqlite_dr, number)
    sqlite_dravm <- dplyr::arrange(sqlite_dravm, number)

    hsql_dr <- dplyr::arrange(hsql_dr, number)
    hsql_dravm <- dplyr::arrange(hsql_dravm, number)


    # bind them together
    postgres <- rbind(postgres_dr,postgres_dravm)
    sqlite <- rbind(sqlite_dr,sqlite_dravm)
    hsql <- rbind(hsql_dr,hsql_dravm)

    dt <- rbind(dt, postgres, sqlite, hsql)
  }
  return(dt)
}

#' FUNCTION: ordering_mutants_per_operator
#'
#' It generates an ordered data frame of mutants (normal type) grouped by each run per operator and its mutation score.
#' @param d Data frame of mutants
#' @return A data frame of ordred mutants and grouped by runs and mutation score per run per operator
#' @importFrom magrittr %>%
#' @export
ordering_mutants_per_operator <- function(d) {
  # Only selecting normal mutants types
  d1 <- d %>% dplyr::filter(type == "NORMAL")
  dt <- NULL

  # get DBMSs
  dbs <- as.vector(distinct(d1, dbms))[[1]]
  # Get Operators
  operators <- as.vector(distinct(d1, operator))[[1]]
  for (selected_operator in operators) {
    for (db in dbs) {
      # Filter Data per operator
      filtered_data <- d1 %>% dplyr::filter(operator == selected_operator, schema != "iTrust", dbms == db) %>% dplyr::group_by(identifier, dbms)
      first_schema <- filtered_data[1,3]
      test <- NULL

      # Get each run for DR
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "directedRandom") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dr_minsitrust <- test %>% dplyr::filter(datagenerator == "directedRandom") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      dr_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "directedRandom", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(dr_itrust) > 0) {
        dr_itrust$number=1:nrow(dr_itrust)
      }
      dr <- rbind(dr_minsitrust, dr_itrust)
      dr <- dr %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-R
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avs") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avs_minsitrust <- test %>% dplyr::filter(datagenerator == "avs") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      avs_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "avs", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(avs_itrust) > 0) {
        avs_itrust$number=1:nrow(avs_itrust)
      }
      avm <- rbind(avs_minsitrust, avs_itrust)
      avm <- avm %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-D
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avsDefaults") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avsd_minsitrust <- test %>% dplyr::filter(datagenerator == "avsDefaults") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      avsd_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "avsDefaults", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(avsd_itrust) > 0) {
        avsd_itrust$number=1:nrow(avsd_itrust)
      }
      avmd <- rbind(avsd_minsitrust, avsd_itrust)
      avmd <- avmd %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for Random
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "random") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      ran_minsitrust <- test %>% dplyr::filter(datagenerator == "random") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      ran_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "random", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(ran_itrust) > 0) {
        ran_itrust$number=1:nrow(avsd_itrust)
      }
      rand <- rbind(ran_minsitrust, ran_itrust)
      rand <- rand %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # save each run per DBMS
      if (db == "Postgres") {
        postgres_dr <- dr
        postgres_avm <- avm
        postgres_avmd <- avmd
        postgres_rand <- rand
      } else if (db == "SQLite") {
        sqlite_dr <- dr
        sqlite_avm <- avm
        sqlite_avmd <- avmd
        sqlite_rand <- rand
      } else if (db == "HyperSQL") {
        hsql_dr <- dr
        hsql_avm <- avm
        hsql_avmd <- avmd
        hsql_rand <- rand
      }

    }
    # Arrange by Runs (numbers)
    postgres_dr <- dplyr::arrange(postgres_dr, number)
    postgres_avm <- dplyr::arrange(postgres_avm, number)
    postgres_avmd <- dplyr::arrange(postgres_avmd, number)
    postgres_rand <- dplyr::arrange(postgres_rand, number)
    sqlite_dr <- dplyr::arrange(sqlite_dr, number)
    sqlite_avm <- dplyr::arrange(sqlite_avm, number)
    sqlite_avmd <- dplyr::arrange(sqlite_avmd, number)
    sqlite_rand <- dplyr::arrange(sqlite_rand, number)
    hsql_dr <- dplyr::arrange(hsql_dr, number)
    hsql_avm <- dplyr::arrange(hsql_avm, number)
    hsql_avmd <- dplyr::arrange(hsql_avmd, number)
    hsql_rand <- dplyr::arrange(hsql_rand, number)

    # bind them together
    postgres <- rbind(postgres_dr, postgres_avm, postgres_avmd, postgres_rand)
    sqlite <- rbind(sqlite_dr, sqlite_avm, sqlite_avmd, sqlite_rand)
    hsql <- rbind(hsql_dr, hsql_avm, hsql_avmd, hsql_rand)

    dt <- rbind(dt, postgres, sqlite, hsql)
  }
  return(dt)
}

#' FUNCTION: ordering_mutants_per_operator
#'
#' It generates an ordered data frame of mutants (normal type) grouped by each run per operator and its mutation score.
#' @param d Data frame of mutants
#' @return A data frame of ordred mutants and grouped by runs and mutation score per run per operator
#' @importFrom magrittr %>%
#' @export
ordering_mutants_per_operator_domino <- function(d) {
  # Only selecting normal mutants types
  d1 <- d %>% dplyr::filter(type == "NORMAL")
  dt <- NULL

  # get DBMSs
  dbs <- as.vector(distinct(d1, dbms))[[1]]
  # Get Operators
  operators <- as.vector(distinct(d1, operator))[[1]]
  for (selected_operator in operators) {
    for (db in dbs) {
      # Filter Data per operator
      filtered_data <- d1 %>% dplyr::filter(operator == selected_operator, schema != "iTrust", dbms == db) %>% dplyr::group_by(identifier, dbms)
      first_schema <- filtered_data[1,3]
      test <- NULL

      # Get each run for DR
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "directedRandom") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dr_minsitrust <- test %>% dplyr::filter(datagenerator == "directedRandom") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      dr_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "directedRandom", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(dr_itrust) > 0) {
        dr_itrust$number=1:nrow(dr_itrust)
      }
      dr <- rbind(dr_minsitrust, dr_itrust)
      dr <- dr %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for DrAVM
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "dravm") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      dravm_minsitrust <- test %>% dplyr::filter(datagenerator == "dravm") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      dravm_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "dravm", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(dravm_itrust) > 0) {
        dravm_itrust$number=1:nrow(dravm_itrust)
      }
      dravm <- rbind(dravm_minsitrust, dravm_itrust)
      dravm <- dravm %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))


      # Get each run for AVM-R
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avs") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avs_minsitrust <- test %>% dplyr::filter(datagenerator == "avs") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      avs_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "avs", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(avs_itrust) > 0) {
        avs_itrust$number=1:nrow(avs_itrust)
      }
      avm <- rbind(avs_minsitrust, avs_itrust)
      avm <- avm %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for AVM-D
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "avsDefaults") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      avsd_minsitrust <- test %>% dplyr::filter(datagenerator == "avsDefaults") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      avsd_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "avsDefaults", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(avsd_itrust) > 0) {
        avsd_itrust$number=1:nrow(avsd_itrust)
      }
      avmd <- rbind(avsd_minsitrust, avsd_itrust)
      avmd <- avmd %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # Get each run for Random
      ids <- filtered_data %>% dplyr::filter(schema== first_schema[[1,1]], datagenerator == "random") %>% dplyr::select(identifier,dbms,schema,operator,type) %>% unique
      ids$number=1:nrow(ids)
      filtered_data %>% left_join(ids, by = c("identifier", "dbms", "schema", "operator", "type"))  %>% dplyr::mutate(number=as.numeric(ifelse(is.na(number),1,number))) %>% ungroup %>% dplyr::mutate(number = cummax(number)) -> test
      ran_minsitrust <- test %>% dplyr::filter(datagenerator == "random") %>% dplyr::group_by(identifier, dbms, datagenerator, number, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      ran_itrust <- d1 %>% dplyr::filter(schema == "iTrust", datagenerator == "random", type == "NORMAL", operator == selected_operator, dbms == db) %>% dplyr::group_by(identifier, dbms, datagenerator, operator) %>% dplyr::summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false")))# %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
      # Chekcing if there is any Itrust mutants
      if (nrow(ran_itrust) > 0) {
        ran_itrust$number=1:nrow(avsd_itrust)
      }
      rand <- rbind(ran_minsitrust, ran_itrust)
      rand <- rand %>% dplyr::group_by(number, datagenerator, dbms, operator) %>% dplyr::summarise(killed_mutants = sum(killed_mutants), total_mutants = sum(total_mutants)) %>% dplyr::mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

      # save each run per DBMS
      if (db == "Postgres") {
        postgres_dr <- dr
        postgres_avm <- avm
        postgres_avmd <- avmd
        postgres_rand <- rand
        postgres_dravm <- dravm
      } else if (db == "SQLite") {
        sqlite_dr <- dr
        sqlite_avm <- avm
        sqlite_avmd <- avmd
        sqlite_rand <- rand
        sqlite_dravm <- dravm
      } else if (db == "HyperSQL") {
        hsql_dr <- dr
        hsql_avm <- avm
        hsql_avmd <- avmd
        hsql_rand <- rand
        hsql_dravm <- dravm
      }

    }

    # Arrange by Runs (numbers)
    postgres_dr <- dplyr::arrange(postgres_dr, number)
    postgres_avm <- dplyr::arrange(postgres_avm, number)
    postgres_avmd <- dplyr::arrange(postgres_avmd, number)
    postgres_rand <- dplyr::arrange(postgres_rand, number)
    postgres_dravm <- dplyr::arrange(postgres_dravm, number)

    sqlite_dr <- dplyr::arrange(sqlite_dr, number)
    sqlite_avm <- dplyr::arrange(sqlite_avm, number)
    sqlite_avmd <- dplyr::arrange(sqlite_avmd, number)
    sqlite_rand <- dplyr::arrange(sqlite_rand, number)
    sqlite_dravm <- dplyr::arrange(sqlite_dravm, number)

    hsql_dr <- dplyr::arrange(hsql_dr, number)
    hsql_avm <- dplyr::arrange(hsql_avm, number)
    hsql_avmd <- dplyr::arrange(hsql_avmd, number)
    hsql_rand <- dplyr::arrange(hsql_rand, number)
    hsql_dravm <- dplyr::arrange(hsql_dravm, number)


    # bind them together
    postgres <- rbind(postgres_dr, postgres_avm, postgres_avmd, postgres_rand, postgres_dravm)
    sqlite <- rbind(sqlite_dr, sqlite_avm, sqlite_avmd, sqlite_rand, sqlite_dravm)
    hsql <- rbind(hsql_dr, hsql_avm, hsql_avmd, hsql_rand, hsql_dravm)

    dt <- rbind(dt, postgres, sqlite, hsql)
  }
  return(dt)
}

#' FUNCTION: dominoR::comparing_sig
#'
#'
#' @return latex of cell
#' @importFrom magrittr %>%
#' @export
comparing_sig <- function(sample1, sample2, effect, result) {
  # U-test AVMR vs DR
  p1 <- wilcox.test(sample1, sample2, alternative = "greater", exact = FALSE)$p.value <= 0.01
  p2 <- wilcox.test(sample1, sample2, alternative = "less", exact = FALSE)$p.value <= 0.01
  r <- result
  # Check one-sided test
  if (p1 == TRUE & p2 == FALSE) {
    r = paste("$\\APLdown$",result,"", sep = "")
  } else if (p1 == FALSE & p2 == TRUE) {
    r = paste("$\\APLup$",result,"", sep = "")
  } else {
  }

  # check effect
  if (effect == "large") {
    r = paste("$^{\\ast}$",r, sep = "")
  }
  # else if (effect == "medium") {
  #   r = paste("$^{\\ast\\ast}$",r, sep = "")
  # } else if (effect == "small") {
  #   r = paste("$^{\\ast}$",r, sep = "")
  # } else {
  #
  # }
  return(r)
}


#' FUNCTION: dominoR::comparing_sig_timing
#'
#'
#' @return latex of cell
#' @importFrom magrittr %>%
#' @export
comparing_sig_timing <- function(sample1, sample2, effect, result) {
  # U-test AVMR vs DR
  p1 <- wilcox.test(sample1, sample2, alternative = "greater", exact = FALSE)$p.value >= 0.01
  p2 <- wilcox.test(sample1, sample2, alternative = "less", exact = FALSE)$p.value >= 0.01
  r <- result
  # Check one-sided test
  if (p1 == TRUE & p2 == FALSE) {
    r = paste("$\\APLup$",result,"", sep = "")
  } else if (p1 == FALSE & p2 == TRUE) {
    r = paste("$\\APLdown$",result,"", sep = "")
  } else {
  }

  # check effect
  if (effect == "large") {
    r = paste("$^{\\ast}$",r, sep = "")
  }
  # else if (effect == "medium") {
  #   r = paste("$^{\\ast\\ast}$",r, sep = "")
  # } else if (effect == "small") {
  #   r = paste("$^{\\ast}$",r, sep = "")
  # } else {
  #
  # }

  return(r)

}
