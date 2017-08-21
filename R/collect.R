#' FUNCTION: read_analysis
#'
#' Read the data file that contains the "analysis" data. This is the data containing all test generation times, coverages,
#' evaluations etc. It is refered in SchemaAnalyst github repo as 'newmutationanalysis.dat'. And it allow us to compare
#' test generation timing and coverages results.
#' @return A Data Frame of analysis
#' @importFrom magrittr %>%
#' @export
read_analysis <- function() {

  # collect all AVM-D
  rnd.avmd <- system.file("extdata", "mutation-analysis.dat",
                               package="ragtag")

  # collect all AVM-R
  HyperSQL.avmr.itrust <- system.file("extdata", "30-HyperSQL-avs-itrust-mutationanalysis.dat",
                               package="ragtag")
  SQLite.avmr.itrust <- system.file("extdata", "30-SQLite-avs-itrust-mutationanalysis.dat",
                             package="ragtag")
  Postgres.avmr.itrust <- system.file("extdata", "30-Postgres-avs-itrust-mutationanalysis.dat",
                               package="ragtag")

  HyperSQL.avmr.minusitrust <- system.file("extdata", "30-HyperSQL-avs-minusitrust-mutationanalysis.dat",
                                      package="ragtag")
  SQLite.avmr.minusitrust <- system.file("extdata", "30-SQLite-avs-minusitrust-mutationanalysis.dat",
                                    package="ragtag")
  Postgres.avmr.minusitrust <- system.file("extdata", "30-Postgres-avs-minusitrust-mutationanalysis.dat",
                                      package="ragtag")

  # collect all DrRan

  HyperSQL.directedRandom.itrust <- system.file("extdata", "30-HyperSQL-directedRandom-itrust-mutationanalysis.dat",
                                      package="ragtag")
  SQLite.directedRandom.itrust <- system.file("extdata", "30-SQLite-directedRandom-itrust-mutationanalysis.dat",
                                    package="ragtag")
  Postgres.directedRandom.itrust <- system.file("extdata", "30-Postgres-directedRandom-itrust-mutationanalysis.dat",
                                      package="ragtag")

  HyperSQL.directedRandom.minusitrust <- system.file("extdata", "30-HyperSQL-directedRandom-minusitrust-mutationanalysis.dat",
                                           package="ragtag")
  SQLite.directedRandom.minusitrust <- system.file("extdata", "30-SQLite-directedRandom-minusitrust-mutationanalysis.dat",
                                         package="ragtag")
  Postgres.directedRandom.minusitrust <- system.file("extdata", "30-Postgres-directedRandom-minusitrust-mutationanalysis.dat",
                                           package="ragtag")

  # collect all DrAVM

  #HyperSQL.dravm.itrust <- system.file("extdata", "30-HyperSQL-dravm-itrust-mutationanalysis.dat",
  #                                    package="ragtag")
  SQLite.dravm.itrust <- system.file("extdata", "30-SQLite-dravm-itrust-mutationanalysis.dat",
                                    package="ragtag")
  #Postgres.dravm.itrust <- system.file("extdata", "30-Postgres-dravm-itrust-mutationanalysis.dat",
  #                                    package="ragtag")

  HyperSQL.dravm.minusitrust <- system.file("extdata", "30-HyperSQL-concentroAVS-minusitrust-mutationanalysis.dat",
                                           package="ragtag")
  SQLite.dravm.minusitrust <- system.file("extdata", "30-SQLite-dravm-minusitrust-mutationanalysis.dat",
                                         package="ragtag")
  Postgres.dravm.minusitrust <- system.file("extdata", "30-Postgres-dravm-minusitrust-mutationanalysis.dat",
                                           package="ragtag")

  #f <- system.file("extdata", "analysis.csv", package="ragtag")
  # Read all :)
  d1 <- readr::read_csv(rnd.avmd)


  d2 <- readr::read_csv(HyperSQL.avmr.minusitrust)
  d3 <- readr::read_csv(SQLite.avmr.minusitrust)
  d4 <- readr::read_csv(Postgres.avmr.minusitrust)

  d5 <- readr::read_csv(HyperSQL.avmr.itrust)
  d6 <- readr::read_csv(SQLite.avmr.itrust)
  d7 <- readr::read_csv(Postgres.avmr.itrust)


  d8 <- readr::read_csv(HyperSQL.directedRandom.minusitrust)
  d9 <- readr::read_csv(SQLite.directedRandom.minusitrust)
  d10 <- readr::read_csv(Postgres.directedRandom.minusitrust)

  d11 <- readr::read_csv(HyperSQL.directedRandom.itrust)
  d12 <- readr::read_csv(SQLite.directedRandom.itrust)
  d13 <- readr::read_csv(Postgres.directedRandom.itrust)

  d17 <- readr::read_csv(HyperSQL.dravm.minusitrust)
  d18 <- readr::read_csv(SQLite.dravm.minusitrust)
  d19 <- readr::read_csv(Postgres.dravm.minusitrust)

  #d14 <- readr::read_csv(HyperSQL.dravm.itrust)
  d15 <- readr::read_csv(SQLite.dravm.itrust)
  #d16 <- readr::read_csv(Postgres.dravm.itrust)

  #allFrames <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19)

  allFrames <- rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d15,d17,d18,d19)

  allFrames <- allFrames %>% dplyr::mutate(casestudy = as.character(gsub("parsedcasestudy.","",casestudy)))
  allFrames$casestudy <- gsub("IsoFlav_R2Repaired", "IsoFlav_R2", allFrames$casestudy)

  allFrames$datagenerator <- gsub("concentroAVS", "dravm", allFrames$datagenerator)

  return(dplyr::tbl_df(allFrames))
}

#' FUNCTION: read_mutants
#'
#' Read the data file that contains the "time-constrained mutation" data. This is the data file containts all mutants,
#' killed or alive. It is refered in SchemaAnalyst github repo as 'mutanttiming.dat'.
#' This file is useful if you are interested in looking at individual mutants.
#' This file contains seven attributes: identifier, dbms, schema, operator, type, killed, time.
#' It allow us to compare mutation scores.
#' @return A Data Frame of mutants
#' @importFrom magrittr %>%
#' @export
read_mutants <- function() {
  # collect all AVM-D
  hypersql.avmdefaults <- system.file("extdata", "hypersql-avmdefaults.dat",
                                      package="ragtag")
  sqlite.avmdefaults <- system.file("extdata", "sqlite-avmdefaults.dat",
                                    package="ragtag")
  postgres.avmdefaults <- system.file("extdata", "postgres-avmdefaults.dat",
                                      package="ragtag")

  hypersql.random <- system.file("extdata", "hypersql-random.dat",
                                 package="ragtag")
  sqlite.random <- system.file("extdata", "sqlite-random.dat",
                               package="ragtag")
  postgres.random <- system.file("extdata", "postgres-random.dat",
                                 package="ragtag")

  # Read data
  d21 <- readr::read_csv(hypersql.avmdefaults)
  d22 <- readr::read_csv(sqlite.avmdefaults)
  d23 <- readr::read_csv(postgres.avmdefaults)
  d24 <- readr::read_csv(hypersql.random)
  d25 <- readr::read_csv(sqlite.random)
  d26 <- readr::read_csv(postgres.random)

  # Adding generator name
  namevector <- c("datagenerator")

  d21[,namevector] <- "avsDefaults"
  d22[,namevector] <- "avsDefaults"
  d23[,namevector] <- "avsDefaults"
  d24[,namevector] <- "random"
  d25[,namevector] <- "random"
  d26[,namevector] <- "random"

  # collect all AVM-R
  HyperSQL.avmr.itrust <- system.file("extdata", "30-HyperSQL-avs-itrust-mutanttiming.dat",
                                      package="ragtag")
  SQLite.avmr.itrust <- system.file("extdata", "30-SQLite-avs-itrust-mutanttiming.dat",
                                    package="ragtag")
  Postgres.avmr.itrust <- system.file("extdata", "30-Postgres-avs-itrust-mutanttiming.dat",
                                      package="ragtag")

  HyperSQL.avmr.minusitrust <- system.file("extdata", "30-HyperSQL-avs-minusitrust-mutanttiming.dat",
                                           package="ragtag")
  SQLite.avmr.minusitrust <- system.file("extdata", "30-SQLite-avs-minusitrust-mutanttiming.dat",
                                         package="ragtag")
  Postgres.avmr.minusitrust <- system.file("extdata", "30-Postgres-avs-minusitrust-mutanttiming.dat",
                                           package="ragtag")

  # Read data
  d2 <- readr::read_csv(HyperSQL.avmr.minusitrust)
  d3 <- readr::read_csv(SQLite.avmr.minusitrust)
  d4 <- readr::read_csv(Postgres.avmr.minusitrust)

  d5 <- readr::read_csv(HyperSQL.avmr.itrust)
  d6 <- readr::read_csv(SQLite.avmr.itrust)
  d7 <- readr::read_csv(Postgres.avmr.itrust)

  # Adding generator name
  d2[,namevector] <- "avs"
  d3[,namevector] <- "avs"
  d4[,namevector] <- "avs"
  d5[,namevector] <- "avs"
  d6[,namevector] <- "avs"
  d7[,namevector] <- "avs"

  # collect all DrRan

  HyperSQL.directedRandom.itrust <- system.file("extdata", "30-HyperSQL-directedRandom-itrust-mutanttiming.dat",
                                                package="ragtag")
  SQLite.directedRandom.itrust <- system.file("extdata", "30-SQLite-directedRandom-itrust-mutanttiming.dat",
                                              package="ragtag")
  Postgres.directedRandom.itrust <- system.file("extdata", "30-Postgres-directedRandom-itrust-mutanttiming.dat",
                                                package="ragtag")

  HyperSQL.directedRandom.minusitrust <- system.file("extdata", "30-HyperSQL-directedRandom-minusitrust-mutanttiming.dat",
                                                     package="ragtag")
  SQLite.directedRandom.minusitrust <- system.file("extdata", "30-SQLite-directedRandom-minusitrust-mutanttiming.dat",
                                                   package="ragtag")
  Postgres.directedRandom.minusitrust <- system.file("extdata", "30-Postgres-directedRandom-minusitrust-mutanttiming.dat",
                                                     package="ragtag")
  # Read data
  d8 <- readr::read_csv(HyperSQL.directedRandom.minusitrust)
  d9 <- readr::read_csv(SQLite.directedRandom.minusitrust)
  d10 <- readr::read_csv(Postgres.directedRandom.minusitrust)

  d11 <- readr::read_csv(HyperSQL.directedRandom.itrust)
  d12 <- readr::read_csv(SQLite.directedRandom.itrust)
  d13 <- readr::read_csv(Postgres.directedRandom.itrust)

  # Adding generator name
  d8[,namevector] <- "directedRandom"
  d9[,namevector] <- "directedRandom"
  d10[,namevector] <- "directedRandom"
  d11[,namevector] <- "directedRandom"
  d12[,namevector] <- "directedRandom"
  d13[,namevector] <- "directedRandom"

  # collect all DrAVM

  #HyperSQL.dravm.itrust <- system.file("extdata", "30-HyperSQL-dravm-itrust-mutanttiming.dat",
  #                                    package="ragtag")
  SQLite.dravm.itrust <- system.file("extdata", "30-SQLite-dravm-itrust-mutanttiming.dat",
                                     package="ragtag")
  #Postgres.dravm.itrust <- system.file("extdata", "30-Postgres-dravm-itrust-mutanttiming.dat",
  #                                    package="ragtag")

  HyperSQL.dravm.minusitrust <- system.file("extdata", "30-HyperSQL-concentroAVS-minusitrust-mutanttiming.dat",
                                            package="ragtag")
  SQLite.dravm.minusitrust <- system.file("extdata", "30-SQLite-dravm-minusitrust-mutanttiming.dat",
                                          package="ragtag")
  Postgres.dravm.minusitrust <- system.file("extdata", "30-Postgres-dravm-minusitrust-mutanttiming.dat",
                                            package="ragtag")

  d17 <- readr::read_csv(HyperSQL.dravm.minusitrust)
  d18 <- readr::read_csv(SQLite.dravm.minusitrust)
  d19 <- readr::read_csv(Postgres.dravm.minusitrust)

  #d14 <- readr::read_csv(HyperSQL.dravm.itrust)
  d15 <- readr::read_csv(SQLite.dravm.itrust)
  #d16 <- readr::read_csv(Postgres.dravm.itrust)

  # Adding generator name
  #d14[,namevector] <- "dravm"
  d15[,namevector] <- "dravm"
  #d16[,namevector] <- "dravm"
  d17[,namevector] <- "dravm"
  d18[,namevector] <- "dravm"
  d19[,namevector] <- "dravm"

  #allFrames <- rbind(d21,d22,d23,d24,d25,d26,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19)

  allFrames <- rbind(d21,d22,d23,d24,d25,d26,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d17,d18,d19,d15)
  allFrames$datagenerator <- gsub("concentroAVS", "dravm", allFrames$datagenerator)


  return(dplyr::tbl_df(allFrames))

}
