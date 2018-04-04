#' FUNCTION: slides_timing_plot
#'
#' Timing Plot for ICST slides
#' @param analysis A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
slides_timing_plot <- function(analysis) {

  analysis %>% filter(dbms == "Postgres",
                      !datagenerator %in% c("dravm", "random"),
                      casestudy %in% c("BrowserCookies", "Flights", "iTrust")) %>%
    ggplot(aes(x = datagenerator, y = testgenerationtime / 1000, fill=datagenerator)) +
    geom_boxplot() +
    facet_grid(dbms ~ casestudy) +
    theme_bw() +
    theme(legend.position="none", text = element_text(size=18, family="Arial")) +
    scale_fill_manual(values=c("#999999", "#009E73", "#0072B2")) +
    labs(x = "Test Data Generators", y = "Test Generation Timing (Seconds)") +
    scale_x_discrete(breaks=c("avs", "avsDefaults", "directedRandom"),
                     labels=c("AVM-R", "AVM-D", "DOMINO"))

}

#' FUNCTION: slides_coverage_plot
#'
#' Coverage Plot for ICST slides
#' @param analysis A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
slides_coverage_plot <- function(analysis) {

  analysis %>% filter(dbms == "Postgres",
                      !datagenerator %in% c("dravm", "random"),
                      casestudy %in% c("BrowserCookies", "Flights", "iTrust")) %>%
    ggplot(aes(x = datagenerator, y = coverage, fill=datagenerator)) +
    geom_boxplot() +
    facet_grid(dbms ~ casestudy) +
    theme_bw() +
    theme(legend.position="none", text = element_text(size=18, family="Arial")) +
    scale_fill_manual(values=c("#999999", "#009E73", "#0072B2")) +
    labs(x = "Test Data Generators", y = "Coverage") +
    scale_x_discrete(breaks=c("avs", "avsDefaults", "directedRandom"),
                     labels=c("AVM-R", "AVM-D", "DOMINO"))

}

#' FUNCTION: slides_mutation_plot
#'
#' Mutation Plot for ICST slides
#' @param mutation A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
slides_mutation_plot <- function(mutation) {
  mutation2 <- dominoReplicate::ordering_mutants_per_schema_others(mutation)
  mutation2 %>% filter(dbms == "Postgres",
                       !datagenerator %in% c("dravm", "random"),
                       schema %in% c("BrowserCookies", "Flights", "iTrust")) %>%
    ggplot(aes(x = datagenerator, y = mutationScore, fill=datagenerator)) +
    geom_boxplot() +
    facet_grid(dbms ~ schema) +
    theme_bw() +
    theme(legend.position="none", text = element_text(size=18, family="Arial")) +
    scale_fill_manual(values=c("#999999", "#009E73", "#0072B2")) +
    labs(x = "Test Data Generators", y = "Mutation Score (%)") +
    scale_x_discrete(breaks=c("avs", "avsDefaults", "directedRandom"),
                     labels=c("AVM-R", "AVM-D", "DOMINO"))

}


#' FUNCTION: slides_timing_hybrid_plot
#'
#' Timing Plot for ICST slides
#' @param analysis A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
slides_timing_hybrid_plot <- function(analysis) {

  analysis %>% filter(dbms == "Postgres",
                      datagenerator %in% c("dravm", "directedRandom"),
                      casestudy %in% c("BrowserCookies", "Flights", "Products")) %>%
    ggplot(aes(x = datagenerator, y = testgenerationtime / 1000, fill=datagenerator)) +
    geom_boxplot() +
    facet_grid(dbms ~ casestudy) +
    theme_bw() +
    theme(legend.position="none", text = element_text(size=18, family="Arial")) +
    scale_fill_manual(values=c("#0072B2", "#D55E00")) +
    labs(x = "Test Data Generators", y = "Test Generation Timing (Seconds)") +
    scale_x_discrete(breaks=c("directedRandom", "dravm"),
                     labels=c("DOMINO", "DOMINO-AVM"))

}

#' FUNCTION: slides_coverage_hybrid_plot
#'
#' Coverage for DominoAVM Plot for ICST slides
#' @param analysis A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
slides_coverage_hybrid_plot <- function(analysis) {

  analysis %>% filter(dbms == "Postgres",
                      datagenerator %in% c("dravm", "directedRandom"),
                      casestudy %in% c("BrowserCookies", "Products", "Flights")) %>%
    ggplot(aes(x = datagenerator, y = coverage, fill=datagenerator)) +
    geom_boxplot() +
    facet_grid(dbms ~ casestudy) +
    theme_bw() +
    theme(legend.position="none", text = element_text(size=18, family="Arial")) +
    scale_fill_manual(values=c("#0072B2", "#D55E00")) +
    labs(x = "Test Data Generators", y = "Coverage") +
    scale_x_discrete(breaks=c("directedRandom", "dravm"),
                     labels=c("Domino", "Domino-AVM"))

}

#' FUNCTION: slides_mutation_hybrid_plot
#'
#' Mutation Plot for ICST slides
#' @param mutation A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
slides_mutation_hybrid_plot <- function(mutation) {
  mutation2 <- dominoReplicate::ordering_mutants_per_schema_domino(mutation)
  mutation2 %>% filter(dbms == "Postgres",
                       datagenerator %in% c("dravm", "directedRandom"),
                       schema %in% c("BrowserCookies", "Flights", "Products")) %>%
    ggplot(aes(x = datagenerator, y = mutationScore, fill=datagenerator)) +
    geom_boxplot() +
    facet_grid(dbms ~ schema) +
    theme_bw() +
    theme(legend.position="none", text = element_text(size=18, family="Arial")) +
    scale_fill_manual(values=c("#0072B2", "#D55E00")) +
    labs(x = "Test Data Generators", y = "Mutation Score (%)") +
    scale_x_discrete(breaks=c("directedRandom", "dravm"),
                     labels=c("DOMINO", "DOMINO-AVM"))

}

#' FUNCTION: demo_timing_plot
#'
#' Timing Plot for ICST slides
#' @param analysis A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
demo_timing_plot <- function(analysis) {

  analysis %>% filter(dbms == "Postgres",
                      datagenerator %in% c("directedRandom", "avs"),
                      casestudy %in% c("iTrust")) %>%
    ggplot(aes(x = datagenerator, y = testgenerationtime / 1000)) +
    geom_boxplot() +
    theme(legend.position="none", text = element_text(size=12)) +
    labs(x = "Test Data Generators", y = "Test Generation Time (Seconds)") +
    scale_x_discrete(breaks=c("avs", "directedRandom"),
                     labels=c("AVM-R", "DOMINO"))

  #ggsave(filename = "figure11.pdf",
  #       plot = a,
  #       path = "/home/abdullah/Documents/papers/demoAndArtifact/demo/figures",
  #       device = "pdf", , width = 4, height = 3, units = "in")

}


#' FUNCTION: demo_timing_plot
#'
#' Timing Plot for ICST slides
#' @param analysis A data frame of read_analysis
#' @return A plot
#' @importFrom magrittr %>%
#' @export
demo_mutation_plot <- function(mutation) {
  mutation2 <- dominoReplicate::ordering_mutants_per_schema_others(mutation)
  mutation2 %>% filter(dbms == "Postgres",
                       datagenerator %in% c("directedRandom", "avs"),
                       schema %in% c("iTrust")) %>%
    ggplot(aes(x = datagenerator, y = mutationScore)) +
    geom_boxplot() +
    theme(legend.position="none", text = element_text(size=12)) +
    labs(x = "Test Data Generators", y = "Mutation Score (%)") +
    scale_x_discrete(breaks=c("avs", "directedRandom"),
                     labels=c("AVM-R", "DOMINO"))

  #ggsave(filename = "figure11.pdf",
  #       plot = a,
  #       path = "/home/abdullah/Documents/papers/demoAndArtifact/demo/figures",
  #       device = "pdf", , width = 4, height = 3, units = "in")

}
