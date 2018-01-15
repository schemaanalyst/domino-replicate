bootstrap <- function() {
  ana <- ragtag::read_analysis()
  mut <- ragtag::read_mutants()

  ragtag::table_generator_coverage_others((ana %>% filter(datagenerator != "dravm")), m = "mean")
  ragtag::table_generator_timing_others_nonRand((ana %>% filter(datagenerator != "dravm")), m = "mean")
  ragtag::table_generator_mutation_score_others_nonRand(mut, m = "mean")
  ragtag::concentro_table_combaining(ana, mut, mm = "mean")
}
