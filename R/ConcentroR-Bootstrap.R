bootstrap <- function() {
  ana <- dominoR::read_analysis()
  mut <- dominoR::read_mutants()

  dominoR::table_generator_coverage((ana %>% filter(datagenerator != "dravm")), m = "mean")
  dominoR::table_generator_timing_nonRandom((ana %>% filter(datagenerator != "dravm")), m = "mean")
  dominoR::table_generator_mutation_score_nonRandom(mut, m = "mean")
  dominoR::domino_table_combined(ana, mut, mm = "mean")
}
