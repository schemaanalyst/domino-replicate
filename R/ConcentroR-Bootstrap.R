bootstrap <- function() {
  ana <- dominoR::read_analysis()
  mut <- dominoR::read_mutants()

  dominoR::table_generator_coverage_others((ana %>% filter(datagenerator != "dravm")), m = "mean")
  dominoR::table_generator_timing_others_nonRand((ana %>% filter(datagenerator != "dravm")), m = "mean")
  dominoR::table_generator_mutation_score_others_nonRand(mut, m = "mean")
  dominoR::concentro_table_combaining(ana, mut, mm = "mean")
}
