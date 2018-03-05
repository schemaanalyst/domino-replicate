bootstrap <- function() {
  ana <- dominoReplicate::read_analysis()
  mut <- dominoReplicate::read_mutants()

  dominoReplicate::table_generator_coverage((ana %>% filter(datagenerator != "dravm")), m = "mean")
  dominoReplicate::table_generator_timing_nonRandom((ana %>% filter(datagenerator != "dravm")), m = "mean")
  dominoReplicate::table_generator_mutation_score_nonRandom(mut, m = "mean")
  dominoReplicate::domino_table_combined(ana, mut, mm = "mean")
}
