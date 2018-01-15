context("test-mutation-scores.R")
### Context: test case for checking mutant operators sample size results

test_that("Sample Size of mutants per schema", {
  setwd('../../')
  # Gettting Mutants
  mutants <- ragtag::read_mutants()
  # Filtering mutants and getting precentage kills SQLite
  mutants <- mutants %>% filter(type == "NORMAL", dbms == "SQLite") %>% select(identifier, dbms, schema, datagenerator, killed) %>% group_by(identifier, dbms, schema, datagenerator) %>% summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false"))) %>% mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
  # getting mutants per generator

  dr <- mutants %>% filter(datagenerator == "directedRandom")
  avmr <- mutants %>% filter(datagenerator == "avs")
  avmd <- mutants %>% filter(datagenerator == "avsDefaults")
  rand <- mutants %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(schema == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(schema == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(schema == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(schema == "StudentResidence")

  # Check 30 runs
  expect_equal(nrow(dr_StudentResidence), 30)
  expect_equal(nrow(avmr_StudentResidence), 30)
  expect_equal(nrow(avmd_StudentResidence), 30)
  expect_equal(nrow(rand_StudentResidence), 30)
})

### Context: test case for checking mutant scores average results
context("mutant-scores-mean-results")

test_that("Averages for mutation socres", {
  setwd('../../')
  # Gettting Mutants
  mutants <- ragtag::read_mutants()
  # Filtering mutants and getting precentage kills SQLite
  mutants <- mutants %>% filter(type == "NORMAL", dbms == "SQLite") %>% select(identifier, dbms, schema, datagenerator, killed) %>% group_by(identifier, dbms, schema, datagenerator) %>% summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false"))) %>% mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
  # getting mutants per generator

  dr <- mutants %>% filter(datagenerator == "directedRandom")
  avmr <- mutants %>% filter(datagenerator == "avs")
  avmd <- mutants %>% filter(datagenerator == "avsDefaults")
  rand <- mutants %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(schema == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(schema == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(schema == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(schema == "StudentResidence")

  # check averages

  averages <- round(sum(dr_StudentResidence$killed_mutants) / sum(dr_StudentResidence$total_mutants) * 100, 2)
  expect_equal(averages, 95.73)

  averages <- round(sum(avmr_StudentResidence$killed_mutants) / sum(avmr_StudentResidence$total_mutants) * 100, 2)
  expect_equal(averages, 96.58)

  averages <- round(sum(avmd_StudentResidence$killed_mutants) / sum(avmd_StudentResidence$total_mutants) * 100, 2)
  expect_equal(averages, 87.18)

  averages <- round(sum(rand_StudentResidence$killed_mutants) / sum(rand_StudentResidence$total_mutants) * 100, 2)
  expect_equal(averages, 77.69)
})

### Context: test case for checking mutant scores U-Test results
context("mutant-scores-u-test-results")

test_that("U-test results for mutation socres", {
  setwd('../../')
  # Gettting Mutants
  mutants <- ragtag::read_mutants()
  # Filtering mutants and getting precentage kills SQLite
  mutants <- mutants %>% filter(type == "NORMAL", dbms == "SQLite") %>% select(identifier, dbms, schema, datagenerator, killed) %>% group_by(identifier, dbms, schema, datagenerator) %>% summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false"))) %>% mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
  # getting mutants per generator

  dr <- mutants %>% filter(datagenerator == "directedRandom")
  avmr <- mutants %>% filter(datagenerator == "avs")
  avmd <- mutants %>% filter(datagenerator == "avsDefaults")
  rand <- mutants %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(schema == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(schema == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(schema == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(schema == "StudentResidence")

  # Testing Siginificant results to the table
  p <- wilcox.test(dr_StudentResidence$mutationScore, avmr_StudentResidence$mutationScore, exact = FALSE)$p.value <= 0.05
  expect_true(p)

  p <- wilcox.test(dr_StudentResidence$mutationScore, avmd_StudentResidence$mutationScore, exact = FALSE)$p.value <= 0.05
  expect_true(p)

  p <- wilcox.test(dr_StudentResidence$mutationScore, rand_StudentResidence$mutationScore, exact = FALSE)$p.value <= 0.05
  expect_true(p)
})

### Context: test case for checking mutant scores Effect Size results
context("mutant-scores-effect-size-results")

test_that("Effect size reuslts for mutation scores", {
  setwd('../../')
  # Gettting Mutants
  mutants <- ragtag::read_mutants()
  # Filtering mutants and getting precentage kills SQLite
  mutants <- mutants %>% filter(type == "NORMAL", dbms == "SQLite") %>% select(identifier, dbms, schema, datagenerator, killed) %>% group_by(identifier, dbms, schema, datagenerator) %>% summarise(killed_mutants = sum(killed == "true"), total_mutants = (sum(killed == "true") + sum(killed == "false"))) %>% mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))
  # getting mutants per generator

  dr <- mutants %>% filter(datagenerator == "directedRandom")
  avmr <- mutants %>% filter(datagenerator == "avs")
  avmd <- mutants %>% filter(datagenerator == "avsDefaults")
  rand <- mutants %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(schema == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(schema == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(schema == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(schema == "StudentResidence")

  # Testing Effect size

  avmr_effectsize <- ragtag::effectsize_accurate(dr_StudentResidence$mutationScore, avmr_StudentResidence$mutationScore)$size
  avmd_effectsize <- ragtag::effectsize_accurate(dr_StudentResidence$mutationScore, avmd_StudentResidence$mutationScore)$size
  rand_effectsize <- ragtag::effectsize_accurate(dr_StudentResidence$mutationScore, rand_StudentResidence$mutationScore)$size

  expect_equal(avmr_effectsize, "small")

  expect_equal(avmd_effectsize, "large")

  expect_equal(rand_effectsize, "large")

})
