context("test-testgenerationtiming.R")

### Context: test case for checking testgenerations sample size
test_that("Sample Size of test generation timings", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoR::read_analysis()
  # getting mutants per generator
  analysis <- analysis %>% filter(dbms == "HyperSQL")
  dr <- analysis %>% filter(datagenerator == "directedRandom")
  avmr <- analysis %>% filter(datagenerator == "avs")
  avmd <- analysis %>% filter(datagenerator == "avsDefaults")
  rand <- analysis %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(casestudy == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(casestudy == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(casestudy == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(casestudy == "StudentResidence")

  # Check 30 runs
  expect_equal(nrow(dr_StudentResidence), 30)
  expect_equal(nrow(avmr_StudentResidence), 30)
  expect_equal(nrow(avmd_StudentResidence), 30)
  expect_equal(nrow(rand_StudentResidence), 30)
})

### Context: test case for checking testgenerations medians
context("testgenerationtiming-median-results")

test_that("Medians for test generation timings resutls", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoR::read_analysis()
  # getting mutants per generator
  analysis <- analysis %>% filter(dbms == "HyperSQL")
  dr <- analysis %>% filter(datagenerator == "directedRandom")
  avmr <- analysis %>% filter(datagenerator == "avs")
  avmd <- analysis %>% filter(datagenerator == "avsDefaults")
  rand <- analysis %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(casestudy == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(casestudy == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(casestudy == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(casestudy == "StudentResidence")

  # check testgenerationtime median

  expect_equal(round(median(dr_StudentResidence$testgenerationtime) / 1000, 2), 0.57)

  expect_equal(round(median(avmr_StudentResidence$testgenerationtime) / 1000, 2), 0.96)

  expect_equal(round(median(avmd_StudentResidence$testgenerationtime) / 1000, 2), 0.78)

  expect_equal(round(median(rand_StudentResidence$testgenerationtime) / 1000, 2), 21.01)
})

### Context: test case for checking testgenerations U-test results
context("testgenerationtiming-u-test-results")

test_that("U-Test resutls for test generation timings", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoR::read_analysis()
  # getting mutants per generator
  analysis <- analysis %>% filter(dbms == "HyperSQL")
  dr <- analysis %>% filter(datagenerator == "directedRandom")
  avmr <- analysis %>% filter(datagenerator == "avs")
  avmd <- analysis %>% filter(datagenerator == "avsDefaults")
  rand <- analysis %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(casestudy == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(casestudy == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(casestudy == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(casestudy == "StudentResidence")

  # Testing Siginificant results to the table
  p <- wilcox.test(dr_StudentResidence$testgenerationtime, avmr_StudentResidence$testgenerationtime, exact = FALSE)$p.value <= 0.05
  # Because of ties
  expect_true(p)

  p <- wilcox.test(dr_StudentResidence$testgenerationtime, avmd_StudentResidence$testgenerationtime, exact = FALSE)$p.value <= 0.05
  expect_true(p)

  p <- wilcox.test(dr_StudentResidence$testgenerationtime, rand_StudentResidence$testgenerationtime, exact = FALSE)$p.value <= 0.05
  expect_true(p)
})

### Context: test case for checking testgenerations effect size results
context("testgenerationtiming-effect-size-results")

test_that("Effect size results for test generation timings", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoR::read_analysis()
  # getting mutants per generator
  analysis <- analysis %>% filter(dbms == "HyperSQL")
  dr <- analysis %>% filter(datagenerator == "directedRandom")
  avmr <- analysis %>% filter(datagenerator == "avs")
  avmd <- analysis %>% filter(datagenerator == "avsDefaults")
  rand <- analysis %>% filter(datagenerator == "random")


  # Checking StudentResidence schema
  dr_StudentResidence <- dr %>% filter(casestudy == "StudentResidence")
  avmr_StudentResidence <- avmr %>% filter(casestudy == "StudentResidence")
  avmd_StudentResidence <- avmd %>% filter(casestudy == "StudentResidence")
  rand_StudentResidence <- rand %>% filter(casestudy == "StudentResidence")

  # Transformation of test generation timing

  #dr_StudentResidence <- dominoR::transform_execution_times_for_threshold(dr_StudentResidence, 1000)
  #avmr_StudentResidence <- dominoR::transform_execution_times_for_threshold(avmr_StudentResidence, 1000)
  #avmd_StudentResidence <- dominoR::transform_execution_times_for_threshold(avmd_StudentResidence, 1000)
  #rand_StudentResidence <- dominoR::transform_execution_times_for_threshold(rand_StudentResidence, 1000)

  # Testing Effect size

  avmr_effectsize <- dominoR::effectsize_accurate(dr_StudentResidence$testgenerationtime, avmr_StudentResidence$testgenerationtime)$size
  avmd_effectsize <- dominoR::effectsize_accurate(dr_StudentResidence$testgenerationtime, avmd_StudentResidence$testgenerationtime)$size
  rand_effectsize <- dominoR::effectsize_accurate(dr_StudentResidence$testgenerationtime, rand_StudentResidence$testgenerationtime)$size

  expect_equal(avmr_effectsize, "large")

  expect_equal(avmd_effectsize, "large")

  expect_equal(rand_effectsize, "large")
})

