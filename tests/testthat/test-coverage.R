context("test-coverage.R")

test_that("Sample size of coverages for StudentResidence schema", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoReplicate::read_analysis()
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


### Context: test case for checking coverage median resutls
context("coverages-median-resutlts")

test_that("Testing Coverage Table Median resutls", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoReplicate::read_analysis()
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

  # check coverage median

  expect_equal(round(median(dr_StudentResidence$coverage), 2), 100.00)

  expect_equal(round(median(avmr_StudentResidence$coverage), 2), 100.00)

  expect_equal(round(median(avmd_StudentResidence$coverage), 2), 100.00)

  expect_equal(round(median(rand_StudentResidence$coverage), 2), 60.00)
})


### Context: test case for checking coverage effect size for one schema per data generator
context("coverages-effect-size")

test_that("Generating effect size for coverage", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoReplicate::read_analysis()
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

  # Testing Effect size

  avmr_effectsize <- dominoReplicate::effectsize_accurate(dr_StudentResidence$coverage, avmr_StudentResidence$coverage)$size
  avmd_effectsize <- dominoReplicate::effectsize_accurate(dr_StudentResidence$coverage, avmd_StudentResidence$coverage)$size
  rand_effectsize <- dominoReplicate::effectsize_accurate(dr_StudentResidence$coverage, rand_StudentResidence$coverage)$size

  expect_equal(avmr_effectsize, "none")

  expect_equal(avmd_effectsize, "none")

  expect_equal(rand_effectsize, "large")
})

### Context: test case for checking coverage U-Test results for one schema per data generator
context("coverages-u-test")

test_that("Generate U-Test results for StudentResidence schema per generator", {
  setwd('../../')
  # Gettting Mutants
  analysis <- dominoReplicate::read_analysis()
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
  p <- wilcox.test(dr_StudentResidence$coverage, avmr_StudentResidence$coverage, exact = FALSE)$p.value <= 0.05
  expect_false(isTRUE(p))

  p <- wilcox.test(dr_StudentResidence$coverage, avmd_StudentResidence$coverage, exact = FALSE)$p.value <= 0.05
  expect_false(isTRUE(p))

  p <- wilcox.test(dr_StudentResidence$coverage, rand_StudentResidence$coverage, exact = FALSE)$p.value <= 0.05
  expect_true(p)
})
