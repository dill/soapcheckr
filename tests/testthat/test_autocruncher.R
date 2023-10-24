# testing script to check autocruncher

# library(testthat)
#
# source("soap_check.R")


context("Autocruncher tests")

test_that("Sissabagama Lake gird checks",{
  # not a list
  expect_warning(
    crunch_ind <- autocruncher(bnd = sissabagama_lake_ls,
                               knots = sissabagama_lake_grid_knots
    )
  )
})
