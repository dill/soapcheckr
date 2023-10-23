# testing script to make sure that the code works

# library(testthat)
#
# source("soap_check.R")


context("Silly checks")

test_that("silly checks",{
  # not a list
  expect_error(soap_check(runif(10)))
  # polygon <4 vertices
  expect_error(soap_check(list(x = c(1, 2),
                               y = c(1, 2))))
})



context("Ramsay horseshoe")

fsb <- list(fs.boundary())
nmax <- 100
## create some internal knots...
knots <- data.frame(x = rep(seq(-0.5, 3, by = 0.5), 4),
                    y = rep(c(-0.6, -0.3, 0.3, 0.6), rep(8, 4)))
## Simulate some fitting data, inside boundary...
set.seed(0)
n <- 600
dat <- data.frame(x = runif(n) * 5 - 1,
                  y = runif(n) * 2 - 1)
x <- dat$x; y <- dat$y
ind <- inSide(fsb, x, y) ## remove outsiders
dat <- dat[ind,]

test_that("boundary works",{
  expect_true(
    soap_check(fsb, plot = FALSE)
  )
})

test_that("knots work",{
  expect_true(
    soap_check(fsb, knots = knots, plot = FALSE)
  )
})

test_that("data works",{
  expect_true(
    soap_check(fsb, data = dat, plot = FALSE)
  )
})

test_that("bad boundary doesn't work...",{
  bnd2 <- fsb
  # start != end
  bnd2[[1]]$x <- bnd2[[1]]$x[-length(bnd2[[1]]$x)]
  bnd2[[1]]$y <- bnd2[[1]]$y[-length(bnd2[[1]]$y)]
  expect_error(
    soap_check(bnd2)
  )

  # overlapping
  #bnd2 <- c(fsb, fsb)
  #expect_error(soap_check(bnd2))

})

# Test bad knots

bad_knots <- rbind(knots, c(-1, -0.5))

test_that("bad knots fail",{

  expect_warning(
    soap_check(fsb, knots = bad_knots, plot = FALSE)
  )
})

context("Ramsay horseshoe (inside out)")

fsb <- fs.boundary()
fsb <- list(fsb,
            list(x = range(fsb$x)[c(1, 1, 2, 2, 1)],#+c(-1,-1,1,1,-1),
                 y = range(fsb$y)[c(1, 2, 2, 1, 1)]))#+c(-1,1,1,-1,-1)))


test_that("check the boundary works",{
  expect_true(
    soap_check(fsb, plot = FALSE)
  )
})

test_that("check when the boundary doesn't work...",{
  bnd2 <- fsb
  # start != end
  bnd2[[1]]$x <- bnd2[[1]]$x[-length(bnd2[[1]]$x)]
  bnd2[[1]]$y <- bnd2[[1]]$y[-length(bnd2[[1]]$y)]
  expect_error(
    soap_check(bnd2)
  )

  # i'm not sure this is testting the rigt thing
    # overlapping
  bnd2 <- c(fsb, fsb)

  # expect_false(
      # soap_check(bnd2)
  #   )

})

test_that("bad knots fail",{
  expect_warning(
    soap_check(fsb, bad_knots, plot = FALSE)
  )
})

test_that("bad data fails",{
  expect_warning(
    soap_check(fsb, data = dat, plot = FALSE)
    )
})

dat <- data.frame(x = runif(n) * 5 - 1,
                  y = runif(n) * 2 - 1)
x <- dat$x; y <- dat$y
ind <- inSide(fsb, x, y) ## remove outsiders
dat <- dat[ind,]

test_that("data works",{
  expect_true(
    soap_check(fsb, data = dat, plot = FALSE)
  )
})

#context("Kodiak islands")

#load("kodiak_example.RData")

context("Overlapping Polygons")

# library(testthat)
# expect_fa

test_that("check when the boundary overlaps",{

  # overlapping

  expect_false(
    expect_warning(
      soap_check(overlap_polygons, plot = FALSE)
    )
  )

})
context("Test Sissabagama Lake Polygons")
test_that("Check when it has multiple polygons",{
  expect_true(
    soap_check(sissabagama_lake_ls)
  )

})
test_that("Check if ID knots to close to boundrary ",{
  expect_true(
    soap_check(sissabagama_lake_ls,
               knots = sissabagama_lake_knots,
               plot = FALSE)
  )
})
test_that("Check if ID knots to outside to boundrary ",{
  expect_true(expect_warning(
    soap_check(sissabagama_lake_ls,
               knots = sissabagama_lake_grid_knots,
               plot = FALSE
               )
  )
  )

})
