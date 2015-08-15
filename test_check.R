# testing script to make sure that the code works

library(testthat)

source("soap_check.R")


context("Silly checks")

test_that("silly checks",{
  # not a list
  expect_error(soap_check(runif(10)))
  # polygon <4 vertices
  expect_error(soap_check(list(x=c(1,2),y=c(1,2))))
})



context("Ramsay horseshoe")

fsb <- list(fs.boundary())
nmax <- 100
## create some internal knots...
knots <- data.frame(x=rep(seq(-.5,3,by=.5),4),
                    y=rep(c(-.6,-.3,.3,.6),rep(8,4)))
## Simulate some fitting data, inside boundary...
set.seed(0)
n<-600
dat <- data.frame(x = runif(n)*5-1,
                  y = runif(n)*2-1)
x <- dat$x; y <- dat$y
ind <- inSide(fsb,x,y) ## remove outsiders
dat <- dat[ind,]

test_that("boundary works",{
  expect_true(soap_check(fsb, plot=FALSE))
})

test_that("knots work",{
  expect_true(soap_check(fsb, knots, plot=FALSE))
})

test_that("data works",{
  expect_true(soap_check(fsb, data=dat, plot=FALSE))
})

test_that("bad boundary doesn't work...",{
  bnd2 <- fsb
  # start != end
  bnd2[[1]]$x <- bnd2[[1]]$x[-length(bnd2[[1]]$x)]
  bnd2[[1]]$y <- bnd2[[1]]$y[-length(bnd2[[1]]$y)]
  expect_error(soap_check(bnd2))

  # overlapping
  #bnd2 <- c(fsb, fsb)
  #expect_error(soap_check(bnd2))

})

bad_knots <- rbind(knots,c(-1,-1))
test_that("bad knots fail",{
  expect_error(soap_check(fsb, bad_knots, plot=FALSE))
})

context("Ramsay horseshoe (inside out)")

fsb <- fs.boundary()
fsb <- list(fsb,
            list(x=range(fsb$x)[c(1,1,2,2,1)],#+c(-1,-1,1,1,-1),
                 y=range(fsb$y)[c(1,2,2,1,1)]))#+c(-1,1,1,-1,-1)))

test_that("check the boundary works",{
  expect_true(soap_check(fsb, plot=FALSE))
})

test_that("check when the boundary doesn't work...",{
  bnd2 <- fsb
  # start != end
  bnd2[[1]]$x <- bnd2[[1]]$x[-length(bnd2[[1]]$x)]
  bnd2[[1]]$y <- bnd2[[1]]$y[-length(bnd2[[1]]$y)]
  expect_error(soap_check(bnd2))

  # overlapping
  #bnd2 <- c(fsb, fsb)
  #expect_error(soap_check(bnd2))

})

# using the knots from the previous example fails
test_that("bad knots fail",{
  expect_error(soap_check(fsb, bad_knots, plot=FALSE))
})

test_that("bad data fails",{
  expect_error(soap_check(fsb, data=dat, plot=FALSE))
})

dat <- data.frame(x = runif(n)*5-1,
                  y = runif(n)*2-1)
x <- dat$x; y <- dat$y
ind <- inSide(fsb,x,y) ## remove outsiders
dat <- dat[ind,]

test_that("data works",{
  expect_true(soap_check(fsb, data=dat, plot=FALSE))
})

#context("Kodiak islands")

#load("kodiak_example.RData")




