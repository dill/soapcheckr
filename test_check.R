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
#nmax <- 100
### create some internal knots...
#knots <- data.frame(v=rep(seq(-.5,3,by=.5),4),
#                    w=rep(c(-.6,-.3,.3,.6),rep(8,4)))
### Simulate some fitting data, inside boundary...
#set.seed(0)
#n<-600
#v <- runif(n)*5-1;w<-runif(n)*2-1
#y <- fs.test(v,w,b=1)
#names(fsb[[1]]) <- c("v","w")
#ind <- inSide(fsb,x=v,y=w) ## remove outsiders
#y <- y + rnorm(n)*.3 ## add noise
#y <- y[ind];v <- v[ind]; w <- w[ind]
#n <- length(y)

test_that("check it works",{
  expect_true(soap_check(fsb, plot=FALSE))
})

test_that("check it doesn't work...",{
  bnd2 <- fsb
  # start != end
  bnd2[[1]]$x <- bnd2[[1]]$x[-length(bnd2[[1]]$x)]
  bnd2[[1]]$y <- bnd2[[1]]$y[-length(bnd2[[1]]$y)]
  expect_error(soap_check(bnd2))

  # overlapping
  #bnd2 <- c(fsb, fsb)
  #expect_error(soap_check(bnd2))

})


context("Ramsay horseshoe -- inside out")

fsb <- fs.boundary()
fsb <- list(fsb,
            list(x=range(fsb$x)[c(1,1,2,2,1)],
                 y=range(fsb$y)[c(1,2,2,1,1)]))

test_that("check it works",{
  expect_true(soap_check(fsb, plot=FALSE))
})



#context("Kodiak islands")

#load("kodiak_example.RData")




