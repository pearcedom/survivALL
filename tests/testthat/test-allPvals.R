library(survivALL)
context("P-value calculations")

###################
#Make example data#
###################
library(survsim)
#Default survival data simulation taken from survsim package reference manual
set.seed(123); sim_dfr <- simple.surv.sim(n=200, foltime=3600, 
                                              dist.ev=c("llogistic" ),
                                              anc.ev=c(0.69978200185280),
                                              beta0.ev=c(5.84298525742252),
                                              anc.cens=1.17783687569519,
                                              beta0.cens=7.39773677281100,
                                              z=list(c("unif", 0.8, 1.2)), 
                                              beta=list(c(-0.4), c(0)), 
                                              x=list(c("bern", 0.5), c("unif", 0.7, 1.3)))
srv <- sim_dfr[c(1, 2, 4)]
measure <- sample(nrow(srv))

#Calculate p-values
test_p <- allPvals(measure, srv, event = "status", time = "stop")

#Because a terminal NA is added so that the length of allPvals() is equall to the length of measure, for some tests we remove this
capless_p <- test_p[-length(test_p)]

###################
#Test example data#
###################
test_that("Ps are numeric", {
          expect_true(is.numeric(test_p))
})

test_that("Ps are between 0 and 1", {
          expect_true(all(capless_p >= 0))
          expect_true(all(capless_p <= 1))
})

test_that("There is 1 P per sample", {
          expect_identical(length(test_p), length(measure))
})

test_that("log-rank is identical between survivALL and survdiff calculations", {

})
