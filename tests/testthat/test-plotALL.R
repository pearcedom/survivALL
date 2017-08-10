

library(survivALL)
context("Check that plotALL produces a ggplot object")

###################
#Make example data#
###################
library(Biobase)
library(ggplot2)
data(nki_subset)
nki_srvcomplete <- nki_subset

gene_vec <- exprs(nki_srvcomplete)["NM_004448", ] #ERBB2 gene id

srvall_p <- plotALL(measure = gene_vec, srv = pData(nki_srvcomplete), time = "t.dmfs", event = "e.dmfs", bs_dfr = c(), measure_name = "ERBB2")

###################
#Test example data#
###################

test_that("plotALL output is a ggplot object", {
              expect_true(is.ggplot(srvall_p))
})
