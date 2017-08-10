
library(survivALL)
context("Arrange the major calculations for all cutpoint analysis")

###################
#Make example data#
###################
library(Biobase)
data(nki_subset)
nki_srvcomplete <- nki_subset
gene_vec <- exprs(nki_srvcomplete)["NM_004448", ] #ERBB2 gene id

srvall_dfr <- survivALL(measure = gene_vec, srv = pData(nki_srvcomplete), time = "t.dmfs", event = "e.dmfs", bs_dfr = c(), measure_name = "ERBB2", n_sd = 1.96)

capless_dfr <- srvall_dfr[-nrow(srvall_dfr),]

expct_cls <- c("numeric", "factor", "integer", "logical", "numeric", "numeric", "numeric", "numeric", "factor", "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "numeric")

###################
#Test example data#
###################

test_that("output structure is as expected", {
              expect_is(srvall_dfr, "data.frame")
              identical(as.vector(sapply(srvall_dfr, class)), expct_cls)
})

test_that("pvalues are within range and that non-significant pvalues are NAs in $log10_p", {
              expect_true(all(capless_dfr$p <= 1))
              expect_true(all(capless_dfr$p >= 0))
              expect_identical(capless_dfr$p > 0.05, is.na(capless_dfr$log10_p))
})

test_that("samples have been reordered correctly", {
              expect_true(all(srvall_dfr$measure == sort(gene_vec)))
              row <- sample(1:nrow(srvall_dfr), 1)
              id <- as.character(srvall_dfr$samples[[row]])
              srvall_id <- srvall_dfr[row,]
              pheno_id <- pData(nki_srvcomplete)[which(nki_srvcomplete$samplename == id),]
              identical(as.character(srvall_id$samples), as.character(pheno_id$samplename))
              identical(as.character(srvall_id$event_time), as.character(pheno_id$t.dmfs))
              identical(srvall_id$event, as.logical(pheno_id$e.dmfs))
})


