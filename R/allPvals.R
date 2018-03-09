#' For all possible separation points for a cohort ordered by a continuous 
#' measurement, perform a uni- or multivariate log-rank test 
#' @param measure A continuous variable used to order survival data. Samples
#' must be ordered exactly as in srv
#' @param srv A dataframe that contains at least two columns, detailing event 
#' and time to event information. 
#' Samples must be ordered exactly as in measure
#' @param time Column name in srv containing time to event information. Must 
#' not contain NAs
#' @param event Column name in srv containing event information coded as 0 (no
#' event) and 1 (event). Must not contain NAs
#' @param multiv Univariate analysis is performed by default, however a vector 
#' of additional variables (corresponding to colnames in srv) can be included
#' @param statistic the statistical test to be used to compute significance.
#' one of "logtest" (likelihood ratio test), "waldtest" (wald statistic) or "sctest"
#' (log-rank test)
#' @return A vector of pvalues calculated from \code{srv} ordered by 
#' \code{measure}
#' @examples
#' library(survivALL)
#' data(nki_subset)
#' library(Biobase)
#' gene_vec <- exprs(nki_subset)["NM_004448", ] #ERBB2 gene id
#'
#' allPvals(measure = gene_vec, 
#'     srv = pData(nki_subset), 
#'     time = "t.dmfs", 
#'     event = "e.dmfs",
#'     statistic = "logtest")
#' @export
allPvals <- function(measure, 
                     srv, 
                     time = "Time", 
                     event = "Event", 
                     multiv = NULL, 
                     statistic = "logtest") {
    # House keeping 
    if (any(is.na(srv[[time]]) | 
            is.na(srv[[event]]))) stop("NAs in survival data")
    if (any(is.na(measure))) stop("NAs in measure data")

    # Order and arrange survival object
    srv_ordered <- srv[order(measure),]
    #srv_dt <- data.table::as.data.table(srv)
    #srv_ordered <- srv_dt[order(measure)]
    srv_obj <- survival::Surv(srv_ordered[[time]], srv_ordered[[event]])

    # Create separation list
    separations <- lapply(1:(nrow(srv) - 1), function(x) {
        rep(c(1, 2), c(x, nrow(srv) - x))
    })

    # Create coxph formula for uni- or multivariate analysis
    if(!is.null(multiv)){
        cox_vars <- paste(c("x", multiv), collapse = " + ")
    } else {
        cox_vars <- "x"
    }

    # Calculate P-values 
    pvals <- sapply(separations, function(x) {
                    cox_form <- stats::as.formula(paste("srv_obj ~", cox_vars))
                    summary(survival::coxph(cox_form, srv))[[statistic]][["pvalue"]]
    })
    #terminal NA makes the result play well with other variables - e.g. the 
    #number of pvalues are n-1 samples, so to align HRs against samples the 
    #additional NA makes this possible
    pvals_out <- c(pvals, NA)
    pvals_out
}

