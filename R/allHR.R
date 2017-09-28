#' For all possible separation points for a cohort ordered by a continuous
#' measurement, calculate hazard ratio
#' @inheritParams allPvals
#' @param log2HR Hazard ratios are returned as log2 values by default
#' @return A vector of hazard ratios calculated from \code{srv} ordered by 
#' \code{measure}
#' @examples
#' library(survivALL)
#' data(nki_subset)
#' library(Biobase)
#' gene_vec <- exprs(nki_subset)["NM_004448", ] #ERBB2 gene id
#' allHR(measure = gene_vec, srv = pData(nki_subset), time = "t.dmfs", 
#'     event = "e.dmfs", log2HR = TRUE)
#' @export
allHR <- function(measure, srv, time = "Time", event = "Event", log2HR = TRUE) {
    # House keeping 
    if (any(is.na(srv[[time]]) | 
            is.na(srv[[event]]))) stop("NAs in survival data")
    if (any(is.na(measure))) stop("NAs in measure data")
    
    # Order and arrange survival object
    srv_ordered <- srv[order(measure),]
    #srv_dt <- data.table::as.data.table(srv)
    #srv_ordered <- srv_dt[order(measure)]
    srv_time <- srv_ordered[[time]]
    srv_event <- srv_ordered[[event]]
    
    # Create separation list
    separations <- lapply(1:(nrow(srv) - 1), function(x) {
        rep(c(1, 2), c(x, nrow(srv) - x))
    })
    
    # Calculate hazard ratios
    hr_vec_base <- suppressWarnings(sapply(separations, function(x) {
        survcomp::hazard.ratio(x, srv_time, srv_event)$hazard.ratio
    }))
    hr_vec <- her_vec_base#removeOutliers(hr_vec_base)
    #a terminal NA makes the result play well with other variables - e.g. the 
    #number of HRs is n-1 samples, so to align HRs against samples the 
    #additional NA makes this possible
    hr_vec <- c(hr_vec, NA) 
    
    
    # Return logged hazard ratios
    if (isTRUE(log2HR)) {
        log2(hr_vec)
    } else {
        hr_vec
    }
}
