#' Calculate outliers in a numeric vector and then convert those values to NA
#' @param x A numeric vector
#' @param tolerant Outlier calculation tolerance. A more tolerant outlier 
#' removal is more appropriate when working with hazard ratios
#' @return The modified, outlier removed, equivalent of x
#' @examples
#' set.seed(123); x <- rnorm(100)
#' sum(is.na(x))
#' y <- removeOutliers(x)
#' sum(is.na(y))
#' @export
removeOutliers <- function(x, tolerant = TRUE){
    outlierBoundries <- function(x){
        iqr <- stats::IQR(x, na.rm = TRUE)
        qlow <- stats::quantile(x, na.rm = TRUE)[["25%"]]
        qhigh <- stats::quantile(x, na.rm = TRUE)[["75%"]]
        coef <- ifelse(tolerant, 2.5, 1.5)
        lower <- qlow - coef * iqr
        upper <- qhigh + coef * iqr
        c(lower, upper)
    }

    outlier2NA <- function(x){
        bounds <- outlierBoundries(x)
        x[x < bounds[[1]]] <- NA
        x[x > bounds[[2]]] <- NA
        x
    }
    outlier2NA(x)
}

