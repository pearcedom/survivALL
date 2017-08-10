#' Calculate per-separation point hazard ratio thresholds
#' @param bs_dfr A matrix of bootstrapped hazard ratio computations as ordered
#' by a random measurement vector. Typically consisting of 5-10,000 repeat 
#' samplings
#' @param n_sd The number of standard deviations used to define threshold width.
#' 95% of random hazard ratio values fall within the thresholds with a standard 
#' deviation of 1.96
#' @return A dataframe of per-separation point mean, upper and lower thresholds
#' @examples
#' data(nki_subset)
#' library(Biobase)
#' library(magrittr)
#' library(ggplot2)
#'
#' #simulate example HR bootstrapped data
#' bs_dfr <- matrix(rnorm(150000), ncol = 1000, nrow = 150)
#'
#' #calculate thresholds
#' thresholds <- bootstrapThresholds(bs_dfr)
#' @export
bootstrapThresholds <- function(bs_dfr, n_sd = 1.96) {
    #Calculate per cutpoint SD
    bs_sd <- apply(bs_dfr, 1, function(x) stats::sd(x, na.rm = TRUE))

    #Calculate per cutpoint mean
    bs_mean <- apply(bs_dfr, 1, function(x) mean(x, na.rm = TRUE))

    #Combine and calculate thresholds
    bstrap.dfr <- data.frame(
                             mean = bs_mean,
                             sd = bs_sd,
                             index = 1:nrow(bs_dfr)
                             )

    bstrap.dfr$sdplus <- bstrap.dfr$mean + (n_sd * bs_sd)
    bstrap.dfr$sdmin <- bstrap.dfr$mean - (n_sd * bs_sd)
    bstrap.dfr$sd <- NULL
    bstrap.dfr
}
