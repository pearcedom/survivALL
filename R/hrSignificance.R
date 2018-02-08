#' Calculate HR significance using bootstrap results
#' @param survivALL_dfr Output of survivALL()
#' @return A per-separation point vector equal to the distance beyond or within
#' the bootstrapped thresholds
hrSignificance <- function(hratios, bs_dfr){
    if(all(bs_dfr == 0)){
        bsp <- NA
    } else {
        n <- ncol(bs_dfr)
        means <- apply(bs_dfr, 1, function(x) mean(x, na.rm = TRUE))
        obs_min_exp <-  abs(hratios - means)
        excd_up <- means + obs_min_exp
        excd_down <- means - obs_min_exp

        up_sums <- rowSums(bs_dfr > excd_up, na.rm = TRUE)
        down_sums <- rowSums(bs_dfr < excd_down, na.rm = TRUE)
        bsp <- (up_sums + down_sums) / n
        bsp[is.na(hratios)] <- NA
    }
    bsp
}
