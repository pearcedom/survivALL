#' Calculate the area above the bootstrapped thresholds but below the HR 
#' distribution and vice versa for each point of separation
#' @param survivALL_dfr Output of survivALL()
#' @return A per-separation point vector equal to the distance beyond or within
#' the bootstrapped thresholds
areaBetweenCurves = function(survivALL_dfr){
    if (is.na(survivALL_dfr$HR)) {
        NA
    } else {
        if (survivALL_dfr$HR >= 0) { #if HR is positive...
         return(survivALL_dfr$HR - survivALL_dfr$sdplus) #compare to upper threshold
        } else { #but if negative...
            return((survivALL_dfr$HR - survivALL_dfr$sdmin) * -1) #compare to lower threshold
        }
    }
}
