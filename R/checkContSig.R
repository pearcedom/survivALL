#' Calculate association between survival and a continuous measure#' @inheritParams allPvals
#' @param measure A continuous variable used to order survival data. Samples
#' must be ordered exactly as in time and event
#' @param time A numeric vector of sample time-to-event data, ordered exactly as 
#' measure. Must not contain NAs
#' @param event A vector of sample event data, ordered exactly as measure. Must not 
#' contain NAs
#' @return p-value of association between measure and survival
#' @examples
#' library(survivALL)
#' library(Biobase)
#' data(nki_subset)
#' 
#' #Calculate p-value for continuous measure SCUBE2
#' srv_dfr <- data.frame(measure = exprs(nki_subset)["NM_020974", ],
#'                       time = nki_subset$t.dmfs,
#'                       event = nki_subset$e.dmfs
#'                       )
#' 
#' checkContSig(srv_dfr$measure, srv_dfr$time, srv_dfr$event)
#' @export
checkContSig <- function(measure, time, event){
    srv_obj <- survival::Surv(time, event)
    broom::tidy(survival::coxph(srv_obj ~ measure))$p.value
}
