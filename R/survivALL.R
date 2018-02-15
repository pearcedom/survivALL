#' Calculate and combine hazard ratio, pvalue, threshold and
#' area-between-curve data as a single dataframe
#' @inheritParams allPvals
#' @inheritParams bootstrapThresholds
#' @param measure_name A descriptive name for the measure used, for example
#' a gene ID
#' @return a dataframe detailing survival, measure, hazard ratio,
#' pvalue, log10 pvalue, threshold and threshold residual information
#' @examples
#' library(survivALL)
#' data(nki_subset)
#' library(Biobase)
#' library(ggplot2)
#'
#' gene_vec <- exprs(nki_subset)["NM_004448", ] #ERBB2 gene id
#' 
#' survivALL_dfr <- survivALL(measure = gene_vec, 
#'    srv = pData(nki_subset), 
#'    time = "t.dmfs", 
#'    event = "e.dmfs")
#' 
#' ggplot(survivALL_dfr, aes_string(x = 'measure', y = 'p')) + 
#'     geom_hline(yintercept = 0.05, linetype = 3) + 
#'     geom_point()
#' @export
survivALL <- function(measure, 
                      srv, 
                      time = "Time", 
                      event = "Event", 
                      bs_dfr = c(), 
                      measure_name = "measure",
                      multiv = NULL,
                      n_sd = 1.96) {
    # In case calculating with no threshold data
    missing_bs <- is.null(bs_dfr)
    if (missing_bs) {
        bs_dfr <- data.frame(matrix(0, ncol = 3, nrow = nrow(srv)))
        row.names(bs_dfr) <- 1:nrow(bs_dfr)
    } else {
        bs_dfr <- bs_dfr
    }

    # Calculate thresholds from bootstrapping results
    threshold_dfr <- bootstrapThresholds(bs_dfr, n_sd = n_sd)

    # Determine sample order
    measure_ord <- order(measure)
    id_ord <- row.names(srv)[measure_ord]

    # Calculate hazard ratios and p-values
    hratios  <- allHR(measure, srv, time = time, event = event)
    pvals <- allPvals(measure, srv, time = time, event = event, multiv = multiv)
    ## For our logged p-values, we also define non-significant points of 
    ## separation as NA - i.e. they will not plot
    pvals_sig <- pvals
    pvals_sig[pvals_sig >= 0.05] <- NA
    ## Bootstrap p-value calculation
    bsp <- hrSignificance(hratios, bs_dfr)
    bsp_sig <- bsp 
    bsp_sig[bsp_sig >= 0.05] <- NA

    # Combine sample, event, measure, hazard ratio, p-value, desirability and 
    # threshold information as a data frame for plotting
    base_dfr <- data.frame(
                           #samples ordered by measure
                           samples = factor(id_ord, levels = id_ord), 
                           #time to event ordered by measure
                           event_time = srv[[time]][measure_ord], 
                           #events ordered by measure
                           event = ifelse(srv[[event]][measure_ord] == 1, 
                                          TRUE, 
                                          FALSE),
                           #the measure itself, for instance a vector of gene 
                           #expression
                           measure = measure[measure_ord], 
                           #calculated hazard ratios
                           HR = hratios, 
                           #calculated p-values
                           p = pvals, 
                           #corrected p-values
                           p_adj = stats::p.adjust(pvals, method = "fdr"),
                           #logged p-values, 
                           log10_p = log10(pvals_sig), 
                           #bootstrap HR p-values
                           bsp = bsp,
                           #corrected bsp
                           bsp_adj = stats::p.adjust(bsp, method = "fdr"),
                           #a ranking index
                           index = 1:nrow(srv), 
                           #the measure under investigation
                           name = measure_name  
                           )

    dfr <- merge(base_dfr, threshold_dfr, by = 'index', all = TRUE)
    #ranking index is converted to a character class by the above step, 
    #convert back
    dfr$index <- as.numeric(dfr$index) 
    # Calculate area between curves
    ## Hazard ratios are considered relative to the bootstrapped thresholds so 
    ## we calculate by how far a hazard ratio either 
    ## exceeds or does not exceed it's threshold value.
    #calculate for every separation point
    threshold_residuals <- sapply(1:nrow(dfr), function(x) 
                                  areaBetweenCurves(dfr[x,])) 
    ## Combine with dfr
    dfr$threshold_residuals <- threshold_residuals
    ## And add another column with only those values that exceed the thresholds
    dfr$threshold_sig <- threshold_residuals
    dfr$threshold_sig[dfr$threshold_sig <= 0] <- NA

    # Reduce hazard ratios and p-values to a single measure: desirability
    ## and use the most desirable cut point to create a classifier
    ## But if no hazard ratios exceed the thresholds, or no p-values are 
    ## significant, desirability = NA
    if(min(dfr$bsp, na.rm = TRUE) > 0.05 & min(dfr$p, na.rm = TRUE) > 0.05){
        dfr$clsf <- NA
    } else {
        n <- nrow(dfr)
        if(missing_bs){
            dichot_index <- which.min(dfr$p)
            dfr$clsf <- rep(c(0, 1), c(dichot_index, n - dichot_index)) 
        } else {
            #we specify three factors - HR, pvalue and distance from flank - and 
            #combine as a single value, desirability
            d_hr <- desiR::d.high(dfr$threshold_residuals, 
                                  cut1 = 0, 
                                  cut2 = max(dfr$threshold_residuals, na.rm = TRUE))
            d_p <- desiR::d.low(dfr$bsp, 
                                cut1 = min(dfr$bsp, na.rm = TRUE), 
                                cut2 = 0.05)
            #d_middle <- ifelse(dfr$index < n/15 | dfr$index > n - n/15, 0, 1)
            dsr <- desiR::d.overall(d_hr, d_p)#, d_middle)
            dsr[dsr == 0] <- NA
            dsr[dsr == 1] <- NA #likely spurious caused by proximety to edges
            #dfr$most_dsr <- dfr$dsr > stats::quantile(dfr$dsr, na.rm = TRUE)[["75%"]]
            #then, using the most desirable point we produce a dichotomous 
            #classifier
            dichot_index <- which.max(dsr)
            dfr$clsf <- rep(c(0, 1), c(dichot_index, n - dichot_index))
        }
    }

    return(dfr)

}

