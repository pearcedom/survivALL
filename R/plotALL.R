#' Calculate and combine hazard ratio, pvalue, threshold and
#' area-between-curve data and plot
#' @inheritParams allPvals
#' @inheritParams bootstrapThresholds
#' @inheritParams survivALL
#' @param title Plot title; as a character string
#' @return Using survival, measure, hazard ratio, pvalue, log10 pvalue, 
#' threshold and threshold residual information, plot the measure-event 
#' relationship
#' @examples
#' data(nki_subset)
#' library(Biobase)
#' 
#' gene_vec <- exprs(nki_subset)["NM_004448", ] #ERBB2 gene id
#' 
#' plotALL(measure = gene_vec, 
#'     srv = pData(nki_subset), 
#'     time = "t.dmfs", 
#'     event = "e.dmfs", 
#'     title = "ERBB2 Example") 
#' 
#' @export
plotALL <- function(measure, 
                    srv, 
                    time = "Time", 
                    event = "Event", 
                    bs_dfr = c(),
                    measure_name = "measure",
                    multiv = NULL,
                    title = "") {
    
    #####################
    ## DATA PROCESSING ##
    #####################
    # In case calculating with no bootstrapping data
    missing_bs <- is.null(bs_dfr)
    if (missing_bs) {
        bs_dfr <- data.frame(matrix(0, ncol = 3, nrow = nrow(srv)))
        row.names(bs_dfr) <- 1:nrow(bs_dfr)
    } else {
        bs_dfr <- bs_dfr
    }

    # Calculate survival statistics with survivALL()
    base_dfr <- survivALL(measure = measure, 
                     srv = srv, 
                     measure_name = measure_name, 
                     time = time, 
                     event = event, 
                     bs_dfr = bs_dfr, 
                     multiv = multiv)
    
    # Calculate thresholds from bootstrapping results
    threshold_dfr <- bootstrapThresholds(bs_dfr)
    # and merge with survivALL() output
    dfr <- merge(base_dfr, threshold_dfr, by = 'index', all = TRUE)
    # ranking index is converted to a character class by the above step, 
    # convert back
    dfr$index <- as.numeric(dfr$index) 

    # To accommodate our event information, we increase the y-axis by 1/20th, from the lower end
    plot_range <- range(dfr$HR, na.rm = TRUE)
    plot_range[1] <- plot_range[1] - diff(plot_range)/20


    ##############
    ## Plotting ##
    ##############
    # Define plot dimensions, labels, horizontal zero-line and thresholds 
    p1 <- ggplot2::ggplot(dfr, ggplot2::aes_string(x = 'index', y = 'HR')) + 
            ggplot2::geom_hline(yintercept = 0, linetype = 9) + #zero-line
            ggplot2::geom_line(ggplot2::aes_string(y = 'sdplus', x = 'index'), 
                               linetype = 2, 
                               colour = "DARK GREY") + #upper threshold
            ggplot2::geom_line(ggplot2::aes_string(y = 'sdmin', x = 'index'), 
                                 linetype = 2, 
                                 colour = "DARK GREY") + #lower threshold
            ggplot2::scale_fill_manual(values = c("BLACK"), guide = FALSE) + 
                ggplot2::scale_alpha(range = c(1, 0.1), guide = FALSE) + 
            ggthemes::theme_pander() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                  legend.position = 'bottom', #determines legend position
                  legend.key.height = ggplot2::unit(0.4, "cm"),
                  axis.ticks.x = ggplot2::element_blank(),
                  axis.title.x = ggplot2::element_blank(),
                  plot.title = ggplot2::element_text(hjust = 0.5)
                  ) +
            ggplot2::ylim(plot_range) + #as defined above
            ggplot2::ggtitle(title)

    # Add hazard ratio layer, coloured by significance
    if(missing_bs){
        ## For non-bootstrap pvalues we fine-tune the colour range
        ### The upper limit of this scale will always be 0.05 (significance)
        colour_scale <- ggplot2::scale_colour_gradientn(
                                         colours = viridis::viridis(option = "D", 
                                                                    begin = 1, 
                                                                    end = 0, 
                                                                    n = 75),
                                                                    limits = range(dfr$log10_p, 
                                                                                   na.rm = TRUE))
        p2 <- p1 + 
            ggplot2::geom_point(ggplot2::aes_string(colour = 'log10_p')) +
            colour_scale
    } else {
        p2 <- p1 + 
            ggplot2::geom_point(ggplot2::aes_string(colour = 'dsr')) +
            ggplot2::scale_colour_gradientn(colours = viridis::viridis(
                                                                       option = "D", 
                                                                       begin = 0, 
                                                                       end = 1, 
                                                                       n = 75)) +
            labs(colour = "Desirability")
    }

    ## Finally we add a visual representation of our events
    p_out <- suppressMessages(
              p2 + 
        ggplot2::geom_tile(data = dfr,
                           ggplot2::aes_string(x = 'index', 
                                               y = plot_range[1] + 
                                                   diff(plot_range)/80, 
                                               fill = 'event'), 
                           colour = "WHITE", 
                           height = diff(plot_range)/40) +
        ggplot2::scale_fill_manual(values = c("WHITE", "BLACK")) + 
        ggplot2::guides(fill = FALSE)
    )
    
    p_out
} 
