#' Calculate and combine hazard ratio, pvalue, threshold and
#' area-between-curve data and plot
#' @inheritParams allPvals
#' @inheritParams bootstrapThresholds
#' @inheritParams survivALL
#' @param plot_range Allows manual specification of the y-axis; 
#' as c(lower, upper)
#' @param scale_upper Allows manual specification of the pvalue colour bar upper
#' limit; e.g. 3
#' @param title Plot title; as a character string
#' @param point_size Hazard ratio point size; e.g. 1.1
#' @param legend Legend position; one of "top", "bottom", "left", "right"
#' @param timeplot_type Determine how event information is displayed; one of 
#' "tiles" or "lollies"
#' @param axis_text Removing axis text allows plots to be more easily compared;
#' logical
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
#' @export
plotALL <- function(measure, 
                    srv, 
                    time = "Time", 
                    event = "Event", 
                    bs_dfr = c(),
                    measure_name = "measure",
                    multiv = NULL,
                    plot_range = "auto", 
                    scale_upper = "auto", 
                    title = "", 
                    point_size = 1.5, 
                    legend = 'bottom', 
                    timeplot_type = "tiles", 
                    axis_text = TRUE) {

    dfr <- survivALL(measure = measure, 
                     srv = srv, 
                     measure_name = measure_name, 
                     time = time, 
                     event = event, 
                     bs_dfr = bs_dfr, 
                     multiv = multiv)

    # The y-axis range and colour scales can be manually specified to facilitate
    # comparisons between plots
    ## Define the plot's y-range
    ### If not specified this will simply be the HR range
    if (plot_range == "auto"){
        plot_range <- c(min(min(dfr$HR, na.rm = TRUE), -2),
                        max(max(dfr$HR, na.rm = TRUE), 2))
    } else {
        plot_range <- plot_range
    }

    ### If tiles are used to display events then we increase the y-axis by 20% 
    ### to accommodate this
    if (timeplot_type == "tiles"){
        plot_range[1] <- plot_range[1] - diff(plot_range)/20
    }
    
    ## Similarly, we fine-tune the p-value colour range
    ### The lower end of this scale will always by 1.30103 (= log10(0.05))
    if (scale_upper == "auto"){
        colour_scale <- ggplot2::scale_colour_gradientn(
                            colours = viridis::viridis(option = "D", 
                                                       begin = 0, 
                                                       end = 1, 
                                                       n = 75),
            limits = range(dfr$log10_p, na.rm = TRUE))
    } else {
        colour_scale <- ggplot2::scale_colour_gradientn(
                            colours = viridis::viridis(option = "D", 
                                                       begin = 0, 
                                                       end = 1, 
                                                       n = 75),
            limits = c(1.30103, scale_upper))
    }

    # Plot
    ## Define plot dimensions, labels, horizontal zero-line and thresholds 
    p1 <- ggplot2::ggplot(dfr, ggplot2::aes_string(x = 'index', y = 'HR')) +
        ggplot2::geom_hline(yintercept = 0, linetype = 9) + #zero-line
        ggplot2::geom_line(ggplot2::aes_string(y = 'sdplus', x = 'index'), 
                             linetype = 2, 
                             pch = "-",
                             colour = "DARK GREY") + #upper threshold
        ggplot2::geom_line(ggplot2::aes_string(y = 'sdmin', x = 'index'), 
                             linetype = 2, 
                             colour = "DARK GREY") + #lower threshold
        ggplot2::scale_fill_manual(values = c("BLACK"), guide = FALSE) + 
            ggplot2::scale_alpha(range = c(1, 0.1), guide = FALSE) + 
        ggthemes::theme_pander() +

        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
              legend.position = legend, #determines legend position
              legend.key.height = ggplot2::unit(0.4, "cm"),
              axis.ticks.x = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5)
              ) +
        colour_scale + #as defined above
        ggplot2::ylim(plot_range) + #as defined above
        ggplot2::ggtitle(title)
    
    ## Alternative axis-less theme. Plotting with no axis text can help when 
    ## comparing plots
    if(axis_text == FALSE){ 
        p1 <- p1 + ggplot2::theme(axis.text = ggplot2::element_blank(), 
                                  axis.ticks = ggplot2::element_blank(), 
                                  axis.title.y = ggplot2::element_blank())
    } else {
        p1 <- p1
    }

    ## Add hazard ratio layer, coloured by significance
   # p2 <- if(all(is.na(dfr$log10_p))){ #remind myself what this check is for
   #   p1 + ggplot2::geom_point(size = point_size, colour = "#737373")
   # } else {
   #   p1 + ggplot2::geom_point(size = point_size, 
   #        ggplot2::aes_string(colour = 'log10_p'))
   # }
    p2 <- p1 + ggplot2::geom_point(size = point_size, 
               ggplot2::aes_string(colour = 'threshold_sig'))

    ## Lolly-plots
    ### Time to event information can be display as barplots, with events 
    ### rising up from y=0 and non-events heading down from y=0
    ### These are plotted and then combined with the main plot however, 
    ### meaning that further modifications cannot be made after the function 
    ### is complete

    ### We first polarise our time to event info, so that non-event times 
    ### become negative (event times remain positive)
    dfr$lolly_time <- sapply(1:nrow(dfr), function(x)
      ifelse(dfr$event[x] == TRUE, dfr$event_time[x], (dfr$event_time[x]*-1)))
    ### Set the range and breaks for this plot
    y_range <- round(range(dfr$lolly_time, na.rm = TRUE))
    y_breaks <- round(seq(y_range[1], y_range[2], diff(y_range)/5))
    
    ## Either construct the lolly plot and combine with p2...
    if(timeplot_type == "lollies") {
        p3 <- ggplot2::ggplot(dfr, ggplot2::aes_string(x = 'index', 
                                                       y = 'lolly_time', 
                                                       colour = 'event')) +
            ggplot2::geom_hline(yintercept = 0, linetype = 1) + 
            ggplot2::geom_bar(stat = 'identity', 
                              width = 0.01, 
                              colour = "BLACK") +
            ggplot2::geom_point(size = 0.7) +
            ggthemes::theme_pander() +
            ggplot2::theme(legend.position = 'none',
                    axis.text.x = ggplot2::element_blank(),
                    axis.title.x = ggplot2::element_blank(),
                    axis.ticks.x = ggplot2::element_blank()) +

            ggplot2::scale_y_continuous(breaks = y_breaks, 
                                        labels = abs(y_breaks)) +
            ggplot2::scale_colour_manual(values = c("#85D4E3", "#F4B5BD")) +
            ggplot2::ylab("Time")
        if(axis_text == FALSE){ #axis text free version
            p3 <- p3 + ggplot2::theme(axis.text = ggplot2::element_blank(), 
                                      axis.ticks = ggplot2::element_blank(), 
                                      axis.title.y = ggplot2::element_blank())
        } else {
            p3 <- p3
        }
        cowplot::plot_grid(p2, p3, ncol = 1, rel_heights = c(4, 1))

    ## or the equivalent plot using tiles to display time-to-event information
    } else if(timeplot_type == "tiles") {
        p3 <- p2 + ggplot2::geom_tile(data = dfr,
                ggplot2::aes_string(x = 'index', 
                                    y = plot_range[1] + 
                                        diff(plot_range)/80, fill = 'event'), 
                                      colour = "WHITE", 
                                      height = diff(plot_range)/40) +
                ggplot2::scale_fill_manual(values = c("WHITE", "BLACK")) + 
                #events are therefore black tiles, non-events are negative space
                ggplot2::guides(fill = FALSE)
    } else {
        cat("please specify how you would like events to be plotted!\n")
    }
    p3
  }

