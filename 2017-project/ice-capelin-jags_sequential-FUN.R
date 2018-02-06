#' plotCredInt
#' plot the credible and prediction intervals
#' @param df 
#' @param yaxis - variable in df for yaxis
#' @param ylab - customize the label of the yaxis
#' @param y_line - predicted values, e.g. y_med
#' @param ci - credible interval
#' @param pi - prediction interval
#' @param model - the model being tested - use a character object
#' @param x location of model text
#' @param y location of model text
#'
#' @return
#' @export
#'
#' @examples plotCredInt(df2, yaxis = "ln_biomass_med", ylab = "ln(capelin)", y_line = y_med, ci_df2, pi_df2, model = txt, x = 2006, y = 8)

plotCredInt <- function(df, yaxis = yaxis, ylab = ylab, y_line = y_line, ci = ci, dpi = dpi, model = model, x = x, y = y){
     
     p <- ggplot()
     #browser()
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5)
     p <- p + xlab("Year") + ylab(paste(ylab))
     p <- p + theme(text = element_text(size=15)) + theme_bw()
     p <- p + geom_line(aes(x = c(df$year, 2018:2019), y = y_line))
     p <- p + geom_ribbon(aes(x = c(df$year, 2018:2019), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5)
     pi_n <- dpi[, (ncol(dpi)-1):ncol(dpi)]
     p <- p + geom_ribbon(aes(x = c(2018:2019), 
                              ymax = pi_n[2, ], 
                              ymin = pi_n[1, ]),
                          alpha = 0.5)
     p <- p + geom_text() +
          annotate("text", label = model, x = x, y = y)
     
     return(p)
     
}