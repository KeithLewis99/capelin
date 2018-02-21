#' plotCredInt---
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
     p <- p + geom_ribbon(aes(x = c(df$year, 2018:2019), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5)
     pi_n <- dpi[, (ncol(dpi)-1):ncol(dpi)]
     p <- p + geom_ribbon(aes(x = c(2018:2019), 
                              ymax = pi_n[2, ], 
                              ymin = pi_n[1, ]))
     p <- p + geom_text() +
          annotate("text", label = model, x = x, y = y, size = 7.5)
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5)
     p <- p + xlab("Year") + ylab(paste(ylab))
     p <- p + theme(text = element_text(size=15)) + theme_bw()
     p <- p + geom_line(aes(x = c(df$year, 2018:2019), y = y_line))
     return(p)
}


#' prioPosterior----
#' Plots the values of the posterior with the prior overlaid on top
#' @param df1 - vector of values for posterior
#' @param df2 - vector of values for prior
#' @param xlim - upper and lowe values for df1 - approximate
#'
#' @return
#' @export
#'
#' @examples priorPosterior(beta_post, dist, x)
priorPosterior <- function(df1, df2, xlim){
     p <- ggplot()
     p <- p + geom_histogram(data = df1, aes(x=V1, y = ..ncount..), colour="black", fill="white") + theme(axis.title = element_blank())
     p <- p + geom_density(data = df2, aes(v2), colour = "red") + xlim(xlim)  + theme(axis.title = element_blank())
     p
}





#' Title
#'
#' @param df 
#' @param df2 
#' @param df3 
#' @param limits 
#' @param x_label 
#' @param priormean 
#' @param priorsd 
#' @param by_bin 
#'
#' @return
#' @export
#'
#' @examples
postPriors <- function(df, df2, df3, limits = limits, x_label = x_label, priormean, priorsd, by_bin = 1){
     #browser()
     p <- ggplot() + theme_bw()
     p <- p + coord_cartesian(xlim = c(limits[1], limits[2])) 
     p <- p + geom_histogram(data = as.data.frame(df), 
                             aes(x = V1, y=..density..), 
                             col="grey", fill="grey",
                             #alpha=.6,
                             breaks = seq(limits[1], limits[2], by=by_bin)) + xlab(x_label) + ylab("Density")
     p <- p + geom_density(data = as.data.frame(df2), aes(x=df2), size = 1, alpha = 0.3) 
     p <- p + geom_rug(data=data.frame(y=df3), aes(x=y), colour = "red", size=3, alpha=1) 
     p <- p + geom_vline(aes(xintercept = 0), colour = "red", size = 2)
     return(p)
}




# not using - not appropriate
DIC_out <- function(df){
     browser()
     dic_mod <-dic.samples(df$model, 1000, "pD")
     dic_val <- sum(dic_mod$deviance) + sum(dic_mod$penalty)
     return(dic_val)
     
}


