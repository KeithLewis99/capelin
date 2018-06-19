#' plotCredInt---
#' plot the credible and prediction intervals
#' @param df - the data frame
#' @param yaxis - variable in df for yaxis
#' @param ylab - customize the label of the yaxis
#' @param y_line - predicted values, e.g. y_med
#' @param ci - credible interval
#' @param pi - prediction interval
#' @param model - the model being tested - use a character object
#' @param x location of model text
#' @param y location of model text
#'
#' @return - a graph of teh capelin biomass with 95% CIs, the fitted values of the model, and the 95% credible and prediction intervals.
#' @export
#'
#' @examples plotCredInt(df2, yaxis = "ln_biomass_med", ylab = "ln(capelin)", y_line = y_med, ci_df2, pi_df2, model = txt, x = 2006, y = 8)

plotCredInt <- function(df, yaxis = yaxis, ylab = ylab, y_line = y_line, ci = ci, dpi = dpi, model = model, x = x, y = y, type = type){
     p <- ggplot()  
     #browser()
     p <- p + geom_ribbon(aes(x = c(df$year, 2018:2019), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5, fill = "grey60")
     pi_n <- dpi[, (ncol(dpi)-1):ncol(dpi)]
     p <- p + geom_ribbon(aes(x = c(2018:2019), 
                              ymax = pi_n[2, ], 
                              ymin = pi_n[1, ]), fill = "grey40")
     #p <- p + geom_text() +
         # annotate("text", label = model, x = x, y = y, size = 7.5)
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5)
    if(!is.na(type)){
         p <- p + geom_errorbar(data = df, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
    } else if (is.na(type)) {
         p
    }
         
     p <- p + xlab("Year") + ylab(paste(ylab))
    # p <- p + theme(axis.title = element_text(size=30), axis.text = element_text(size = 20)) 
     #p <- p + theme(text = element_text(size=25)) + theme_bw()
     p <- p + geom_line(aes(x = c(df$year, 2018:2019), y = y_line))
    
     p <- p + theme_bw(base_size = 30) + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
     return(p)
}

# as above but just for 2018
plotCredInt_a <- function(df, yaxis = yaxis, ylab = ylab, y_line = y_line, ci = ci, dpp = dpp, dpi = dpi, model = model, x = x, y = y, type = type){
     p <- ggplot()  
     #browser()
     p <- p + geom_ribbon(aes(x = c(df$year, 2018), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5, fill = "grey60")
     pi_n <- dpi[, ncol(dpi)]
    # p <- p + geom_ribbon(aes(x = c(2018), 
     #                         ymax = pi_n[2], 
      #                        ymin = pi_n[1]), fill = "grey40")
     
     p <- p + geom_errorbar(aes(x = c(2018), 
                                ymax = pi_n[2], 
                                ymin = pi_n[1]))
     #p <- p + geom_point(aes(x = c(2018), 
      #                       y = last(dpp)),
       #                  colour = "black",
        #                 shape = 17, size = 2)
     #p <- p + geom_text() +
     # annotate("text", label = model, x = x, y = y, size = 7.5)
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5)
     if(!is.na(type)){
          p <- p + geom_errorbar(data = df, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     } else if (is.na(type)) {
          p
     }
     
     p <- p + xlab("Year") + ylab(paste(ylab))
     # p <- p + theme(axis.title = element_text(size=30), axis.text = element_text(size = 20)) 
     #p <- p + theme(text = element_text(size=25)) + theme_bw()
     p <- p + geom_line(aes(x = c(df$year, 2018), y = y_line))
     
     p <- p + theme_bw(base_size = 30) + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
     return(p)
}

#' priorPosterior----DEPRECATED
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


#' postPriors()----
#'
#' @param df - the chain of values of the parameter
#' @param df2 - the distribution of the prior based on randomly generated data
#' @param df3 - the chain values within the credible intervals of the parameter
#' @param limits - an ojbect to set the limits on the grid 
#' @param x_label - label the graph with the name of the parm
#' @param priormean -  the mean value of the prior
#' @param priorsd - the SD of the prior
#' @param by_bin - 
#'
#' @return - figure of the posterior distribution, the prior distribution and the 95% credible interval as a rug.  Abline = 0
#' @export
#'
#' @example p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)
postPriors <- function(df, df2, df3, limits = limits, x_label = x_label, priormean, priorsd, by_bin = 1, vline="yes"){
     #browser()
     p <- ggplot() + theme_bw(base_size = 20)
     p <- p + coord_cartesian(xlim = c(limits[1], limits[2])) 
     p <- p + geom_histogram(data = as.data.frame(df), 
                             aes(x = V1, y=..density..), 
                             col="grey", fill="grey",
                             #alpha=.6,
                             breaks = seq(limits[1], limits[2], by=by_bin)) + xlab(x_label) + ylab("Density")
     p <- p + geom_density(data = as.data.frame(df2), aes(x=df2), size = 1, alpha = 0.3) 
     p <- p + geom_rug(data=data.frame(y=df3), aes(x=y), colour = "red", size=3, alpha=1) 
     if(vline == "yes"){
          p <- p + geom_vline(aes(xintercept = 0), colour = "red", size = 2)
          }
     return(p)
}


#' Posterior_fig----
#' helper function for extracting credible interval for a given parameter
#' @param ls - the output of the JAGS model - extact the chain for a given parameter
#'
#' @return - list with the chain values for the parameter, the 95% CI, and the chain values within credible intervals
#' @export
#'
#' @examples - alpha <- posterior_fig(out$sims.list$alpha)
posterior_fig <- function(ls){
     #browser()
     df <- ls
     df_quant <- quantile(ls, c(0.025, 0.975))
     df_cred <- subset(df, df > df_quant[1] & df < df_quant[2])
     return(list(df = df, df_quant = df_quant, df_cred = df_cred))
}

#' Posterior_fig1----
#' helper function for extracting credible interval for a given parameter
#' @param ls - the output of the JAGS model - extact the chain for a given parameter
#'
#' @return - list with the chain values for the parameter, the 95% CI, and the chain values within credible intervals
#' @export
#'
#' @examples - gamma <- posterior_fig1(out$sims.list$gamma, transform = "yes", parm = "gamma")
posterior_fig1 <- function(ls, transform = transform, parm = parm){
     #browser()
     if(transform == "no"){
          df <- ls
          df_quant <- quantile(ls, c(0.025, 0.975))
          df_cred <- subset(df, df > df_quant[1] & df < df_quant[2])
     } else if (transform == "yes" & parm == "slope") {
               df <- ls/100
               df_quant <- quantile(df, c(0.025, 0.975))
               df_cred <- subset(df, df > df_quant[1] & df < df_quant[2])
     } else if(transform == "yes" & parm == "width") {
          df <- ls*100
          df_quant <- quantile(df, c(0.025, 0.975))
          df_cred <- subset(df, df > df_quant[1] & df < df_quant[2])
     }
     return(list(df = df, df_quant = df_quant, df_cred = df_cred))
}

##' capelin_data---
#' function to switch between data types
#' @param capelin_data_set - data set: either biomass or abundance
#'
#' @return - returns the capelin data set desired for the analysis
#' @export
#'
#' @examples cap <- capelin_data(capelin_data_set)
capelin_data <- function(capelin_data_set){
# source: "Age disaggregate abundance for Keith Lewis - 2017 added_v1.xlsx"
     #browser()
     if(capelin_data_set == "biomass"){
          cap <- read_csv('data/capelin-2017.csv')
          cap$ln_abun_med <- log(cap$abundance_med)
          cap$ln_ab_lci <- log(cap$ab_lci)
          cap$ln_ab_uci <- log(cap$ab_uci)
          cap$ln_biomass_med <- log(cap$biomass_med)
          cap$ln_bm_lci <- log(cap$bm_lci)
          cap$ln_bm_uci <- log(cap$bm_uci)
          return(cap)
     }
     else if (capelin_data_set == "age") {
#from "capelin_age_disaggregate_abundance.xlsx" worksheet:Age disagg acoustic index
          capelinAbun <- read_csv('data/capelin_age_disaggregate_abundance.csv')
          str(capelinAbun)
          capelinAbun$age2_log <- log(capelinAbun$age2)
          capelinAbun$age2_log10 <- log10(capelinAbun$age2)
          capelinAbun$Ln_adult_abun <- log(sum(capelinAbun$age3, capelinAbun$age4, capelinAbun$age5, capelinAbun$age6, na.rm=T) + capelinAbun$age2*capelinAbun$age2PerMat)
          capelinAbun$adult_abun <- sum(capelinAbun$age3, capelinAbun$age4, capelinAbun$age5, capelinAbun$age6, na.rm=T) + capelinAbun$age2*capelinAbun$age2PerMat
          return(capelinAbun)
     }
}

##' condition_data---
#' function to switch between condition data types (index v residuals)
#' @param cond what type measure of condtion is being used, the fitted v. residuals (cond) or just the plain resids as provided by Fran
#' @param df name of the data set
#'
#' @return - the condition data set desired for the analysis
#' @export
#'
#' @examples - cond <- condition_data(cond, 'data/condition_ag1_MF_out.csv')
condition_data <- function(cond, df){
     browser()
     if(cond == "age"|cond == "biomass"){
          cond <- read_csv(paste(df))
          #lag data
          cond$meanCond_lag <- lag(cond$meanCond, 1)
          cond$medCond_lag <- lag(cond$medCond, 1)
          return(cond)

     } else if(cond == "resids") {
          #Fran's original
          condResids <- read_csv('data/archive/condition_resids.csv')
               condResids[19:20,2] <- NA
               condResids[c(19, 20), 1] <- matrix(c(2016, 2017), ncol = 1) 
               condResids$resids_lag <- lag(condResids$resids, 1)
               condResids$resids_lag[20] <- condResids$resids_lag[19]
               #View(condResids)
               return(condResids)
     }
}


##' make_direct---DEPRECATED - 
#' KEEP THIS just because the code may be useful
#' @param x a switch to start making folders
#'
#' @return
#' @export
#'
#' @examples
make_direct <- function(x){
     if(x==
        "yes"){
          if(!dir.exists("Bayesian/age2/mortality_0")) dir.create("Bayesian/age2/mortality_0")
     if(!dir.exists("Bayesian/age2/mortality_1")) dir.create("Bayesian/age2/mortality_1")
     if(!dir.exists("Bayesian/age2/mortality_2")) dir.create("Bayesian/age2/mortality_2")
     if(!dir.exists("Bayesian/age2/mortality_3")) dir.create("Bayesian/age2/mortality_3")
     if(!dir.exists("Bayesian/age2/recruitment_0")) dir.create("Bayesian/age2/recruitment_0")
     if(!dir.exists("Bayesian/age2/recruitment_1")) dir.create("Bayesian/age2/recruitment_1")
     if(!dir.exists("Bayesian/age2/recruitment_2")) dir.create("Bayesian/age2/recruitment_2")
     if(!dir.exists("Bayesian/age2/rm_1")) dir.create("Bayesian/age2/rm_1")
     if(!dir.exists("Bayesian/age2/rm_2")) dir.create("Bayesian/age2/rm_2")
     if(!dir.exists("Bayesian/age2/rm_3")) dir.create("Bayesian/age2/rm_3")
     if(!dir.exists("Bayesian/age2/rm3_1p_high")) dir.create("Bayesian/age2/rm3_1p_high")
     if(!dir.exists("Bayesian/age2/rm3_1p_med")) dir.create("Bayesian/age2/rm3_1p_med")     
     }
}

# a vector of folder names to feed "folder_path1" 
folder_names <- c("mortality_0", 
                  "mortality_1", 
                  "mortality_2", 
                 # "mortality_3", 
                  "recruitment_0", 
                  "recruitment_1", 
                  "recruitment_2", 
                  "rm_1", 
                  "rm_2", 
                  "rm_3",
                  "sensitivity"
                  #"rm3_1p_high",
                  #"rm3_1p_med"
                  )


#' make_direct1----
#' function that makes directories automatically if they don't exist
#' @param folder_names from vector above
#' @param folder_path_gen - folder path that indicates parent folder where child folders are to be created
#'
#' @return - a folder with sub-folders as per "folder_names"
#' @export
#'
#' @examples - make_direct1(folder_names, folder_path1)
make_direct1 <- function(folder_names, folder_path_gen){
     #browser()
     for(i in folder_names){
          if(!dir.exists(paste0(folder_path_gen, i))) dir.create(paste0(folder_path_gen, i))
     }
}



##' name_pairPlot()---
#' makes a smaller data frame with just the variables needed for the pairs-plot
#' 
#' @param x - dataframe
#' @param y - character object indicating "biomass" or "age"
#'
#' @return - returns a pairs plot used in EDA phase
#' @export
#'
#' @examples name_pairPlot(df1, "biomass")
name_pairPlot <- function(x, capelin_data_set){
     #browser()
     df_name <- x
     df_name <- subset(df_name, year >= 2003)
     if(capelin_data_set == "biomass"){
          names(df_name)[names(df_name) == "ln_biomass_med"] <- "ln capelin biomass"
          #names(df_name)[names(df_name) == "resids_lag"] <- "condition l1"
          names(df_name)[names(df_name) == "meanCond_lag"] <- "condition l1"
          
          names(df_name)[names(df_name) == "surface_tows_lag2"] <- "larval abundance l2"
          names(df_name)[names(df_name) == "ps_meanTot_lag2"] <- "zooplankton abun l2"
          return(df_name)
     } else if (capelin_data_set == "age"){
          names(df_name)[names(df_name) == "age2"] <-  "log10_age2 index"
          #names(df_name)[names(df_name) == "resids_lag"] <- "condition l1"
          names(df_name)[names(df_name) == "meanCond_lag"] <- "condition l1"
          
          names(df_name)[names(df_name) == "surface_tows_lag2"] <- "larval abundance l2"
          names(df_name)[names(df_name) == "ps_meanTot_lag2"] <- "zooplankton abun l2"
          return(df_name)
     }
     
}


# DEPRECATED - SEE PlotCredInt2
plotCredInt1 <- function(df, yaxis = yaxis, ylab = ylab, y_line = y_line, ci = ci, dpi = dpi, model = model, x = x, y = y){
     p <- ggplot()  
     #browser()
     p <- p + geom_ribbon(aes(x = c(df$year, 2015), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5, fill = "grey60")
     pi_n <- dpi[, ncol(dpi)]
     p <- p + geom_errorbar(aes(x = insert$year, 
                                ymax = pi_n[2], 
                                ymin = pi_n[1]), fill = "grey40")
     #p <- p + geom_text() +
     # annotate("text", label = model, x = x, y = y, size = 7.5)
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5)
     p <- p + geom_errorbar(data = df, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     p <- p + xlab("Year") + ylab(paste(ylab))
     # p <- p + theme(axis.title = element_text(size=30), axis.text = element_text(size = 20)) 
     #p <- p + theme(text = element_text(size=25)) + theme_bw()
     p <- p + geom_line(aes(x = c(df$year, 2015), y = y_line))
     p <- p + geom_point(data = insert, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5,
                         colour = "red", 
                         adj = 0.5)
     p <- p + geom_errorbar(data = insert, width = 0.3, colour = "red", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     
     p <- p + theme_bw() + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
     return(p)
}


##' plotCredInt2---
#'  This is a crude approach to leave one out
#' @param df - the data frame
#' @param insert - this is the value of the year to insert into the graph
#' @param yaxis - variable in df for yaxis 
#' @param ylab - customize the label of the yaxis
#' @param y_line - predicted values, e.g. y_med
#' @param ci - credible interval 
#' @param dpp - median value of the predition interval
#' @param dpi - 2.5% and 97.5% prediction interval
#' @param insert_year - value of year to align with credible interval
#' @param type - this indicates if the capelin index has a confidence interval: 
#'
#' @return - a graph of teh capelin biomass with 95% CIs, the fitted values of the model, and the 95% credible and prediction intervals.
#' @export
#'
#' @examplesplotCredInt2(df3, insert, yaxis = yaxis1, ylab = ylab1, y_line = y_med, ci = ci_df3, dpp = p_med, dpi = pi_df3, insert_year = insert_year, type = NA)

plotCredInt2 <- function(df, insert, yaxis = yaxis, ylab = ylab, y_line = y_line, ci = ci, dpp = dpp, dpi = dpi, insert_year = x, type=type){
     p <- ggplot()  
     browser()
     # plot credible interval
     p <- p + geom_ribbon(aes(x = c(df$year, insert_year), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5, fill = "grey60")
     pi_n <- dpi[, ncol(dpi)]
# predition interval for point
     p <- p + geom_errorbar(aes(x = insert$year, 
                                ymax = pi_n[2], 
                                ymin = pi_n[1]))
     p <- p + geom_point(aes(x = insert$year, 
                                y = last(dpp)),
                         colour = "black",
                         shape = 17, size = 2)
     
# ci for capelin
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5,
                         colour = "black")
     
     if(!is.na(type)){
          p <- p + geom_errorbar(data = df, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     } else if (is.na(type)) {
          p
     }
     
     p <- p + xlab("Year") + ylab(paste(ylab))
     p <- p + geom_line(aes(x = c(df$year, insert_year), y = y_line))
# ci for knockout year
     #browser()
     pd1 <- position_dodge(width = 2)
     p <- p + geom_point(data = insert, # insert is df3[x,]
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5,
                         colour = "red", 
                         position=pd1)
     if(!is.na(type)){
     p <- p + geom_errorbar(data = insert, width = 0.3, colour = "red", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     } else if (is.na(type)) {
          p
     }
     p <- p + theme_bw() + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
     return(p)
}


##' plotCredInt3---
#'  This is a crude approach to leave one out
#' @param df - the data frame
#' @param insert - this is the value of the year to insert into the graph
#' @param yaxis - variable in df for yaxis 
#' @param ylab - customize the label of the yaxis
#' @param y_line - predicted values, e.g. y_med
#' @param ci - credible interval 
#' @param dpp - median value of the predition interval
#' @param dpi - 2.5% and 97.5% prediction interval
#' @param insert_year - value of year to align with credible interval
#' @param type - this indicates if the capelin index has a confidence interval: 
#'
#' @return - a graph of teh capelin biomass with 95% CIs, the fitted values of the model, and the 95% credible and prediction intervals.
#' @export
#'
#' @examplesplotCredInt2(df3, insert, yaxis = yaxis1, ylab = ylab1, y_line = y_med, ci = ci_df3, dpp = p_med, dpi = pi_df3, insert_year = insert_year, type = NA)

plotCredInt3 <- function(df, insert, yaxis = yaxis, ylab = ylab, y_line = y_line, ci = ci, dpp = dpp, dpi = dpi, insert_year = x, type=type){
     p <- ggplot()  
     #browser()
     # plot credible interval
     p <- p + geom_ribbon(aes(x = c(df$year, insert_year), 
                              ymax = ci[2, ], 
                              ymin = ci[1, ]),
                          alpha = 0.5, fill = "grey60")
     pi_n <- dpi[, ncol(dpi)]
     # predition interval for point
     p <- p + geom_errorbar(aes(x = insert$year, 
                                ymax = pi_n[4], 
                                ymin = pi_n[1]))
     p <- p + geom_point(aes(x = insert$year, 
                             y = last(dpp)),
                         colour = "black",
                         shape = 17, size = 2)
     
     # ci for capelin
     p <- p + geom_point(data = df, 
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5,
                         colour = "black")
     
     if(!is.na(type)){
          p <- p + geom_errorbar(data = df, width = 0.3, colour = "black", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     } else if (is.na(type)) {
          p
     }
     
     p <- p + xlab("Year") + ylab(paste(ylab))
     p <- p + geom_line(aes(x = c(df$year, insert_year), y = y_line))
     # ci for knockout year
     #browser()
     pd1 <- position_dodge(width = 2)
     p <- p + geom_point(data = insert, # insert is df3[x,]
                         aes_string(y = yaxis, x = "year"),
                         shape = 16, 
                         size = 1.5,
                         colour = "red", 
                         position=pd1)
     if(!is.na(type)){
          p <- p + geom_errorbar(data = insert, width = 0.3, colour = "red", aes(x = year, min=ln_bm_lci, ymax=ln_bm_uci))
     } else if (is.na(type)) {
          p
     }
     p <- p + theme_bw() + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"))
     return(p)
}

insert_year <- c(2003:2017)     

# not using - not appropriate
DIC_out <- function(df){
     browser()
     dic_mod <-dic.samples(df$model, 1000, "pD")
     dic_val <- sum(dic_mod$deviance) + sum(dic_mod$penalty)
     return(dic_val)
     
}

#' rsq_bayes()----
#'
#' @param ypred - vector values of y from Bayesian chain 
#' @param out - list of output values from Bayesian analysis
#'
#' @return - Bayesian R-squared
#' @export
#'
#' @examples R1_r2 <- rsq_bayes(ypred = y_pred, out = run_recruit)
# see jags_bayes_R2 for proof that this matches the Gelmen code
rsq_bayes <- function(ypred = ypred, out=out){
          #browser()
          # variance of predicted values
          y_var = apply(ypred, 1, 'var')

          # variance of residuals
          res_pred = out$BUGSoutput$sims.list$Res
          #res_med = apply(res_pred,2,'median')
          res_var = apply(res_pred,1,'var')
          
          # Distribution of Bayesian R-squared
          bay_rsq <- y_var/(res_var + y_var)
          return(bay_rsq)
     }
