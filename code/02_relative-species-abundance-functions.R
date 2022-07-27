# Please note: This script is obtained/modified from: https://github.com/GotelliLab/FAMA/blob/1.0/Scripts/BetaDirFunctions.R
# This script was a supplement to Gotelli et al. 2021

# Source file of utility functions for working with beta and dirichlet distributions

# Preliminaries -----------------------------------------------------------
library(gtools) # for dirichlet
library(tidyverse)
library(betareg)
library(patchwork)
# set.seed(100) # for repeatable results

# -------------------------------------------------------------------------


# FUNCTION beta_stats -----------------------------------------------------
# inputs: win= number of successful trials
#         lose= number of unsuccessful trials
#         interval=span of confidence interval (0.95 is default)
# outputs: a list with the following elements
#          win: number of successful trials
#          lose: number of unsuccessful trials
#          freq_mean: frequentist expectation of probability of success
#          bayes_mean: bayesian expectation of probability of success
#          bayes_low: 95% lower confidence bound
#          bayes_high: 95% upper confidence bound
# -------------------------------------------------------------------------
beta_stats <- function(win=sample(x=0:10,,size=1),
                       lose=sample(x=0:10,size=1),
                       interval=0.95){
  
  p_low <- (1 - interval)/2
  p_high <- interval + p_low
  
  freq_mean <- ifelse(win + lose==0,NA,win/(win + lose))
  bayes_mean <- (win + 1)/(win + lose +2)
  bayes_low <- qbeta(p=p_low,shape1=win + 1, shape2=lose+1)
  bayes_high <- qbeta(p=p_high,shape1=win + 1, shape2=lose+1)
  
  beta_stats_out <- list(win=win,
                         lose=lose,
                         freq_mean=freq_mean,
                         bayes_mean=bayes_mean,
                         bayes_low=bayes_low,
                         bayes_high=bayes_high)
  
  return(beta_stats_out)
}

# -------------------------------------------------------------------------


# beta_stats Example ------------------------------------------------------
# beta_stats(0:3,0:3)
# beta_stats()


# FUNCTION dirch_stats ----------------------------------------------------
# input: x=a vector of counts of different species (may include zeroes)
#        optionally, x can be a named vector with species names
# output: dirch_sad = a dataframe with the following columns. Each row 
#         corresponds to the statistics for each species
#          species: species names, if present in input vector
#          freq_mean: frequentist expectation of probability of success
#          bayes_mean: bayesian expectation of probability of success
#          bayes_low: 95% lower confidence bound
#          bayes_high: 95% upper confidence bound
#
dirch_stats <- function(x=c(0,sample(1:100,10)),
                        interval=0.95,
                        replicates=10000){
  
  p_low <- (1 - interval)/2
  p_high <- interval + p_low
  
  # bootstrap the SAD; rows=reps, columns=species
  sad <- rdirichlet(n=replicates,alpha=x+1) 
  
  # convert matrix to data frame
  list_data <- as.data.frame((sad))          
  
  # little function to summarize output
  sum_fun <- function(z){list(bayes_mean=mean(z),
                              bayes_low=quantile(z,p_low),
                              bayes_high=quantile(z,p_high))}
  # map the function and bind the rows of output into a data frame
  df <- bind_rows(map(list_data,sum_fun))       
  
  # calculate the frequentist mean and include it as the first column
  z <- list(counts=x,freq_mean=x/sum(x))
  
  # if the input vector has species names, turn those into a tibble column
  df <- bind_cols(list(species=names(x)),z,df)
  return(df)
}

# -------------------------------------------------------------------------

# dirch_stats Example -----------------------------------------------------
# named_vec <- c(a=1,b=2,c=4)
# dirch_stats(x=named_vec)
# dirch_stats()

# -------------------------------------------------------------------------
# --------------------------------------
# FUNCTION dirch_reg
# description: regression stats from dirichlet estimates of FA
# inputs: x = field proportional abundance of each species
#         y = museum proportional abundance of each species
# outputs: after log-log transformation of both axes
#         intercept = regression intercept
#         slope = regression slope
#         p_val = p value
#         r2 = r squared
#         cutpoint = a/(1 - b)
########################################
dirch_reg <- function(data="MockData",
                      x = runif(100),
                      y = runif(100)) {
  my_model <- lm(log10(y) ~ log10(x))
  
  # test null hypothesis that slope = 1
  # see https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r
  my_model2 <- lm(log10(y) ~ log10(x) +
                    offset(1*log10(x)))
  
  
  slope1_test <- summary(my_model2)$coefficients[2,4]
  intercept <- summary(my_model)$coefficients[1,1]
  slope <-summary(my_model)$coefficients[2,1]
  p_val <- summary(my_model)$coefficients[2,4]
  r2 <- summary(my_model)$r.squared
  cutpoint <- 10^(intercept/(1 - slope))
  
  # find minimum and maximum field and museum proportions
  
  z <- which(x %in% min(x))
  field_min <- x[z]
  museum_min <- y[z]
  min_dif <-  - museum_min - field_min
  min_ratio <- museum_min/field_min
  
  z <- which(x %in% max(x))
  field_max <- x[z]
  museum_max <- y[z]
  max_dif <- field_max - museum_max
  max_ratio <- field_max/museum_max
  
  return(list(data=data,
              intercept=intercept,
              slope=slope,
              p_val=p_val,
              r2=r2,
              cutpoint=cutpoint,
              slope1_test=slope1_test,
              min_dif=min_dif,
              min_ratio=min_ratio,
              max_dif=max_dif,
              max_ratio=max_ratio))
  
} # end of dirch_reg
# --------------------------------------
dirch_reg()


# --------------------------------------
# FUNCTION dirch_plots
# description: Creates dirichet plots of FAMA data
# inputs: df_field and df_musee output data frames from dirch_stats, file name
# each data frame has columns for:
# species: species name
# freq_mean: frequentist mean proportions
# bayes_mean: dirichlet mean proportions
# bayes_low: dirichlet 95% low estimate
# bayes_high: dirichlet 95% high estimate
# outputs: output_description
########################################
dirch_plots <- function(df_field,
                        df_musee,
                        dataset="MockData") {
  
  # Combine data frames ------------------
  names(df_field) <- paste("field_",names(df_field),sep="")
  names(df_musee) <- paste("musee_",names(df_musee),sep="")
  df <- cbind(df_field,df_musee)  
  # create Figure 1 untransformed -----------------
  Figure1 <- ggplot(data=df,
                    aes(x=field_bayes_mean,y=musee_bayes_mean)) +
    geom_point() +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, linetype="solid") +
    # scale_y_log10() +
    geom_abline(slope=1, intercept=0,linetype="dashed") +
    xlab("Field RA") +
    ylab("Museum RA") +
    theme_bw(base_size=25) 
  print(Figure1)
  # Create Figure 2 log-log plot ------------------
  Figure2 <- ggplot(data=df,
                    aes(x=field_bayes_mean,y=musee_bayes_mean)) +
    
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 1) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, linetype="solid") +
    scale_y_log10() +
    geom_abline(slope=1, intercept=0,linetype="dashed") +
    scale_x_log10()+ 
    xlab("Field RA") +
    ylab("Museum RA") +
    labs(title=dataset) +
    theme_bw(base_size=25) +
    theme(plot.title = element_text(size=10)) 
  
  plot(Figure2)
  
  Figure2 <- Figure2 + labs(title=NULL)
  
  # create Figure 3 log-log vertical error bars ------------------
  Figure3 <- ggplot(data=df,
                    aes(x=field_bayes_mean,y=musee_bayes_mean)) +
    geom_segment(aes(x=field_bayes_mean, y=musee_bayes_low,
                     xend=field_bayes_mean,yend=musee_bayes_high))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 1) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, linetype="solid") +
    scale_y_log10() +
    geom_abline(slope=1, intercept=0,linetype="dashed") +
    scale_x_log10()+ 
    xlab("Field RA") +
    ylab("Museum RA") +
    theme_bw(base_size=25)
  print(Figure3)
  
  # Figure 4 log-log plot with vertical, horizontal error bars ------------------  
  # 
  
  Figure4 <- ggplot(data=df,
                    aes(x=field_bayes_mean,y=musee_bayes_mean)) +
    geom_segment(aes(x=field_bayes_mean, y=musee_bayes_low,
                     xend=field_bayes_mean,yend=musee_bayes_high))  +
    scale_y_log10() +
    geom_abline(slope=1, intercept=0,linetype="dashed") +
    scale_x_log10()+ 
    xlab("Field RA") +
    ylab("Museum RA") +
    geom_segment(aes(x=field_bayes_low, y=musee_bayes_mean, 
                     xend=field_bayes_high,yend=musee_bayes_mean)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 1) +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95, linetype="solid") +
    theme_bw(base_size=25) 
  print(Figure4)
  
  
  composite_plot <- (Figure1 | Figure2)/(Figure3 | Figure4) +
    plot_annotation(title=dataset, 
                    tag_levels = 'a',
                    tag_suffix = '.') & 
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))
  print(composite_plot)
  ggsave(plot=composite_plot, filename=paste("Output/",dataset,"_c.pdf",sep=""),width=11,height=8.55,units="in",device="pdf")  
  ggsave(plot=Figure2, filename=paste("Output/",dataset,".pdf",sep=""),width=11,height=8.55,units="in",device="pdf")  
  ggsave(plot=composite_plot, filename=paste("Output/",dataset,"_c.jpeg",sep=""),width=11,height=8.55,units="in",device="jpeg")  
  ggsave(plot=Figure2, filename=paste("Output/",dataset,".jpeg",sep=""),width=11,height=8.55,units="in",device="jpeg")  
} # end of dirch_plots
# --------------------------------------
# field <- c(0,0,5,10,20,100)
# names(field) <- letters[1:length(field)]
# musee <- c(1,2,4,6,10,20)
# names(musee) <- names(field)
# f <- dirch_stats(field)
# m <- dirch_stats(musee)
# dirch_plots(f,m)
# dirch_reg(x=f$bayes_mean,m$bayes_mean)

# FUNCTION beta_reg -------------------------------------------------------
# input: x_var= a continuous predictor variable (e.g., time or latitude)
#        y_var= a probability variable 0.0<y_var<1.0
#        weights= a continuous variable of regression weights (e.g. sample size)
# output: a list with the following elements:
#         slope= parameter from beta regression
#         p_val= statistical significance of slope
#         pseudo_r2= pseudo-r-squared from regression model
#         plot_reg= boolean variable to call ggplot
# if needed, other functions can be built from this
# that have two or more predictor variables

beta_reg <- function(x_var=1970:1986,
                     y_var=runif(length(x_var)),
                     weights=rep(1,length(x_var)),
                     plot_reg=FALSE){
  
  m <- betareg(y_var~x_var,weights=weights)
  
  
  slope <- summary(m)$coefficients$mean[2]
  p_val <- summary(m)$coefficients$mean[8]
  pseudo_r2 <- summary(m)$pseudo.r.squared
  my_stats <- list(slope=slope,p_val=p_val,psuedo_r2=pseudo_r2)
  
  # kick out optional fitted line
  if(plot_reg){
    print("show the plot")
    df <- data.frame(x_var,y_var)
    my_plot <- ggplot(df, aes(x = x_var, y = y_var)) +
      geom_point(size=3) +
      scale_fill_grey() +
      geom_line(aes(y = predict(m, df))) +
      theme_bw()
    print(my_plot) 
  }
  return(my_stats)
  
}
# -------------------------------------------------------------------------
# --------------------------------------
# FUNCTION share_stats
# description: gets unique and shared species
# inputs: two vectors of species counts
# outputs: unique S1, unique S2, sharedS
########################################
share_stats <- function(a=rbinom(n=10,size=1,
                                 prob=0.7),
                        b=rbinom(n=10,size=1,prob=0.7)){
  s1_unique <- sum(a>0 & b==0)
  s2_unique <- sum(b>0 & a==0)
  shared <- sum(a>0 & b>0)
  
  ## pull out abundance regression stats
  counts_reg <- lm(b~a)
  p_val <- summary(counts_reg)$coefficients[2,4]
  r2 <- summary(counts_reg)$r.squared
  return(list(s1=s1_unique,
              s2=s2_unique,
              shared=shared,
              p_val_counts=p_val,
              r2_counts=r2))

} # end of share_stats
# --------------------------------------
share_stats()
# beta_reg Example --------------------------------------------------------
#suppressWarnings(beta_reg(weights=runif(17),
#                           plot_reg=TRUE))

# -------------------------------------------------------------------------


# Toy Example of using functions ------------------------------------------
# This shows how to use the two functions to get a regression model
# for a single hexbin
# The functions in the code above this line should probably be set up as a separate script
# you would first source() that script,
# and then run the code below for the actual function



# first create some fake hexbin data:

# year <- seq(1970,1990,by=5) # 5 year time slices
# aliens <- c(0,3,10,50,5) # number of records of alien ant species
# natives <- c(10,50,30,100,10) # number of records of native ant species


# get beta stats on counts of natives and aliens in each year
# . <- beta_stats(win=natives,lose=aliens)
# print(.)

# bundle them up as a data frame and include date
# df <- as.data.frame(.)
# df <- cbind(year,df)
# print(df)


# initial unweighted run

# unweighted_run <- suppressWarnings(beta_reg(x_var=df$year, 
#                            y_var=df$bayes_mean,
#                            plot_reg=TRUE)) 
# print(unweighted_run)
# 
# now use inverse confidence interval as weights
# weighted_run <- suppressWarnings(beta_reg(x_var=df$year, 
#                            y_var=df$bayes_mean,
#                            weights=1/(df$bayes_high - df$bayes_low),
#                            plot_reg=TRUE)) 
# print(weighted_run)


# --------------------------------------
# --------------------------------------
# FUNCTION final_plot
# description: construct publication summary plot for each data set
# inputs: field RA, museum RA
# outputs: ggplot object
########################################
final_plot <- function(df=data.frame(ran_x=runif(10),ran_y=runif(10))) {
  
  # function body
  p_final <- ggplot(data=df) +
    aes(x=df[,1],y=df[,2]) +
    geom_point() +
    scale_y_log10() +
    geom_abline(slope=1, 
                intercept=0,
                linetype="dashed") +
    scale_x_log10()+
    labs(x = "Field RA", 
         y = "Museum RA") +
    geom_point(shape = 21, 
               colour = "black", 
               fill = "white", 
               size = 2, 
               stroke = 1) +
    geom_smooth(method="lm", 
                se=TRUE, 
                fullrange=FALSE, 
                level=0.95, 
                linetype="solid") +
    theme_bw(base_size=20) 
  return(p_final)
  
} # end of final_plot
# --------------------------------------
# function to create individual plots for final multipanel
fig_gen <- function(ob) {
  df <- readRDS(ob)
  f <- dirch_stats(df$x)[,3]
  m <- dirch_stats(df$y)[,3]
  return(final_plot(df=data.frame(f,m)))
}