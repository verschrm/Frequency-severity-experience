##### Credentials #####

### Author:       R.M. (Robert) Verschuren
### Institution:  Amsterdam School of Economics, 
###               University of Amsterdam
### Date:         17/08/2022

##### Description #####

### This function creates the ordered Lorenz curve
### based on a reference (benchmark) model and a 
### candidate (alternative) model, and computes 
### its Gini index with the corresponding variance.
### The function accepts the following arguments:

# Candidate  = The risk premia for the candidate model;
# Reference  = The risk premia for the reference model;
# Actual     = Actual observed losses;
# Figure     = Whether to plot the (ordered) Lorenz curve.

##### Function rGini() #####
rGini <- function(Candidate, Reference, Actual, Figure = FALSE) {
  ### Make sure the required packages are loaded
  Packages <- c('ggplot2')
  invisible(lapply(Packages, library, character.only = TRUE))
  
  ### Calculate and order the relativities
  ORel <- as.matrix(cbind(Reference, Actual, Candidate / Reference))
  ORel <- ORel[order(ORel[, 3], decreasing = FALSE), ]
  
  ### Determine the empirical cumulative distribution function
  H <- length(Actual)
  F_H <- ecdf(ORel[, 3])
  ORel_H <- as.matrix(F_H(ORel[, 3]))
  Grid_ORel <- as.matrix(unique(ORel_H))
  Rep <- sapply(Grid_ORel, function(z) length(which(ORel_H <= z)))
  Rep <- diff(c(0, Rep))
  
  ### Compute the ordered premium and loss distribution
  F_P <- sapply(cumsum(Rep), function(z) sum(ORel[1:z, 1])) / sum(ORel[, 1])
  F_P <- rep(F_P, Rep)
  F_L <- sapply(cumsum(Rep), function(z) sum(ORel[1:z, 2])) / sum(ORel[, 2])
  F_L <- rep(F_L, Rep)
  
  ### Calculate the empirical Gini index
  Gini <- 1 - (F_P - c(0, F_P[-H])) %*% (F_L + c(0, F_L[-H]))
  
  ### Determine an estimate of the variance of the Gini index
  mu_L <- sum(Actual) / H
  mu_P <- sum(Reference) / H
  h <- (mu_L * Reference * F_L + Actual * mu_P * (1 - F_P)) / 2
  mu_h <- sum(h) / H
  Sigma_L <- ((Actual - mu_L) %*% (Actual - mu_L)) / H
  Sigma_P <- ((Reference - mu_P) %*% (Reference - mu_P)) / H
  Sigma_h <- ((h - mu_h) %*% (h - mu_h)) / H
  Sigma_hL <- ((h - mu_h) %*% (Actual - mu_L)) / H
  Sigma_hP <- ((h - mu_h) %*% (Reference - mu_P)) / H
  Sigma_LP <- ((Actual - mu_L) %*% (Reference - mu_P)) / H
  Gini <- c(Gini, 
            4 / (mu_L^2 * mu_P^2) * (4 * Sigma_h + mu_h^2 / mu_L^2 * Sigma_L + 
            mu_h^2 / mu_P^2 * Sigma_P - 4 * mu_h / mu_L * Sigma_hL -
            4 * mu_h / mu_P * Sigma_hP + 2 * mu_h^2 / (mu_L * mu_P) * Sigma_LP) / H)
  
  ### Store output
  Output <- list(Gini) 
  names(Output) <- c('Gini')
  
  ### Plot the ordered Lorenz curve
  if (Figure) {
    df_fig <- data.frame(F_P, F_L)
    labs <- ifelse(rep(all.equal(as.numeric(Reference), rep(1, length(Reference))), 2) == TRUE, 
                   c('Policy distribution', 'Lorenz curve'), c('Premium distribution', 'Ordered Lorenz curve'))
    colnames(df_fig) <- c('P', 'L')
    ### Construct the ordered Lorenz curve
    Lorenz <- ggplot(df_fig, aes(x = P, y = L)) + 
                     geom_line(aes(x = P, y = P, color = rgb(0, 0, 0, alpha = 0.6)), size = 0.7) +
                     geom_line(aes(color = rgb(248/256, 118/256, 109/256, alpha = 0.80)), size = 1) +
                     xlab(labs[1]) + ylab('Loss distribution') + theme_classic() + 
                     scale_y_continuous(labels = scales::percent, breaks = c(seq(0.00, 1.00, by = 0.20))) + 
                     scale_x_continuous(labels = scales::percent, breaks = c(seq(0.00, 1.00, by = 0.20))) +
                     theme(plot.title = element_text(size = 12, face = 'bold', hjust = 0.5), 
                     legend.position = c(0.12, 0.92), axis.title.x = element_text(size = 9, face = 'italic'), 
                     axis.title.y = element_text(size = 9, face = 'italic'), strip.text.x = element_blank()) +
                     theme(panel.grid.major = element_line(color = rgb(0, 0, 0, alpha = 0.1), size = 0.5)) + 
                     scale_color_manual(NULL, values = c(rgb(0, 0, 0, alpha = 0.6),
                     rgb(248/256, 118/256, 109/256, alpha = 0.80)), 
                     labels = c('Line of equality', labs[2]),
                     guide = guide_legend(override.aes = list(linetype = c('solid',  'solid'), pch = 
                     c(NA, NA), size = c(1,  1), fill = c(rgb(0, 0, 0, alpha = 0.6),
                     rgb(248/256, 118/256, 109/256, alpha = 0.80)))))
    
    ### Store output
    Output <- list(Output, Lorenz)
    names(Output) <- c(names(Output)[1], 'Lorenz')
  }
  
  ### Return output  
  return(Output)
}
