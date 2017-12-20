postResampleSpectro <- function(pred, obs){
  # If preds and obs are void
  if (length(obs) + length(pred) == 0) {
    out <- rep(NA, 2)
  } else {
    if (length(unique(pred)) < 2 || length(unique(obs)) < 2) {
      resamplCor <- NA
    }
    else {
      # Compute preds-obs correlation
      resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
      # Manage error
      if (class(resamplCor) == "try-error") 
        resamplCor <- NA
    }
    
    n <- length(obs)
    
    # Compute R2
    r2 <- resamplCor^2
    
    # Compute MSE/SEP2
    mse <- mean((pred - obs)^2)
    
    # Standard error of prediction
    # (SEP/RMSEP)
    rmsep <- sqrt(mse)
    
    # Bias
    bias <- mean(pred) - mean(obs)
    
    # Standard error
    #     se <- sd(pred - obs)
    sse <- sum((pred - obs)^2)
    se <- sqrt(sse/(n - 1))
    
    # Residual  variance
    #     sep2c <- sqrt(sum(((pred - bias - obs)^2) / n))
    
    # ratio of performance to deviation
    rpd <- sd(obs)/rmsep
    
    # Ratio of performance to interquartile distance
    qs <- quantile(obs, probs = seq(0, 1, 0.25), names = FALSE)
    iq <- qs[4] - qs[2]
    rpiq <- iq/rmsep
    
    # Lin's CCC
    ccc <- as.numeric(epi.ccc(obs, pred)$rho.c[1])
    
    out <- c(rmsep, r2, rpd, rpiq, ccc, bias, se)
  }
  
  # Manage var names
  names(out) <- c("RMSE", "Rsquared", "RPD", "RPIQ", "CCC", "Bias", "SE")
  
  # Manage NAs
  if (any(is.nan(out))) 
    out[is.nan(out)] <- NA
  out  
}

#' @title Calculates performance indictors across resamples
#' @name postResampleSpectro
#' @aliases postResampleSpectro spectroSummary
#' @description Given two numeric vectors of data, the root mean squared error, the R-squared, the bias, the RPD, the RPIQ, the CCC and the standard error are calculated. For two factors, the overall agreement rate and Kappa are determined. 
#' @usage 
#' postResampleSpectro(pred, obs)
#' spectroSummary(data, lev = NULL, model = NULL)
#' @param pred A vector of numeric data
#' @param obs A vector of numeric data
#' @param data a data frame or matrix with columns obs and pred for the observed and predicted outcomes
#' @param lev a character vector of factors levels for the response. In regression cases, this would be NULL.
#' @param model a character string for the model name 
#' @details This function extends \code{postResample} in the \code{caret} package.
#' @author Pierre Roudier, adapted from code from Max Kuhn
#' @examples
#' predicted <-  matrix(rnorm(50), ncol = 5)
#' observed <- rnorm(10)
#' apply(predicted, 2, postResampleSpectro, obs = observed)
spectroSummary <- function (data, lev = NULL, model = NULL) {
  if (is.character(data$obs)) 
    data$obs <- factor(data$obs, levels = lev)
  postResampleSpectro(data[, "pred"], data[, "obs"])
}