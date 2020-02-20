#' Internal helper function for the splice correction
#' Computes the splice correction matrix between the 
#' VIS and SWIR1 sensors
.correct_splice_1 <- function(w, nir, vertex, reference) {
  
  # Getting index of wavelengths
  idx_vtx <- which(w == vertex)
  idx_ref <- which(w == reference)
  # List all wavelengths between vertex and reference
  wls <- w[idx_vtx:idx_ref]
  
  # Process each row of the spectra matrix
  res <- laply(1:nrow(nir), function(i) {
    correction_top_left <- (wls - vertex + 1)^2 
    correction_top_right <- nir[i, idx_ref + 1] - nir[i, idx_ref]
    correction_top <- correction_top_left * correction_top_right
    
    correction_bottom <- nir[i, idx_ref] * (reference - vertex + 1)^2
    
    1 + correction_top / correction_bottom
  })
  
  res
}

#' Internal helper function for the splice correction
#' Computes the splice correction matrix between the 
#' SWIR1 and SWIR2 sensors
.correct_splice_2 <- function(w, nir, vertex, reference) {
  
  # Getting index of wavelengths
  idx_vtx <- which(w == vertex)
  idx_ref <- which(w == reference)
  
  # List all wavelengths between vertex and reference
  wls <- w[idx_ref:idx_vtx]
  
  # Process each row of the spectra matrix
  res <- laply(1:nrow(nir), function(i) {
    correction_top_left <- (wls - vertex)^2 
    correction_top_right <- nir[i, idx_ref] - nir[i, idx_ref + 1]
    correction_top <- correction_top_left * correction_top_right
    
    correction_bottom <- nir[i, idx_ref + 1] * (reference - vertex)^2
    
    1 + correction_top / correction_bottom
  })
  
  res
}

# Method for a matrix
# w a vector of the wavelengths at which spectra are collected
# nir a matrix of spectra
# locations a list of length 2 describing the reference points and vertex points
# for the parabolic offset correction.
.splice.numeric <- function(w, nir, locations = list(c(700, 1000), c(1830, 1950))) {
  
  # Initialisation of the correction matrix as identity
  correction <- nir * 0 + 1
  
  # Correct left/first splice
  correction_1 <- .correct_splice_1(
    w = w, nir = nir, 
    vertex = min(locations[[1]]), 
    reference = max(locations[[1]])
  ) 
  idx_correction_1 <- which(
    w >= min(locations[[1]]) & w <= max(locations[[1]])
  )
  correction[, idx_correction_1] <- correction_1
  
  # Correct second/right splice
  correction_2 <- .correct_splice_2(
    w = w, nir = nir, 
    vertex = max(locations[[2]]), 
    reference = min(locations[[2]])
  ) 
  idx_correction_2 <- which(
    w >= min(locations[[2]]) & w <= max(locations[[2]])
  ) + 1
  correction[, idx_correction_2] <- correction_2
  
  nir * correction
}

# Method for a Spectra*
# s an object of class Spectra*
# locations a list of length 2 describing the reference points and vertex points
# for the parabolic offset correction.
.splice.Spectra <- function(s, locations = list(c(750, 1000), c(1830, 1950))) {
  res <- .splice.numeric(w = wl(s), nir = spectra(s), locations = locations)
  spectra(s) <- res
  return(s)
}

if (!isGeneric("splice"))
  setGeneric("splice", function(x, ...)
    standardGeneric("splice"))

#' @title Splice correction of a spectra collected using ASD hardware
#' @name splice 
#' @aliases splice,Spectra-method
#' @description This is the correction method available in the ViewSpec Pro
#' software from ASD, which aims at correcting steps in the data (see details).
#' 
#' @details The SWIR1 part of the spectrum (1000-1800 nm) is taken as a reference 
#' for corrections as it is stable to the instrument sensitivity drift (Beal and Eamon, 2010)
#' 
#' This is based on a description of the splice correction algorithm described in:
#' 
#' Beal, D. and Eamon, M., 1996. Dynamic, Parabolic Linear Transformations of 
#'   'Stepped' Radiometric Data. Analytical Spectral Devices Inc., Boulder, CO.
#' 
#' @param x a \code{Spectra} object
#' @param locations the wavelengths to cut out and interpolate
#' @return an object of same class as x
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' oz_spliced <- splice(australia)
#' plot(oz_spliced)
#' 
setMethod("splice", 
  signature(x = "Spectra"),
  function(x, locations = list(c(750, 1000), c(1830, 1950))) {
    .splice.Spectra(x, locations)
  }
)
