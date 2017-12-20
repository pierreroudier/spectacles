#' @include SpectraDataFrame-methods.R

#' @title Australia spectra library data set
#' @name australia
#' @aliases australia oz
#' 
#' @description The \code{australia} data set gathers 100 soil spectra, along side with their organic
#' carbon, pH and clay content values. This data set has been collected by
#' CSIRO.The \code{oz} dataset is the first 5 samples of \code{australia}, and is here to keep examples fast to run.
#' 
#' The \code{data.frame} contains the following columns: 
#' \itemize{
#'  \item{sr_no}{a unique identifier for each spectrum}
#'  \item{carbon}{soil organic carbon values}
#'  \item{ph}{soil pH values}
#'  \item{clay}{soil clay values}
#'  \item{X350, X351, \dots{},X2499, X2500}{reflectance in wavelengths 350 to 2500nm}
#' }
#' 
#' @docType data
#' @author Data kindly contributed by Raphael Viscarra Rossel. To reduce the
#' size of the package, the original dataset has been reduced to 100 spectra
#' using the Kennard-Stone algorithm.
#' @references Viscarra Rossel, R. and Behrens, T. 2010.
#' @keywords datasets
#' @examples
#' 
#' data(australia)
#' big.head(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' data(oz)
#' big.head(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' 
NULL

#' @title Load the australia dataset
#' @name load_oz
#' @description Loads the australia dataset as a SpectraDataFrame 
#' @param n the number of spectra to return. By default, it returns all 100 spectra in the dataset.
#' @return a SpectraDataFrame
#' @author Pierre Roudier
#' @examples 
#' oz <- load_oz()
#' 
load_oz <- function(n = NULL) {
  # Load data as SPC
  data("australia", envir = environment())
  
  if (!is.null(n)) {
    stopifnot(n > 0)
    australia <- australia[sample(1:nrow(australia), size = n),]
  }
  
  spectra(australia) <- sr_no ~ ... ~ 350:2500
  
  australia
}
