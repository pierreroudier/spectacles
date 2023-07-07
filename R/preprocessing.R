#' @title Apply a function on the spectra of a Spectra* object
#' @name apply_spectra
#' @description Aggregates spectral and data information of a \code{Spectra} object using a
#' user-defined function.
#' 
#' Apply a function and update the spectra of a \code{Spectra} object. This
#' function is particularly interesting for pre-processing, e.g to derive the
#' spectra, or apply pre-processing functions such as \code{snv}.
#' 
#' @usage apply_spectra(obj, fun, ...)
#' @details The philosophy of this function is to let the user free to use any function
#' to pre-process a spectra collection, using either functions from the stats
#' package, functions from other packages such as \code{signal}, or personal
#' functions.
#' 
#' @param obj an object inheriting from class \code{Spectra}
#' @param fun an aggregation function
#' @param ... expressions evaluated in the context of \code{fun}
#' @return An object of the same class as \code{obj}
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{\link{aggregate_spectra}}, \code{\link{snv}},
#' \code{\link{rnv}}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Second derivative
#' r <- apply_spectra(australia, diff, 2)
#' plot(r)
#' 
#' # Smoothing kernel
#' k <- kernel("daniell", 20) # define a kernel
#' r <- apply_spectra(australia, kernapply, k)
#' plot(r)
#' 
#' \dontrun{
#' # Savitzky-Golay filter (from the signal package)
#' library(signal)
#' r <- apply_spectra(australia, sgolayfilt, n = 33, p = 4)
#' plot(r)
#' }
#' 
#' @export apply_spectra
apply_spectra <- function(obj, fun, ...) {
  nir <- aaply(spectra(obj), 1, fun, ...)

  # Managing case where only one spectra
  if (nrow(obj) == 1)
    nir <- matrix(nir, nrow = 1, dimnames = list(ids(obj), names(nir)))

  spectra(obj) <- nir
  obj
}

#' @title Standard and Robust Normal Variate transformations
#' @name snv
#' @aliases snv rnv
#' @description Standard and Robust Normal Variate transformations are often used in
#' chemometrics to normalise a spectra collection and remove the baseline
#' effect.
#' 
#' The Standard Normal Variate transformation (SNV, Barnes et al., 1989) is a
#' common method to reduce within-class variance.
#' 
#' The Robust Normal Variate transformation (RNV, Guo et al., 1999) is a
#' modification of the SNV to make it more robust to closure problems.
#' 
#' These function are to be used in conjonction with
#' \code{apply_spectra}.
#' 
#' @param x a vector of numeric values
#' @param r the percentile to use in the RNV computation
#' @return A vector of numeric values
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @references 
#' \itemize{ 
#'   \item Barnes, R.J., Dhanoa, M.S., Lister, S.J. 1989.
#' Standard normal variate transformation and detrending of near-infra-red
#' diffuse reflectance spectra. Applied Spectroscopy 43, 772--777.
#'   \item Guo, Q., Wu, W., Massar, D.L. 1999. The robust normal variate
#' transform for pattern recognition with near-infrared data. Analytica Chimica
#' Acta 382:1--2, 87--103.  
#' }
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Standard Normal Variate transform
#' s <- apply_spectra(australia[1:10,], snv)
#' plot(s)
#' 
#' # The scale function in the base package is actually doing
#' # the same thing!
#' s <- apply_spectra(australia[1:10,], scale, center = TRUE, scale = TRUE)
#' plot(s)
#' 
#' # Robust Normal Variate transform
#' s <- apply_spectra(australia[1:10,], rnv, r = 0.25)
#' plot(s)
#' 
#' @export snv rnv
snv <- function(x){
  (x - mean(x))/sd(x)
}

#' @rdname snv
rnv <- function(x, r){
  pct <- as.numeric(quantile(x = x, probs = r, na.rm = TRUE))
  (x - pct)/(sd(x[x <= pct]))
}


## Baseline using the baseline package

# if (!isGeneric('base_line'))
  setGeneric('base_line', function(object, ...)
    standardGeneric('base_line')
)

#' @title Baseline correction using the baseline package#' 
#' @name base_line
#' @aliases base_line base_line,Spectra-method
#' @docType methods
#' @description Estimates baselines for the spectra in the \code{obj} object, using the
#' algorithm named in 'method'.
#' 
#' @details The baseline package implements various algorithms for the baseline
#' correction. The following methods are available:
#' 
#' \itemize{ \item 'als': Baseline correction by 2nd derivative constrained
#' weighted regression \item 'fillPeaks': An iterative algorithm using
#' suppression of baseline by means in local windows \item 'irls' (default): An
#' algorithm with primary smoothing and repeated baseline suppressions and
#' regressions with 2nd derivative constraint \item 'lowpass': An algorithm for
#' removing baselines based on Fast Fourier Transform filtering \item
#' 'medianWindow': An implementation and extention of Mark S. Friedrichs'
#' model-free algorithm \item 'modpolyfit': An implementation of Chad A. Lieber
#' and Anita Mahadevan-Jansen's algorithm for polynomial fiting \item
#' 'peakDetection': A translation from Kevin R. Coombes et al.'s MATLAB code
#' for detecting peaks and removing baselines \item 'rfbaseline': Wrapper for
#' Andreas F. Ruckstuhl, Matthew P. Jacobson, Robert W. Field, James A. Dodd's
#' algorithm based on LOWESS and weighted regression \item 'rollingBall': Ideas
#' from Rolling Ball algorithm for X-ray spectra by M.A.Kneen and H.J.
#' Annegarn. Variable window width has been left out }
#' 
#' See baseline package documentation for more information and references.
#' 
#' Additionally, the baseline package provides a nice GUI that helps choosing
#' the good baseline method and the good parametrisation. This GUI can be used
#' with the \code{inspectr} package. This is demonstrate in the Examples
#' section.
#' 
#' @param object an object inheriting from class \code{Spectra}
#' @param ... additional arguments to be passed to the \code{baseline} function
#' in the baseline package. The main option would be \code{'method'}, to switch
#' between the several baseline methods presented in teh details section.
#' @return An object of the same class as \code{obj} with the continuum removed
#' from its spectra.
#' @author Interface to the baseline package by Pierre Roudier
#' \email{pierre.roudier@@gmail.com}, baseline package authored by Kristian Hovde
#' Liland and Bjorn-Helge Mevik
#' @seealso \code{\link{continuum_removal}}, \code{\link{snv}},
#' \code{\link{rnv}}
#' @references Kristian Hovde Liland and Bjorn-Helge Mevik (2011). baseline:
#' Baseline Correction of Spectra. R package version 1.0-1.
#' http://CRAN.R-project.org/package=baseline
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Subsample for demo purposes
#' australia <- australia[1:10,]
#' 
#' # Correction using the default method (irls)
#' bl <- base_line(australia)
#' plot(bl)
#' 
#' # Specifying another method for baseline calculation
#' bl2 <- base_line(australia, method = "modpolyfit")
#' plot(bl2)
#' 
#' # Using the baseline package independently
#' # (useful to plot the corrections)
#' \dontrun{
#' library(baseline)
#' bl3 <- baseline(spectra(australia), method = 'irls')
#' class(bl3) # this is a baseline object
#' plot(bl3)
#' # Affecting the baseline-corrected spectra back
#' # to the SpectraDataFrame object
#' spectra(australia) <- getCorrected(bl3)
#' plot(australia)
#' 
#' # Using the baselineGUI with inspectr
#' baselineGUI(spectra(australia))
#' ## When happy with a configuration, clik "Apply to all" and 
#' ## save the results under a name, e.g. "corrected.spectra"
#' spectra(australia) <- getCorrected(corrected.spectra)
#' plot(australia)
#' }
#' 
setMethod('base_line', 'Spectra', function(object, ...) {
  nir <- spectra(object)
  new_nir <- baseline(nir, ...)
  spectra(object) <- getCorrected(new_nir)
  object
})

#' Continuum removal
#' 
#' Operates a continuum removal on a vector.
#' 
#' This operation is commonly done to normalize reflectance spectra and allow
#' comparison of individual absorption features from a common baseline. The
#' removal is based on the upper convex hull of the spectra.
#' 
#' This function is working on vectors. It may applied on matrix or data.frames
#' using the \code{apply} function, or on \code{Spectra*} objects using the
#' \code{apply_spectra} function.
#' 
#' @param x a numeric vector
#' @param wl the wavelengths of the spectra
#' @param upper if TRUE, removes the upper convex hull from the spectra, if
#' FALSE, takes the lower convex hull
#' @return A numeric vector with its continuum removed.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}, based on code from
#' Raphael Viscarra-Rossel.
#' @seealso \code{\link{baseline}}, \code{\link{snv}}, \code{\link{rnv}}
#' @references Clark, R.N., and Roush, T.L. 1984. Reflectance spectroscopy:
#' Quantitative analysis techniques for remote sensing applications. Journal of
#' Geophysical Research 89, 6329--6340.
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' s <- apply_spectra(australia, continuum_removal)
#' plot(s)
#' 
#' s <- apply_spectra(australia, continuum_removal, upper = FALSE)
#' plot(s)
#' 
#' @export continuum_removal
continuum_removal <- function(x, wl = as.numeric(names(x)), upper = TRUE){
  
  if (!upper) x <- -1*x
  
  # Compute convex hull
  ch_idx <- chull(x = wl, y = x)
  # Close the polygon
  ch_idx <- c(ch_idx, ch_idx[1])
  # Put in a data.frame
  ch <- data.frame(wl = wl[ch_idx], nir = x[ch_idx])
  
  # We don't want the polygon, just the line above the spectra
  # We will select only the bit that goes from wl min to wl max
  
  idx_min <- which(ch$wl == min(ch$wl))[1]
  idx_max <- which(ch$wl == max(ch$wl))[1]
  
  if (idx_min > idx_max) {
    idx_select <- c(idx_min:nrow(ch), 1:idx_max)
  } else {
    idx_select <- idx_min:idx_max
  }
  ch <- ch[idx_select,]
  
  # Linear interpolation
  res <- approx(x = ch$wl, y = ch$nir, xout = wl)

  # Remove continuum
  x - res$y
}
