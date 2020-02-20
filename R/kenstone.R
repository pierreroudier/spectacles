if (!isGeneric("kenstone"))
  setGeneric("kenstone", function(x, size, ...)
    standardGeneric("kenstone"))

#' @title Kennard-Stone algorithm for optimal calibration set selection.
#' @name kenstone
#' @description An implemnentation of the Kennard-Stone algorithm for calibration set selection.
#' @aliases kenstone kenstone,Spectra-method
#' @usage kenstone(x, size, ...)
#' @param x a \code{Spectra} object
#' @param size a positive number, the number of items to choose from
#' @param ... ignored
#' @return A vector of length size giving the indices of the selected
#' individuals in x.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @references 
#' Kennard, L.A. Stone, Technometrics 11 (1969) 137.
setMethod("kenstone", "Spectra", 
  function(x, size, ...) {
    idx <- .kenstone.matrix(spectra(x), size = size, ...)
    idx
  }
)

## Kennard-Stone algorithm for calibration set selection
##
## -- use it on PCs to compress the data !! --
##
#
# 1. Select the object which is closest to the data mean
# and add it to the subset
#
# 2. Calculate the dissimilarity between the remaining objects
# in the data set and the objects in the subset
#
# 3. Select the object, which is the most dissimilar to the ones
# already included in the subset
#
# 4. return to 2 until the desired number of objects in the subset
# is reached
#
.kenstone.matrix <- function(x, size, progress = TRUE, ...){
  # Initialisation
  k <- 1
  sub <- vector(mode = 'numeric', length = size)
  sub_mat <- matrix(NA, ncol = ncol(x), nrow = size)
  
  # Compute data mean
  mn <- colMeans(x)
  # Get closest point to the mean
  idx <- which.min(.distToMean(x, mn))
  # Put point into subset and remove point from x
  sub[k] <- idx
  sub_mat[k,] <- x[idx,]
  k <- k + 1
  #   x <- x[-1*idx,]
  x[idx, ] <- NA
  
  # Get farest point from the first selected subset
  idx <- which.max(.distToMean(x, sub_mat[1, ]))
  # Put point into subset and remove point from x
  sub[k] <- idx
  sub_mat[k, ] <- x[idx, ]
  k <- k + 1
  #   x <- x[-1*idx,]
  x[idx, ] <- NA
  
  if (progress) {
    i_pb <- 1
    pb <- txtProgressBar(min = 1, max = size, style = 3)
  }      
  
  # Sequentially select points
  while(k <= size) {
    #
    # (1) compute the distance between each point in remaining in x
    # and each point already in the subset
    # (2) for each point in x, select the MINIMUM distance with each
    # point already in the subset
    # (3) the algorithm is selecting the point in x that have the 
    # BIGGEST minimum value (ie is most dissimilar to points already
    # in the subset)
    #
    # I know it sounds crazy but it goes better with a pen and paper
    #
    idx <- which.max(apply(.distToSubset(x, sub_mat[1:(k-1),]), 2, min))
    sub[k] <- idx
    sub_mat[k,] <- x[idx,]
    k <- k + 1
    #     x <- x[-1*idx,]
    x[idx, ] <- NA
    
    if (progress) {
      setTxtProgressBar(pb, k)
    }
  }
  
  if (progress)
    close(pb)
  
  sub
}

# Distance of each point of a matrix x to
# a unique value
.distToMean <- function(x, mn){
  apply(x, 1, function(x) dist(rbind(mn, x)))
}

# Distance of each point in a matrix x to
# each point in another matrix sub
.distToSubset <- function(x, sub) {
  apply(x, 1, function(y) apply(sub, 1, function(x) dist(rbind(x, y))))
}
