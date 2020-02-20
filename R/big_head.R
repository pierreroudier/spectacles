#' Return the First or Last Part of an Object
#' 
#' Return the First or Last Part of an Object
#' 
#' Returns the first or last rows of a data frame like head() and tail(), but
#' also only returns the first and last columns. This has been implemented to
#' check big data frames.
#' 
#' @aliases big.head big.tail big.head big.tail
#' @param x a \code{"data.frame"} or a \code{"matrix"} object
#' @param n a single, positive integer, number of rows for the object to return
#' @param l a single, positive integer, the number of columns to include on the left
#' @param r a single, positive integer, the number of columns to include on the right
#' @return An object (usually) like 'x' but generally smaller.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{\link{head}}, \code{\link{tail}}
#' @examples
#' 
#' big.head(mtcars)
#' big.tail(mtcars)
#' big.tail(mtcars, 10)
#' big.head(mtcars, 10, 2, 4)
#' big.head(mtcars, , , 1)
#' 
#' data(australia)
#' big.head(australia)
#' 
#' @export big.head
big.head <- function(x, n = 5, l = 5, r = 5){
  
  # n, l, and r must be of length 1
  stopifnot(length(n) == 1L)
  stopifnot(length(l) == 1L)
  stopifnot(length(r) == 1L)
  
  # Correcting n, l, and r if they are bigger than the number of rows
  # or columns
  n <- min(nrow(x), n)
  l <- min(ncol(x) - 1, l)
  r <- min(ncol(x), r)
  
  # Get index of columns on left
  idx_l <- seq_len(l)
    
  # Get index of columns on right
  idx_r <- seq(ncol(x) - r + 1, ncol(x))
  
  # Get index of columns represented by dots
  idx_dots <- setdiff(1:ncol(x), union(idx_l, idx_r))
  
  # if no dots, there might be an overlap between l and r
  if (length(idx_dots) == 0) {
    # No dots necessary, just return the cropped rows
    res <- x[seq_len(n), , drop = FALSE]
  } else {
    # We need to include dots fill-up
    x1 <- x[seq_len(n), seq_len(l), drop = FALSE]
    xdots <- rep('...', length.out = n)
    x2 <- x[seq_len(n), seq(ncol(x) - r + 1, ncol(x)), drop=FALSE]
    res <- data.frame(x1, xdots, x2)
    names(res)[l + 1] <- "..."
  }
  
  res
}

big.tail <- function(x, n = 5, l = 5, r = 5){
  
  # n, l, and r must be of length 1
  stopifnot(length(n) == 1L)
  stopifnot(length(l) == 1L)
  stopifnot(length(r) == 1L)
  
  # Correcting n, l, and r if they are bigger than the number of rows
  # or columns
  n <- min(nrow(x), n)
  l <- min(ncol(x) - 1, l)
  r <- min(ncol(x), r)
  
  # Get index of columns on left
  idx_l <- seq_len(l)
  
  # Get index of columns on right
  idx_r <- seq(ncol(x) - r + 1, ncol(x))
  
  # Get index of columns represented by dots
  idx_dots <- setdiff(1:ncol(x), union(idx_l, idx_r))
  
  # if no dots, there might be an overlap between l and r
  if (length(idx_dots) == 0) {
    # No dots necessary, just return the cropped rows
    res <- x[seq.int(to = nrow(x), length.out = n), , drop = FALSE]
  } else {
    # We need to include dots fill-up
    x1 <- x[seq.int(to = nrow(x), length.out = n), seq_len(l), drop = FALSE]
    xdots <- rep('...', length.out = n)
    x2 <- x[seq.int(to = nrow(x), length.out = n), seq(ncol(x) - r + 1, ncol(x)), drop=FALSE]
    res <- data.frame(x1, xdots, x2)
    names(res)[l + 1] <- "..."
  }
  
  res
}
