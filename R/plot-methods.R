# Default colour scheme for plot
#
# It's from colorbrewer. I use this fixed thing so I
# do not need to add a dep on RColorBrewer, but users are
# encouraged to make use of it in the doc.
#
# Currently teh default scheme is defined by:
# library(RcolorBrewer)
# brewer.pal(n=5, name ="Set1")
#
.defaultSpectraColours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

#' @title Plots an object inheriting from the Spectra class
#' @name plot-Spectra
#' @rdname plot-Spectra
#' @description The philosophy of this plotting routine is to provide a "quick'n'dirty" way
#' to plot your spectra collection. For advanced visualisations, the use of
#' \code{\link{melt_spectra}} alongside with ggplot2 or lattice is encouraged.
#' 
#' @aliases plot plot.Spectra plot,Spectra,ANY-method 
#' 
#' @usage \method{plot}{Spectra}(x,gg,gaps,attr,\dots)
#' 
#' @param x an object of class \code{"Spectra"} or inheriting from this class
#' @param gg if TRUE, uses the \code{ggplot2} package to plot the data, if
#' FALSE uses \code{matplot} from base graphics (much faster)
#' @param gaps if TRUE, gaps in the spectra are not plotted
#' @param attr attribute against which lines are coloured (only for \code{gg =
#' TRUE}
#' @param ... additional parameters passed to \code{matplot}
#' @author Pierre Roudier \url{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Default plotting method
#' plot(australia[1:5,])
#' 
#' # Default plot using ggplot2
#' plot(australia[1:5,], gg = TRUE)
#' 
#' \dontrun{
#' 
#' # Managing gaps in the spectra
#' s <- cut(australia, wl =c(-1*450:500, -1*1800:2050))
#' plot(s, gaps = TRUE)
#' plot(s, gaps = FALSE)
#' 
#' # passing various options to matplot
#' plot(
#'   australia, 
#'   lty = 1:5, 
#'   col = 'blue', 
#'   xlab = 'foo', ylab = 'bar', 
#'   ylim = c(0.4,0.6), 
#'   main = 'my plot'
#' )
#' 
#' # Using colour ramps
#' plot(
#'   australia, 
#'   lty = 1, 
#'   col = rainbow(10), 
#'   main = "It is possible to create really ugly visualisations"
#' )
#' 
#' # Example using colours given by ColorBrewer (http://colorbrewer2.org/)
#' library(RColorBrewer)
#' plot(australia, lty = 1, col = brewer.pal(n = 8, name = "Set2"))
#' 
#' # Using an attribute to group spectra
#' 
#' # Generate some kind of factor
#' australia$fact <- sample(LETTERS[1:3], size = nrow(australia), replace = TRUE) 
#' 
#' s <- aggregate_spectra(australia, fun = mean, id = 'fact')
#' plot(s, gg = TRUE, attr = 'fact')
#' }
NULL

#' @rdname plot-Spectra
#' @export plot
plot.Spectra <- function(x, gg = FALSE, gaps = TRUE, attr = NULL, ...) .plot.Spectra(x, gg, gaps, attr, ...)

.plot.Spectra <- function(x, gg = FALSE, gaps = TRUE, attr = NULL, ...){

  # Show gaps in the data?
  if (gaps) {
    x <- fill_spectra(x, fill = NA, ...)
  }

  if (gg) {
    
    if (is.null(attr)) s.melt <- melt_spectra(x)
    else s.melt <- melt_spectra(x, attr = attr)

    # force id colname
    names(s.melt)[1] <- 'id'

    p <- ggplot(s.melt) 

    if (is.null(attr)) {
      p <- p + geom_line(aes_string(x = 'wl', y = 'nir', group = 'id'))
    }
    else {
      p <- p + geom_line(aes_string(x = 'wl', y = 'nir', group = 'id', colour = attr))
    }
    p <- p +
      labs(x = paste("Wavelength (", wl_units(x), ")", sep = ""), y = "Value") +
      theme_bw()
    return(p)
  }
  else {
    # Fast implenmentation using matplot

    # inspect dots to check matplot args: type, lty, ylab, xlab, ylim
    dots <- list(...)
    nm_dts <- names(dots)

    # insert default values if no matplot args given by user
    if (!("type" %in% nm_dts)) dots$type <- 'l'
    if (!("lty" %in% nm_dts)) dots$lty <- 1
    if (!("ylab" %in% nm_dts)) dots$ylab <- "Value"
    if (!("xlab" %in% nm_dts)) dots$xlab <- paste("Wavelength (", wl_units(x), ")", sep = "")
    if (!("xlim" %in% nm_dts)) dots$xlim <- range(wl(x))
    if (!("ylim" %in% nm_dts)) dots$ylim <- range(spectra(x), na.rm = TRUE, finite = TRUE)
    if (!("col" %in% nm_dts)) dots$col <- .defaultSpectraColours
    
    # insert x and y values 
    dots$x <- wl(x)
    dots$y <- t(spectra(x))
    
    do.call("matplot", dots)
  }
}

if (!isGeneric("plot_summary")) {
    setGeneric("plot_summary", function(x, fun = mean, se = TRUE, ...)
        standardGeneric("plot_summary"))
}

#' @title Summary plot of a collection of spectra
#' @aliases plot_summary plot_summary,Spectra-method
#' @name plot_summary
#' @description Creates a summary plot of a collection of Spectra
#' @usage plot_summary(x, fun = mean, se = TRUE, ...)
#' @param x an object of class \code{"Spectra"} or inheriting from this class
#' @param fun an aggregation function
#' @param se if TRUE, plots the standard deviation around the summary spectra
#' (computed by function as given by \code{fun}). If FALSE, does not plot
#' dispersion information. If a function, uses this function instead of
#' \code{sd}.
#' @param ... additional parameters, currently ignored
#' @author Pierre Roudier
#' @examples 
#'  oz <- load_oz()
#'  plot_summary(oz)
setMethod("plot_summary", signature('Spectra'), 
  function(x, fun = mean, se = TRUE, ...) {

        # If sd is given as TRUE or FALSE
    if (is.logical(se)) {
      if (se) {
        plot.se <- TRUE
        fun.se <- sd
      }
      else {
        plot.se <- FALSE
      }
    }
    # If sd is being given a function
    else {
      # If the function is valid
      if (is.function(se)) {
        plot.se <- TRUE
        fun.se <- se
      }
      # Else stop
      else {
        stop('The se = ... option must evaluate to either logical or function.')
      }
    }

    s.melt <- melt_spectra(x)
    
  #   s.summary <- ddply(s.melt, 'wl', fun, ...)
    s.summary <- ddply(s.melt, 'wl', function(x) {do.call(fun, list(x$nir))})
    names(s.summary)[2] <- 'nir'

    if (plot.se) {
      # initiate dummy vars to pas R CMD check
      wl <- nir <- nir_se <- NULL

      s.se <- ddply(s.melt, 'wl', function(x) {do.call(fun.se, list(x$nir))})
      names(s.se)[2] <- 'nir_se'
      
      s <- join(s.summary, s.se, by = 'wl')

      p <- ggplot() + 
        geom_line(data = s, aes(x = wl, y = nir - nir_se), linetype = 2) +
        geom_line(data = s, aes(x = wl, y = nir + nir_se), linetype = 2) +
        geom_line(data = s, aes(x = wl, y = nir)) +
        labs(x = paste("Wavelength (", wl_units(x), ")", sep = ""), y = "Value") +
        theme_bw()
    } 
    else {
      p <- ggplot() + 
        geom_line(data = s.summary, aes(x = wl, y = nir)) +
        labs(x = paste("Wavelength (", wl_units(x), ")", sep = ""), y = "Value") +
        theme_bw()
    }

    p
  }
)

if (!isGeneric("plot_stack")) {
    setGeneric("plot_stack", function(x)
        standardGeneric("plot_stack"))
}

#' @title Stacked plot of a collection of spectra
#' @name plot_stack
#' @aliases plot_stack plot_stack,Spectra-method
#' @description Creates a stacked plot of a collection of Spectra
#' @usage plot_stack(x)
#' @param x an object of class \code{"Spectra"} or inheriting from this class
#' @author Pierre Roudier
#' @examples 
#'  oz <- load_oz(3)
#'  plot_stack(oz)
setMethod("plot_stack", signature('Spectra'), 
  function(x){
    
    m <- melt_spectra(x)
    idnm <- names(m)[1]
    m[[idnm]] <- as.factor(m[[idnm]])
    form_grid <- as.formula(paste(idnm, '~.'))
    ggplot(m) + 
      geom_line(aes_string(x = 'wl', y = 'nir', colour = idnm)) + 
      facet_grid(form_grid) + 
      labs(x = paste("Wavelength (", wl_units(x), ")", sep = ""), y = "Value") +
      theme_bw()
  }
)

if (!isGeneric("plot_offset")) {
    setGeneric("plot_offset", function(x, offset = 1)
        standardGeneric("plot_offset"))
}

#' @title Offset plot of a collection of spectra
#' @name plot_offset
#' @aliases plot_offset plot_offset,Spectra-method
#' @description Creates an offset plot of a collection of Spectra
#' @usage plot_offset(x, offset = 1)
#' @param x an object of class \code{"Spectra"} or inheriting from this class
#' @param offset Offset between spectra
#' @author Pierre Roudier
#' @examples 
#'  oz <- load_oz(3)
#'  plot_offset(oz)
#'  plot_offset(oz, 0.3)
#'  plot_offset(oz, 2)
setMethod("plot_offset", signature('Spectra'), 
  function(x, offset = 1){
    
    # offsets values
    offsets <- (seq_len(nrow(x)) - 1) * offset
    
    # affect spectra with offset values
    spectra(x) <- aaply(1:length(offsets), 1, function(i) {
      spectra(x)[i,] + offsets[i]
    })
    
    # melt spectra for data visualisation with ggplot2
    m <- melt_spectra(x)
    
    idnm <- names(m)[1]
    m[[idnm]] <- as.factor(m[[idnm]])
    
    ggplot(m) + 
      geom_line(aes_string(x = 'wl', y = 'nir', colour = idnm, group = idnm)) +
      labs(x = paste("Wavelength (", wl_units(x), ")", sep = ""), y = "") +
      theme_bw()

  }
)

## Code for adding NAs to potentially removed WLs
##
# ref reference wavelengths
# fill value to fill missing WLs with
#



if (!isGeneric("fill_spectra")) {
  setGeneric("fill_spectra", function(obj, ...)
    standardGeneric("fill_spectra"))
}

.fill_spectra <- function(obj, ref = NULL, fill = NA, ...) {

  if (is.null(ref)) {
    # Trying to get the most common resolution values
    r <- as.numeric(names(which.max(table(diff(wl(obj))))))
    nb_gaps <- length(table(diff(wl(obj))))
    if (nb_gaps > 2)
      warning("Sorry, at this stage removing gaps does not work well with irreguarly spaced wavelengths. Results might be odd.")
    ref <- seq(from = min(wl(obj)), to = max(wl(obj)), by = r)
  }
  
  # Detect missing WLs
  missing_wl <- setdiff(ref, wl(obj))

  # If tehre is gaps in the data, we add these as NAs values
  if (length(missing_wl) > 0) {
    # Create matrix of NAs for the missing WLs
    new_nir <- matrix(fill, ncol = length(missing_wl), nrow = nrow(obj))
    colnames(new_nir) <- missing_wl
    
    # Collate the NA matrix with the rest of the spectra
    new_nir <- cbind(spectra(obj), new_nir)
    # Re-order the spectra matrix
    idx <- order(as.numeric(colnames(new_nir)))
    new_nir <- new_nir[, idx, drop = FALSE]
    
    spectra(obj) <- new_nir
  }

  obj
}

#' @title Fill missing wavelengths of a Spectra* object with a given value
#' @name fill_spectra
#' @aliases fill_spectra fill_spectra,Spectra-method
#' @description Fill missing wavelengths of a Spectra* object with a given value. This is
#' mostly usefull to include NA values in the spectra in order to show missing
#' bits in plots.
#' 
#' @details At this stage removing gaps does not work well with irreguarly spaced
#' wavelengths. Results might be odd for binned spectra.
#' 
#' @param obj an object inheriting from class \code{Spectra}
#' @param ref a numeric vector, giving the reference wavelengths (ie the entire
#' range of wavelengths expected to be in the spectra before some waveleng5ths
#' have been cut out). If NULL, the function is trying to guess it.
#' @param fill values to fill gaps in the data with
#' @param ... ignored
#' @return An object of the same class as \code{obj}
#' @author Pierre Roudier \url{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Cut wavelengths out of the collection
#' oz <- cut(australia, wl=-1*c(355:400, 2480:2499))
#' big.head(spectra(oz), , 7)
#' 
#' # Giving the wavelengths at which I want data
#' oz_filled <- fill_spectra(oz, ref = 350:2500, fill = NA)
#' big.head(spectra(oz_filled), , 7)
#' plot(oz_filled)
#' 
#' # Trying to guess ref values
#' oz_filled <- fill_spectra(oz, fill = -999)
#' big.head(spectra(oz_filled), , 7)
#' plot(oz_filled)
#' 
#' @export fill_spectra
setMethod("fill_spectra", signature('Spectra'), .fill_spectra)