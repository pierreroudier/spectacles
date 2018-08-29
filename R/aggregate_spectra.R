
if (!isGeneric("aggregate_spectra"))
  setGeneric("aggregate_spectra", function(obj, fun = mean, ...)
    standardGeneric("aggregate_spectra"))

#' @title Aggregates spectral and data information
#' @name aggregate_spectra
#' @aliases aggregate_spectra,Spectra-method aggregate_spectra,SpectraDataFrame-method aggregate_spectra
#' @docType methods
#' @description Aggregates spectral and data information of a \code{Spectra} object using a
#' user-defined function
#' 
#' @details There is two distinct function for \code{Spectra} and \code{SpectraDataFrame} classes. For \code{SpectraDataFrame} objects, associated data is also aggregated using the function provided by the \code{fun} option. Additionally, the method for \code{SpectraDataFrame} has an \code{id} option that allows to specify an attribute which will be used to split the object, apply sequentially the \code{fun} function, and recombine the results in an unique object.
#' @param obj see below
#' @param fun see below
#' @param id see below
#' @param ... see below
#' @section Methods: \describe{
#' 
#' \bold{x=Spectra}
#' 
#' \code{aggregate_spectra(obj, fun=mean, ...)}
#' 
#' \tabular{rll}{ \tab \code{obj} \tab A \code{Spectra} object \cr \tab
#' \code{fun} \tab An aggregation function \cr \tab \code{...} \tab Expressions
#' evaluated in the context of \code{fun} \cr }
#' 
#' \bold{x=SpectraDataFrame}
#' 
#' \code{aggregate_spectra(obj, fun=mean, id=NULL, ...)}
#' 
#' \tabular{rll}{ \tab \code{obj} \tab A \code{SpectraDataFrame} object \cr
#' \tab \code{fun} \tab An aggregation function \cr \tab \code{id} \tab
#' Attribute(s) to split the object (character vector) \cr \tab \code{...} \tab
#' Expressions evaluated in the context of \code{fun} \cr }
#' 
#' }
#' @return An object of the same class as \code{obj}
#' @author Pierre Roudier \url{pierre.roudier@@gmail.com}
#' @seealso \code{\link{apply_spectra}}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Aggregation on the whole collection
#' m <- aggregate_spectra(australia, fun = mean)
#' plot(m)
#' 
#' # Aggregation factor-wise
#' 
#' # Generate some kind of factor
#' australia$fact <- sample(LETTERS[1:3], size = nrow(australia), replace = TRUE)
#' m <- aggregate_spectra(australia, fun = mean, id = 'fact')
#' plot(m)
setMethod("aggregate_spectra", "Spectra",
  function(obj, fun = mean, ...){
    
    # making up an id name from the aggregation function
    id_fun <- as.character(substitute(fun, env = parent.frame()))[1]
    id_obj <- as.character(substitute(obj, env = parent.frame()))
    id <- paste(id_fun, id_obj, sep = '.')
  
    # applying the function to the spectra
    nir <- aaply(.data = spectra(obj), .margins = 2, .fun = fun, ...)

    # Create and return Spectra object
    Spectra(wl = wl(obj), nir = nir, id = id, units = wl_units(obj))
  }
)

# In the case of a SDF, an id can be given to split the SDF and apply fun
#' @rdname aggregate_spectra
setMethod("aggregate_spectra", "SpectraDataFrame",
  function(obj, fun = mean, id = NULL, ...){
    
    # No split --> the whole data is aggregated together
    if (is.null(id)) {
      # making up an id name from the aggregation function
      id_fun <- as.character(substitute(fun, env = parent.frame()))[1]
      id_obj <- as.character(substitute(obj, env = parent.frame()))
      
      # Select and paste only alphanumeric chars
      id_obj <- paste(id_obj[grep(x = id_obj, pattern = '[[:alnum:]]')], collapse = '.')
      # Combine object name and function name into an id
      id <- paste(id_fun, id_obj, sep = '.')
  
      # applying the function to the spectra
      nir <- apply(spectra(obj), 2, fun, ...)
      
      res <- Spectra(wl = wl(obj), nir = nir, id = id, units = wl_units(obj))
      
      data <- sapply(features(obj), fun, ...)
            
      res <- SpectraDataFrame(res, data = data.frame(matrix(data, nrow = 1, dimnames = list(id, names(obj)))))
    }
    
    # There is a variable against which the data will be aggregated
    else {
     if (id %in% names(features(obj))) {
                
        # Col index of the splitting variable
        idx <- paste("(names(features(obj)) == '", id, "')", sep = "")
        idx <- paste(idx, collapse="|")
        idx <- which(eval(parse(text=idx)))
                
        # Creating spectra splits
        s <- data.frame(id = features(obj)[, idx, drop = FALSE], spectra(obj))
        colnames(s) <- gsub("id.", "", colnames(s))
        s <- na.omit(ddply(s, id, colwise(fun), ...))
              
        # Remove id used to split data.frame
        s <- s[, -(1:(length(id)))]
                
        # new data slot
        feat <- merge(select(features(obj), as.vector(id)), (select_if(features(obj), is.numeric)), by="row.names",  suffixes = "")
        rownames(feat) <- feat[,1]
        feat <- feat[,-1]
        d <- ddply(feat, id, colwise(fun), .drop=TRUE, ...)
        d <- na.omit(d[colSums(!is.na(d)) > 0])
        
        # recompose the object
        res <- SpectraDataFrame(wl = wl(obj), nir = s, units = wl_units(obj), data = d)
      }
      else
        stop('Bad aggregation identifier.')
    }

    res
  }
)
