#' @include Spectra-methods.R

#' @title Constructor for the SpectraDataFrame class.
#' @name SpectraDataFrame
#' @description Constructor for the SpectraDataFrame class. Creates a SpectraDataFrame
#' object, either from scratch, or from an existing Spectra object.
#' 
#' 
#' @param ... an object inheriting from \code{"Spectra"}
#' @param wl a numeric vector giving the wavelengths at with the spectra have
#' been measured
#' @param nir a \code{"matrix"} or a \code{"data.frame"} object giving the
#' spectra values for each sample
#' @param id a vector giving the unique id of each sample in the collection
#' @param units a character giving the unit in which the wavelengths values are
#' expressed
#' @param data object of class \code{"data.frame"} containing the attribute
#' data
#' @return a new \code{"SpectraDataFrame"} object
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{\link{spectra}}, \code{\link{wl}},
#' \code{\link{Spectra-class}}
#' @examples
#' 
#'   # Creating a SpectraDataFrame object from scratch
#'   my.wl <- 350:2500
#'   my.id <- c("A", "B")
#'   my.nir <- matrix(runif(2*length(my.wl)), nrow=2)
#'   my.data <- data.frame(foo = runif(2), bar = LETTERS[1:2])
#'   my.sdf <- SpectraDataFrame(wl = my.wl, nir = my.nir, id = my.id, data = my.data)
#' 
#'   # Creating a SpectraDataFrame object from an existing Spectra object
#'   my.s <- Spectra(wl = my.wl, nir = my.nir, id = my.id)
#'   my.sdf <- SpectraDataFrame(my.s, data = my.data)
#' 
#' @export SpectraDataFrame
"SpectraDataFrame" <- function(..., wl=numeric(), nir=matrix(), id=as.character(NA), units="nm", data=data.frame()) {

  dotargs <- list(...)

  # Initialisation from Spectra object(s)
  if (any(sapply(dotargs, inherits, "Spectra"))) {
    id_spectra <- which(sapply(dotargs, inherits, "Spectra"))
    # if there's more than one Spectra object
    if (length(id_spectra) > 1) {
      ss <- dotargs[id_spectra]
      s <- ss[[1]]
      for (i in 2:length(id_spectra))
	s <- rbind(s, ss[[i]])
    }
    # if theres only one Spectra object
    else
      s <- dotargs[[1]]

    wl <- wl(s)
    nir <- spectra(s)
    id <- ids(s, as.vector = FALSE)
    units <- wl_units(s)
  }

  else {
    # if the wl are given as an integer vector they are translated into a numeric vector
    # for clarity (only one type to manage)
    if (is(wl, "integer"))
      wl <- as.numeric(wl)
    if (is(nir, 'data.frame'))
      nir <- as.matrix(nir)
#     if (!is(id, "character"))
#       id <- as.character(id)
    if (!is(id, "data.frame"))
      id <- data.frame(id = id)

    # If no id is given
    if (all(is.na(id))) {
      # If the object is void
      if (length(nir) == 1)
        id <- data.frame(NULL)
      # if a matrix is here
      else
        id <- data.frame(id = as.character(seq(1, nrow(nir))))
    }
    else {
      # Test of inconsistent ids when id is specified by the user
      if (is.null(nrow(nir))) { # if theres only one spectra
        if (nrow(id) != 1)
          stop("number of individuals and number of rows in the spectra matrix don't match")
        if ((length(wl) > 1) & (length(nir) != length(wl)))
          stop("number of columns in the spectra matrix and number of observed wavelengths don't match")
        nir <- matrix(nir, nrow=1)
      }
      else {
        if (nrow(nir) != nrow(id))
          stop("number of individuals and number of rows in the spectra matrix don't match")
        if ((length(wl) > 1) & (ncol(nir) != length(wl)))
          stop("number of columns in the spectra matrix and number of observed wavelengths don't match")
        colnames(nir) <- wl
        rownames(nir) <- as.vector(do.call('rbind', id))
      }
    }
  }
  
  if (is(data, "numeric") | is(data, "integer"))
    data <- as.data.frame(data)

  # We only use rownames if there's no duplicated IDs. rownames are purely aesthetic.
  if (!any(duplicated(id[, 1]))) {
    rownames(data) <- id[, 1]
  } else {
    rownames(data) <- NULL
  }
  
  new("SpectraDataFrame", wl=wl, nir=nir, id=id, units=units, data=data)
}

## coercition methods

as.data.frame.SpectraDataFrame = function(x, ..., expand = TRUE, exclude_id = FALSE)  {
  
  data <- features(x, exclude_id = exclude_id)

  if (expand) {
    df <- data.frame(data, spectra(x))
    names(df) <- c(names(data), wl(x))
  }
  else {
    df <- data.frame(data, NIR = I(spectra(x)))
  }

  df
}

setAs("SpectraDataFrame", "data.frame", function(from) {
	as.data.frame.SpectraDataFrame(from)
  })

## Getting the data

# if (!isGeneric("features"))
  setGeneric("features", function(object, exclude_id = TRUE)
  # setGeneric("features", function(object, ...)
    standardGeneric("features"))

#' @title Retrieves or sets the data slot of a SpectraDataFrame object.
#' @name features
#' @aliases features features<- features-methods 
#' features,SpectraDataFrame-method features<-,Spectra-method
#' features<-,SpectraDataFrame-method
#' @docType methods
#' @description Either retrieves the attributes values from the data slot of a
#' SpectraDataFrame object, or upgrades a Spectra object to a SpectraDataFrame
#' object by initialising its data slot by a suitable \code{"data.frame"}
#' object.
#' @usage 
#' \S4method{features}{SpectraDataFrame}(object,exclude_id)
#' \S4method{features}{Spectra}(object,safe,exclude_id,key,append) <- value
#' 
#' @param object a \code{Spectra} object
#' @param value see below
#' @param safe see below
#' @param key see below
#' @param exclude_id see below
#' @param append see below
#' @return The \code{features} methods return a \code{data.frame} object, while
#' the \code{"features<-"} methods return a \code{SpectraDataFrame} object.
#' @section Methods: \describe{
#' 
#' \bold{x=Spectra}
#' 
#' \code{features(object, safe=TRUE, key=NULL, exclude_id=TRUE) <- value}
#' 
#' \tabular{rll}{ 
#' 
#' \tab \code{object} \tab A \code{Spectra} object \cr \tab
#' \code{safe} \tab Logical. If TRUE, data is being added to the object using a
#' SQL join (using a key field  given by the \code{key} option), otherwise it is
#' assumed the order of the rows is consitent with the order of the rows in
#' \code{object} \cr \tab \code{key} \tab Character, name of the column of the
#' data.frame storing the ids for the SQL join. Ignored if \code{safe} is
#' \code{FALSE}. \cr \tab \code{exclude_id} \tab Logical, if \code{TRUE}, ids
#' used for the SQL join are removed from the data slot after the join.\cr 
#' }
#' 
#' \bold{x=SpectraDataFrame}
#' 
#' \code{features(obj, exclude_id=TRUE)}
#' 
#' \code{features(obj, safe=TRUE, key=NULL, exclude_id=TRUE, append=TRUE) <-
#' value}
#' 
#' \tabular{rll}{ \tab \code{object} \tab A \code{SpectraDataFrame} object \cr
#' \tab \code{safe} \tab Logical. If TRUE, data is being added to the object
#' using a SQL join (using a key field given by the \code{key} option),
#' otherwise it is assumed the order of the rows is consitent with the order of
#' the rows in \code{object} \cr \tab \code{key} \tab Character, name of the
#' column of the data.frame storing the ids for the SQL join. Ignored if
#' \code{safe} is \code{FALSE}. \cr \tab \code{exclude_id} \tab Logical. For
#' the \code{features} method, if \code{TRUE}, the spectra ids are added to the
#' \code{data.frame} that is returned. For the \code{"features<-"} method, If
#' \code{TRUE}, ids used for the SQL join are removed from the data slot after
#' the join. \cr \tab \code{append} \tab Logical, if \code{TRUE}, the data is
#' appended to any existing data. if FALSE, the data provided is erasing any
#' existing data. \cr }
#' 
#' }
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{\link{spectra}}, \code{\link{wl}},
#' \code{\link{SpectraDataFrame-class}}
#' @examples
#' 
#' # Loading example data
#' data(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' 
#' # Printing available data
#' features(oz)
#' 
#' # Promoting a Spectra to a SpectraDataFrame object
#' s <- as(oz, "Spectra")
#' 
#' # Generating dummy data
#' d <- data.frame(
#'   id = ids(oz), 
#'   foo = runif(nrow(oz)), 
#'   bar = sample(LETTERS[1:5], size = nrow(oz), replace = TRUE)
#' )
#' head(d)
#' 
#' # Affecting data to Spectra object
#' features(s) <- d
#' features(s)
#' 
#' # Adding data to an existing SpectraDataFrame object
#' 
#' features(oz) <- d
#' features(oz)
#' 
setMethod("features", "SpectraDataFrame",
  function(object, exclude_id = TRUE) {
  # function(object, ..., exclude_id = TRUE) {
    if (!exclude_id) {
      res <- data.frame(ids(object, as.vector = FALSE), object@data) 
    } else {
      res <- object@data
    }
    res
  }
)

## Append or replace data

# if (!isGeneric('features<-'))
  # setGeneric('features<-', function(object, value, safe = TRUE, key = NULL, exclude_id = TRUE, append = TRUE)
  setGeneric('features<-', function(object, value)
    standardGeneric('features<-')
)

setReplaceMethod(
  "features", 
  signature("Spectra", "ANY"),
  function(object, value) {
    SpectraDataFrame(object, data = value)
  }
)
                   
                   
setReplaceMethod("features", signature("SpectraDataFrame", "ANY"),
  # safe enables id check
  # key gives the column name of the ids in the data.frame
  function(object, value) {
  # function(object, value, safe = TRUE, key = NULL, exclude_id = TRUE, append = TRUE) {
    
    if (!inherits(value, "data.frame"))
      stop('data must be provided as a data.frame object')

    # Check whether the number of records in the Spectra object and the data.frame are the same
    
    # if (safe) {
    #   if (is.null(key))
    #     stop("In the safe mode, you need to provide either the column name of the sample ids to the key option.")
    #   if (length(key) != 1)
    #     stop("Please provide only ONE id column.")
    # 
    #   if (is.numeric(key)) {
    #     key <- names(value)[key]
    #   }
    # 
    #   if (append) {
    #     # Actual ID sanity check
    #     d <- data.frame(ids(object, as.vector = FALSE), features(object))
    #     # Using the "key" name for ids
    #     names(d)[1] <- key
    #   }
    #   else {
    #     d <- ids(object, as.vector = FALSE)
    #     # Using the "key" name for ids
    #     names(d) <- key
    #   }
    # 
    #   # Safety: to avoid headaches with factors,
    #   # Both id columns are forced to be characters
    #   d[[key]] <- as.character(d[[key]])
    #   value[[key]] <- as.character(value[[key]])
    # 
    #   # Put data together
    #   data <- join(d, value,  by = key, type = "left", match = "first")
    #   # removing the id column      
    #   if (exclude_id)
    #     data <- data[, -1*which(names(data) == key), drop = FALSE]
    # }
    # else {
      warning("Sample ID check has been disabled. This mode assumes you made sure the order of the rows in your data is consistent with the order in which these samples appear in the Spectra object.")
      
      # if (append) data <- data.frame(features(object), value)
      # else data <- value
      
    data <- cbind(features(object), value)
    # }
    
    SpectraDataFrame(object, data = data)
  }
)

setMethod("$", "SpectraDataFrame",
  definition=function(x, name) x@data[[name]]
)

setReplaceMethod("$", "Spectra",
  definition=function(x, name, value) {
    # For SpectraDataFrame
    if ('data' %in% slotNames(x)) {
      x@data[[name]] <- value
    }
    # Else promoting Spectra to SpectraDataFrame
    else {
      data <- data.frame(value)
      names(data) <- name
      x <- SpectraDataFrame(x, data = data)
    }
    x
  }
)

setMethod("[[", c("SpectraDataFrame", "ANY", "missing"),
  function(x, i, j, ...) {
    if (!("data" %in% slotNames(x)))
      stop("no [[ method for object without attributes")
    x@data[[i]]
  }
)

setReplaceMethod("[[", c("Spectra", "ANY", "missing", "ANY"),
  function(x, i, j, value) {
    # For SpectraDataFrame
    if ('data' %in% slotNames(x)) {
      x@data[[i]] <- value
    }
    # Else promoting Spectra to SpectraDataFrame
    else {
      data <- data.frame(value)
      names(data) <- i
      x <- SpectraDataFrame(x, data = data)
    }
    x
  }
)

setMethod("[", c("SpectraDataFrame", "ANY", "ANY", "missing"),
  function(x, i, j, ..., k, drop = FALSE) {
    
    .bracket <- function(x, i, j, k, ..., drop = FALSE) {

      missing.i <- missing(i)
      missing.j <- missing(j)
      missing.k <- missing(k)

      # ROWS
      if (missing.i) {
        i <- TRUE
      } 
      else {
        # throws an error if trying to index rows using NAs
        if (any(is.na(i))) {
          stop("NAs not permitted in row index")
        }
        # in the case indexing rows by ids
        if (is.character(i)) {
          i <- which(x@id %in% i)
        }
      }

      # COLS
      if (missing.j) {
        j <- TRUE
      } 
      else {
        if (is.numeric(j)) {
          # If the indices are all negative, cols are removed
          if (all(j < 0)) {
            j <- setdiff(1:ncol(x), abs(j))
          }
        }
        else {
          j <- which(names(x) %in% j)
        }
      }

      # WAVELENGTHS
      if (missing.k) {
        k <- TRUE
      } 
      else {
        # If the indices are all negative, cols are removed
        if (all(k < 0)) {
          k <- setdiff(as.numeric(wl(x)), abs(k))
        }
        k <- which(as.numeric(wl(x)) %in% k)
      }

      SpectraDataFrame(wl = x@wl[k], nir = x@nir[i, k, drop = FALSE], id = x@id[i, , drop = FALSE], data = features(x)[i, j, drop = FALSE])
    }
  
  .bracket(x, i, j, k, ..., drop = drop)
  }
)

names.SpectraDataFrame <- function(x) names(x@data)

"names<-.SpectraDataFrame" <- function(x, value) {
  names(x@data) <- value
  x
}

if (!isGeneric("melt_spectra"))
  setGeneric("melt_spectra", function(obj, attr, ...)
    standardGeneric("melt_spectra"))

#' @rdname melt_spectra
setMethod("melt_spectra", "SpectraDataFrame", function(obj, attr = NULL, ...){
  
  id.nm <- names(ids(obj, as.vector = FALSE))

  if (!is.null(attr)) {
    data <- subset(features(obj), select = attr)
    x <- data.frame(ids(obj, as.vector = FALSE), data, spectra(obj))
    names(x) <- c(id.nm, attr, wl(obj))
  }
  else {
    x <- data.frame(ids(obj, as.vector = FALSE), spectra(obj))
    names(x) <- c(id.nm, wl(obj))
  }
  
  res <- melt(x, id.vars = c(id.nm, attr), variable.name = 'wl', value.name = "nir")
  res$wl <- as.numeric(as.character(res$wl))
  res
})

.subset.SpectraDataFrame <- function(x, subset, select, drop = FALSE, ...) {
  # adapted from subset.data.frame
  df <- features(x)
  if (missing(subset))
        r <- TRUE
  else {
    e <- substitute(subset)
    r <- eval(e, df, parent.frame())
    if (!is.logical(r))
	stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
  }
  if (missing(select))
    vars <- TRUE
  else {
    nl <- as.list(seq_along(df))
    names(nl) <- names(df)
    vars <- eval(substitute(select), nl, parent.frame())
  }
  df_sub <- df[r, vars, drop = drop]
  # remove unused factors
  df_sub <- droplevels(df_sub)
  id_selected <- which(rownames(df) %in% rownames(df_sub))
  SpectraDataFrame(wl = wl(x), nir = spectra(x)[id_selected, , drop = FALSE], id = ids(x, as.vector = FALSE)[id_selected, 1, drop = FALSE], units = wl_units(x), data = df_sub)
}

#' @title Subset SpectraDataFrame object
#' @name subset
#' @description Returns subsets of a SpectraDataFrame object.
#' @aliases subset subset.SpectraDataFrame subset,SpectraDataFrame-method
#' @docType methods
#' @usage subset(x, ...)
#' @param x SpectraDataFrame object
#' @param ... see details below
#' @details 
#' Additional parameters:
#' \describe{
#'   \item{subset}{logical expression indicating elements or rows to keep: missing values are taken as false}
#'   \item{select}{expression, indicating columns to select from the data slot}
#'   \item{drop}{passed on to "[" indexing operator}
#'   \item{...}{Additional arguments}
#' }
#' @return SpectraDataFrame object
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{\link{mutate}}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Subset on attributes
#' s <- subset(australia, carbon < 5)
#' 
#' # Subset and selection of attributes
#' s <- subset(australia, carbon < 5, select = 1)
setMethod("subset", "SpectraDataFrame", .subset.SpectraDataFrame)
