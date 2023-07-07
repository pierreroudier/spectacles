#' @include AAA-Classes.R

#' @title Constructor for the Spectra class.
#' @aliases Spectra
#' @rdname Spectra
#' 
#' @description Constructor for the Spectra class. Creates a Spectra object from scratch.
#' 
#' @param wl a numeric vector giving the wavelengths at with the spectra have
#' been measured
#' @param nir a \code{"matrix"} or a \code{"data.frame"} object giving the
#' spectra values for each sample
#' @param id a vector giving the unique id of each sample in the collection
#' @param units a character giving the unit in which the wavelengths values are
#' expressed
#' @return a new \code{"Spectra"} object
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{spectra}, \code{wl},
#' \code{Spectra-class}, \code{SpectraDataFrame}
#' @examples
#' 
#'   wls <- 350:2500
#'   id <- c("A", "B")
#'   nir <- matrix(runif(2*length(wls)), nrow = 2)
#'   s <- Spectra(wl = wls, nir = nir, id = id, units = "nm")
#' 
#' @export Spectra
"Spectra" <- function(wl=numeric(), nir=matrix(), id=as.character(NA), units="nm") {
  
  # if the wl are given as an integer vector they are translated into a numeric vector
  # for clarity (only one type to manage)
  if (is(wl, "integer"))
    wl <- as.numeric(wl)

  if (is(nir, 'data.frame'))
    nir <- as.matrix(nir)
  
  if (!is(id, "data.frame"))
    id <- data.frame(id = id)

  # If no id is given
  if (all(is.na(id))) {
    # If the object is void
    if (length(nir) == 1)
      id <- as.character(NULL)
    # if a matrix is here
    else
      id <- data.frame(id = as.character(seq(1, nrow(nir))))
  }
  # if ids are actually given by the user
  else {
    # Test of inconsistent ids when id is specified by the user

    # if theres only one spectra
    if (is.null(nrow(nir))) {
      if (nrow(id) != 1)
	stop("number of individuals and number of rows in the spectra matrix don't match")
      if ((length(wl) > 1) & (length(nir) != length(wl)))
        stop("number of columns in the spectra matrix and number of observed wavelengths don't match")
      nir <- matrix(nir, nrow=1)
    }

    # if theres more than one specta
    else {
      if (nrow(nir) != nrow(id))
        stop("number of individuals and number of rows in the spectra matrix don't match")
      if ((length(wl) > 1) & (ncol(nir) != length(wl)))
        stop("number of columns in the spectra matrix and number of observed wavelengths don't match")
      colnames(nir) <- wl
      rownames(nir) <- as.vector(do.call('rbind', id))
    }
  }

  # consistency nimber of wl/number of cols in the NIR matrix
  if ((length(wl) > 1) & (ncol(nir) != length(wl)))
    stop("number of columns in the spectra matrix and number of observed wavelengths don't match")
  
  # Making sure wavelengths are increasing
  if (!identical(wl, sort(wl))) {
    order_wl <- order(wl)
    # Re-order wavelengths
    wl <- wl[order_wl]
    # Re-order NIR matrix
    nir <- nir[, order_wl]    
  }
  
  rownames(nir) <- as.vector(do.call('rbind', id))
  colnames(nir) <- wl
  
  new("Spectra", wl = wl, nir = nir, id = id, units = units)
}

## SUMMARY

# if (!isGeneric("summary"))
#   setGeneric("summary", function(object, ...)
#     standardGeneric("summary"))

#' @title Summary 
#' @name summary
#' @description Summarize a Spectra* object.
#' @aliases summary.Spectra print.summary.Spectra
#' @usage \\method{summary}{Spectra}(object, ...)
#' @param object an object of class \code{Spectra} or \code{SpectraDataFrame}
#' @param x a result of \code{summary}
#' @param ... Additional arguments passed to \code{summary}
#' @return A \code{"summary.Spectra"} object
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' data(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' summary(oz)
#' 
#' @export summary.Spectra
#' @export
#' 
summary.Spectra <- function (object, ...){
    obj = list()
    obj[["class"]] = class(object)
    obj[["wl"]] = object@wl
    obj[["id"]] = object@id
    obj[["nir"]] = object@nir
    obj[["units"]] = object@units
    if ("data" %in% slotNames(object)) {
        if (ncol(object@data) > 1)
            obj[["data"]] = summary(object@data)
        else obj[["data"]] = summary(object@data[[1]])
    }
    else obj[["data"]] = NULL
    class(obj) = "summary.Spectra"
    obj
}

#' @return \code{NULL}
#'
#' @rdname summary
#' @method print summary.Spectra
#' @export print.summary.Spectra
#' @export
print.summary.Spectra <- function(x, ...) {
    cat(paste("Object of class ", x[["class"]], "\n", sep = ""))
    cat("Set of ", nrow(x[['id']])," spectra\n", sep = "")
    if (nrow(x[['id']]) > 0){
      cat("Wavelength range: ")
      cat(min(x[["wl"]], na.rm=TRUE), " to ", max(x[["wl"]], na.rm=TRUE)," ", x[["units"]], "\n", sep="")
      SpectralResolution <- res(x[["wl"]])
      if (length(SpectralResolution) > 1)
	cat("Spectral resolution: irregular wavelength spacing\n")
      else {
	if (length(SpectralResolution) == 0)
	  cat("Spectral resolution: NA\n")
	else
	  cat("Spectral resolution: ", SpectralResolution , " ",  x[["units"]], "\n", sep="")
      }
      if (!is.null(x$data)) {
	cat("Data attributes:\n")
	print(x$data)
      }
    }
    invisible(x)
}

## PRINT

setMethod(
  f='show',
  signature='Spectra',
  definition=function(object){
    cat(paste("Object of class ", class(object), "\n", sep = ""))
    cat("Set of ", nrow(object@id)," spectra\n", sep='')
    if (nrow(object@id) > 0){
      cat("Wavelength range: ", min(object@wl, na.rm=TRUE),"-",max(object@wl, na.rm=TRUE)," ", object@units, "\n", sep="")
      SpectralResolution <- res(object)
      if (is.null(SpectralResolution) > 1)
        cat("Spectral resolution: irregular wavelength spacing\n")
      else {
        if (length(SpectralResolution) == 0)
          cat("Spectral resolution: NA\n")
        else
          cat("Spectral resolution: ", SpectralResolution , " ", object@units, "\n", sep="")
      }
    }
    if ("data" %in% slotNames(object)) {
      cat("Data attributes:\n")
      print((object@data))
    }
  }
)

as.data.frame.Spectra <- function(x, ..., exclude_id = FALSE)  {
  df <- as.data.frame(spectra(x))
  names(df) <- wl(x)
  if (!exclude_id) {
    df <- data.frame(ids(x, as.vector = FALSE), df)
  }
  df
}

setAs("Spectra", "data.frame", function(from)
	as.data.frame.Spectra(from))

## Accessing data

# Getting the spectra matrix
# if (!isGeneric("spectra")) {
  setGeneric(
    name = "spectra", 
    def = function(object, ...) standardGeneric("spectra")
  )
# }

#' @title Retrieves or sets the spectra of a \code{Spectra*} objects.
#' @name spectra-methods
#' @rdname spectra-methods
#' @aliases spectra spectra<- spectra,data.frame-method spectra,Spectra-method spectra,SpectraDataFrame-method
#' @description Either retrieves the spectra matrix from a \code{Spectra*} object, or
#' creates a \code{Spectra*} object from a \code{"data.frame"} object different
#' interfaces detailed below.When applied to a \code{Spectra*} object, this functions simply returns the spectra it is storing.
#' @param object an object of class \code{"Spectra"} or inheriting from this class
#' @param ... see details below
#' @param value see details below
#' @details If applied on a \code{"data.frame"} object, it is an helper function to
#' create a \code{Spectra*} object. Two kind of interfaces are then available.
#' \code{value} can be: 
#' 
#' \describe{
#'   \item{a vector:}{Similarly to \code{wl}, the wavelengths of the spectra can be passed by a \code{"numeric"} vector. Alternatively, the names of the columns that contain the spectra information can be passed as a \code{"character"} vector.}
#'   \item{a formula:}{This interface is specific to inspectr. It follows a
#' scheme where differents parts can be observed, the id column, the attributes
#' columns, and the spectra columns, described by the wavelengths at which it
#' has been measured:}
#' 
#' \itemize{
#'   \item \bold{Placeholders:} 
#'     \itemize{
#'       \item \code{...} placeholder for all the columns of your data.frame object except those that have been already used in other parts of the formula. This can lead to errors. E.g. if \code{object} has data one every wavelength between 350 and 2500 nm, \code{spectra(object) <- id_field ~ ... ~ 500:2500} will stores the columns corresponding to the wavelengths 350-499 nm in its data slot!
#'       \item \code{id} For the creation of a \code{SpectraDataFrame}, it is important to always specify an id field in the formula. If no id column is present, the \code{id} placeholder will create one for you.
#'     }
#'   \item \code{spectra(object) <- ~ 350:2500} will build a \code{Spectra} object from the wavelengths between 350 and 2500, based on the column names.
#'   \item \code{spectra(object) <- ~ 350:2:2500} will build a \code{Spectra} object from the wavelengths in \code{seq(350, 2500, by = 2)}
#'   \item \code{spectra(object) <- ~ 500:2350} will build a \code{Spectra} object from the wavelengths between 500 and 2350, even though other wavelengths are present (they will be dropped)
#'   
#'   In the three later cases, the id field has been dropped (it will be automatically created). If you want to use a column of \code{"data.frame"} as an id filed, you can still use the first part of the formula:
#'   
#'   \item \code{spectra(object) <- id_field ~ 350:2500} 
#'   \item \code{spectra(object) <- id_field ~ 350:5:2500}
#'   
#'   Some data can also be added to the object, which will then be of \code{SpectraDataFrame} class:
#'   
#'   \item \code{spectra(object) <- id_field ~ property1 ~ 500:2300} will create a \code{SpectraDataFrame} with ids from the id_field column, data from the property1 column, and spectral information between 500 and 2300 nm. That means that data property2, and all spectral information from bands < 500 nm or > 2300 nm will be dropped
#'   
#' You can also combine the placeholders:
#' 
#'   \item \code{spectra(object) <- id_field ~ ... ~ 350:2500} will create a \code{SpectraDataFrame} object with ids from the id_field column, all spectral bands between 350 and 2500 nm. The data slot is given all the remaining columns.  
#'  }
#' }
#' 
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(oz)
#' class(oz) # this is a simple data.frame
#' # structure of the data.frame: it is rowwise-formatted
#' big.head(oz) 
#' 
#' ## CREATING Spectra OBJECTS
#' ##
#' 
#' # Using spectra() to initiate a Spectra from 
#' # the data.frame
#' spectra(oz) <- sr_no ~ 350:2500
#' 
#' # It is possible to select wavelengths using the formula interface
#' data(oz)
#' spectra(oz) <- sr_no ~ 350:5:2500
#' 
#' data(oz)
#' spectra(oz) <- sr_no ~ 500:1800
#' 
#' ## CREATING SpectraDataFrame OBJECTS
#' ##
#' 
#' # Using spectra() to initiate a SpectraDataFrame from 
#' # the data.frame
#' data(oz)
#' spectra(oz) <- sr_no ~ carbon + ph + clay ~ 350:2500
#' 
#' # Selecting data to be included in the SpectradataFrame object
#' data(oz)
#' spectra(oz) <- sr_no ~ carbon ~ 350:2500
#' 
#' # Forcing the creation of new ids using the id keyword in the 
#' # formula interface
#' data(oz)
#' spectra(oz) <- id ~ carbon ~ 350:2500
#' ids(oz, as.vector = TRUE)
#' 
#' # Using the "..." short-hand to select all the remaining columns
#' data(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' 
#' ## CREATING Spectra OBJECTS FROM
#' ## BY-COLS-FORMATTED DATA
#' ##
#' 
#' # For data formatted in the colwise format, 
#' # use the "colwise" mode
#' 
#' # Transforming data into colwise format
#' # for demonstration's sake
#' #
#' m <- melt_spectra(oz)
#' oz_by_col <- reshape2::acast(m, ... ~ sr_no)
#' oz_by_col <- data.frame(
#'   wl = rownames(oz_by_col), 
#'   oz_by_col, 
#'   check.names = FALSE)
#' 
#' # Here's colwise-formatted data 
#' big.head(oz_by_col)
#' 
#' # Convert it into Spectra object
#' spectra(oz_by_col, mode = "colwise") <- wl ~ ...
#' 
#' # Then data can be added to promote it as a SpectraDataFrame
#' my.data <- features(oz)
#' features(oz_by_col) <- my.data
#' 
setMethod("spectra", "Spectra",
  function(object) {
    res <- object@nir
    colnames(res) <- object@wl
    res
  }
)

# Getting the wavelengths
# if (!isGeneric("wl"))
  setGeneric("wl", function(object, ...)
    standardGeneric("wl"))

#' @title Retrieves or sets the wavelengths of a \code{Spectra*} object.
#' 
#' @description Either retrieves the wavelengths from a \code{Spectra*} object, or creates a
#' \code{Spectra*} object from a \code{"data.frame"} object by setting some of
#' its columns as the wavelengths.
#' 
#' When applied to a \code{Spectra*} object, this functions simply returns the
#' wavelengths of the spectra it is storing.
#' 
#' If applied on a \code{"data.frame"} object, it is an helper function to
#' create a \code{Spectra*} object. It then needs to be indicated the
#' wavelengths at which the spectra values are measured. The assumption is that
#' each row of the \code{"data.frame"} is a spectra, and the column names of
#' the \code{"data.frame"} contain the wavelengths values.
#' 
#' If all the columns are used to create the \code{Spectra*} object, a
#' \code{Spectra} object is created. If some attributes are left, they will be
#' used to generate a \code{SpectraDataFrame} object.
#' 
#' @name wl
#' @aliases wl wl<- wl,Spectra-method wl<-,Spectra-method
#' wl<-,data.frame-method
#' @docType methods
#' @param object a \code{"data.frame"} or an object inheriting from class
#' \code{Spectra}
#' @param value the wavelengths of the \code{Spectra*} object to create
#' @return If applied on a \code{"data.frame"}, either a \code{Spectra} or a
#' \code{SpectraDataFrame} object. If applied on a \code{Spectra*} object, a
#' vector.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @seealso \code{spectra}, \code{Spectra-class},
#' \code{SpectraDataFrame-class}
#' @examples
#' 
#' # Loading example data
#' data(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' 
#' # Retrieving wavelengths from Spectra* object
#' wl(oz)
#' 
#' # Replacing wavelength values - USE WITH CAUTION!
#' wl(oz) <- 1:length(oz)
#' wl(oz)
#' 
#' # Use to initiate a Spectra* object from a data.frame
#' data(oz)
#' wl(oz) <- 350:2500
#' ids(oz) <- ~ sr_no
#' 
#' @export wl
setMethod("wl", "Spectra",
  function(object)
    object@wl
)

# Getting the ids
# if (!isGeneric("ids"))
  setGeneric("ids", function(object, ...)
    standardGeneric("ids"))

#' @title Retrieves or sets the ids of a \code{Spectra*} object.
#' 
#' @description Either retrieves the wavelengths from a \code{Spectra*} object, or creates a
#' \code{Spectra*} object from a \code{"data.frame"} object by setting some of
#' its columns as the wavelengths. The \code{"ids<-"} method for 
#' \code{SpectraDataFrame} objects allows to use a formula interface to use a 
#' column in its \code{data} slot as the object IDs (see the last example provided 
#' in the Examples section).
#' 
#' @name ids
#' @rdname ids
#' @aliases ids ids<- ids,Spectra-method ids<-,Spectra-method
#' ids<-,SpectraDataFrame-method
#' @docType methods
#' @usage 
#' ids(object, ...)
#' ids(object) <- value
#' @param object an object of class \code{"Spectra"} or inheriting from this class
#' @param ... \describe{
#'   \item{\code{as.vector}}{Controls whether the IDs are returned as a vector or as a data.frame (defaults to TRUE)}
#' }
#' @param value character vector for new IDs
#' @return The \code{ids} methods return a vector if \code{as.vector} is TRUE,
#' a \code{data.frame} otherwise. The \code{"ids<-"} method return a
#' \code{SpectraDataFrame} object (or a \code{Spectra} object if the column in
#' the data slot that has been used to initiate the IDs was the only
#' attribute).
#' @section Methods: \itemize{
#'  \item \code{ids(object, ..., as.vector = TRUE)} 
#'  \item \code{ids(object) <- value} 
#' }
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' 
#' # Retrieving ids
#' ids(oz)
#' 
#' # Setting ids using a vector of values
#' ids(oz) <- seq_len(nrow(oz))
#' ids(oz)
#' 
#' # Setting ids using an attribute
#' oz$new_id <- seq_len(nrow(oz)) + 1000
#' ids(oz) <- ~ new_id
#' ids(oz)
#' 
#' @export ids
setMethod("ids", "Spectra",
  function(object, ..., as.vector = TRUE) {
    if (as.vector) {
      as.character(res <- object@id[[1]])
    } else {
      res <- object@id
    }
    res
  }
)

# Getting the units
# if (!isGeneric("wl_units"))
  setGeneric("wl_units", function(object)
    standardGeneric("wl_units"))

#' @title Wavelengths of Spectra* objects
#' @name wl_units
#' @description Retrieves the wavelengths units from \code{Spectra*} object
#' @aliases wl_units wl_units<- wl_units,Spectra-method wl_units<-,Spectra-method 
#' @docType methods
#' @usage wl_units(object)
#' @param object an object inheriting from class \code{Spectra}
#' @param value a character string
#' @return A vector
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' # Loading example data
#' data(oz)
#' spectra(oz) <- sr_no ~ ... ~ 350:2500
#' 
#' # Print wavelength information
#' wl(oz)
#' range(wl(oz))
#' 
#' # Manipulate wavelength information
#' wl(oz) <- wl(oz) + 1000
#' wl(oz)
#' 
setMethod("wl_units", signature = "Spectra",
  function(object)
    object@units
)

# if (!isGeneric('wl_units<-'))
  setGeneric('wl_units<-', function(object, value)
    standardGeneric('wl_units<-'))

#' @rdname wl_units
setReplaceMethod("wl_units", "Spectra",
  function(object, value) {
    if (!is.character(value) | length(value) != 1)
      stop("Units have to be passed as a single character string.")
    object@units <- value
    object
  }
)

#' @title Retrieve dimensions of Spectra* objects
#' @name dimensions
#' @aliases nrow,Spectra-method ncol,Spectra-method dim length dim.Spectra
#' length.Spectra
#' @docType methods
#' @description Retrieves the wavelengths units and the spectral resolution from
#' \code{Spectra*} objects.
#' 
#' @param x For \code{nrow}, \code{length}, \code{dim}, a \code{Spectra}
#' object. For \code{ncol}, a \code{SpectraDataFrame} object.
#' @return \code{nrow}, \code{ncol}, \code{nwl}, and \code{length}, return an
#' \code{integer}, \code{dim} returns a vector of length 2 or 3 (see Details section).
#' 
#' @details 
#' * Methods for \code{Spectra} objects
#' 
#' \code{nrow} returns the number of individuals in the collection
#' \code{length} returns the number of wavelengths in the collection
#' \code{ncol} returns NULL \code{dim} returns a vector containing (1) the
#' number of individuals and (2) in the number of wavelengths in the collection
#' 
#' * Methods for \code{Spectra} objects
#' 
#' \code{nrow} returns the number of individuals in the collection
#' \code{length} returns the number of wavelengths in the collection
#' \code{ncol} returns the number of attributes in the collection \code{dim}
#' returns a vector containing (1) the number of individuals, (2) in the number
#' of wavelengths, and (3) the number of attributes in the collection
#' 
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' nrow(australia)
#' ncol(australia)
#' length(australia)
#' dim(australia)
#' 
length.Spectra <- function(x)
    ncol(x@nir)

#' Returns the number of samples in the object
#' @rdname dimensions
setMethod(f='nrow', signature='Spectra',
definition = function(x)
  nrow(ids(x, as.vector = FALSE))
)

#' @rdname dimensions
setMethod(f='ncol', signature='Spectra',
definition = function(x) {
  if ("data" %in% slotNames(x)) {
    n <- ncol(x@data)
  } else {
    n <- NULL
  }
  n
}
)

#' @rdname dimensions
dim.Spectra <- function(x) {
  r <- c(nrow(x), length(x))
  if ('data' %in% slotNames(x)) {
    r <- c(r, ncol(features(x)))
  }
  r
}

## Returns spectral resolution of the wavelengths

# if (!isGeneric("res"))
  setGeneric("res", function(x)
    standardGeneric("res"))

.res.numeric <- function(x){
  unique(round(diff(x), digits = 10)) # round - otherwise diff() picks some unsignificant values
}

#' @rdname res
setMethod("res", "numeric", .res.numeric)
#' @rdname res
setMethod("res", "integer", .res.numeric)

.res.Spectra <- function(x) {
  r <- unique( round( diff(wl(x)), digits = 10) )
  if (length(r) > 1) r <- NULL # if resolution is not regular
  r
}

#' @title Spectral resolution
#' @description Returns the spectral resolution of an object
#' @name res
#' @aliases res,Spectra-method res,SpectraDataFrame-method 
#' @param x object an object inheriting from \code{Spectra}
#' @return a vector
#'
#' @method resolution Spectra
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
setMethod("res", "Spectra", .res.Spectra)

## overloads

#' @title Extracting and replacing parts of Spectra* objects
#' @name extraction-methods
#' @description These methods emulates classic base methods '[', '[[' and '$' to extract or replace parts of Spectra* objects.
#' 
#' @aliases [ [<- [[ [[<- $ $<- [,Spectra-method [[,Spectra-method [[<-,Spectra-method [,Spectra,ANY,ANY,missing-method [,SpectraDataFrame,ANY,ANY,missing-method [[,SpectraDataFrame,ANY,missing-method [[<-,Spectra,ANY,missing-method [<-,SpectraDataFrame-method $<-,Spectra-method $,SpectraDataFrame-method  [[<-,Spectra,ANY,missing,ANY-method 
#' 
#' 
#' @docType methods
#' 
#' @usage 
#' \\S4method{[}{Spectra}(x,i,j,\dots,drop=FALSE)
#' \\S4method{[[}{Spectra}(x,i,j,\dots)
#' \\S4method{$}{SpectraDataFrame}(x,name)
#' \\S4method{$}{Spectra}(x,name) <- value
#' \\S4method{[[}{Spectra}(x,i,j,\dots) <- value
#' 
#' @param x an object of class \code{Spectra} or \code{SpectraDataFrame}
#' @param i,j,... indices specifying elements to extract or replace
#' @param drop currently ignored
#' @param name A literal character string or a name (possibly backtick quoted)
#' @param value typically an array-like R object of a similar class as x
#' 
#' @return These methods either return an object of the same class as \code{x},
#' or can promote a \code{Spectra} object to a \code{SpectraDataFrame} object
#' by adding data ("[[<-" and "$<-" methods).
#' @section Methods: \describe{
#' 
#' \bold{x=Spectra}
#' 
#' \code{x[i, j, ..., drop = FALSE]}
#' 
#' \code{x$name <- value}
#' 
#' \code{x[[name]] <- value}
#' 
#' \tabular{rll}{ \tab \code{x} \tab A \code{Spectra} object \cr \tab \code{i}
#' \tab Row index of the selected individuals \cr \tab \code{j} \tab Selected
#' wavelengths \cr \tab \code{name} \tab A literal character string or a name
#' \cr \tab \code{...} \tab Ignored \cr \tab \code{drop} \tab Ignored \cr }
#' 
#' \bold{x=SpectraDataFrame}
#' 
#' \code{x[i, j, k, ..., drop = FALSE]}
#' 
#' \code{x[[name]]}
#' 
#' \code{x[[name]] <- value}
#' 
#' \code{x$name}
#' 
#' \code{x$name <- value}
#' 
#' \tabular{rll}{ \tab \code{x} \tab A \code{SpectraDataFrame} object \cr \tab
#' \code{i} \tab Row index of the selected individuals \cr \tab \code{j} \tab
#' Selected wavelengths \cr \tab \code{k} \tab Selected columns in the @@data
#' slot\cr \tab \code{name} \tab A literal character string or a name \cr \tab
#' \code{...} \tab Ignored \cr \tab \code{drop} \tab Ignored \cr }}
#' 
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Getting features information from SpectraDataFrame
#' australia$carbon
#' australia[['carbon']]
#' 
#' # Creating new features
#' australia$foo <- runif(nrow(australia))
#' australia[['bar']] <- runif(nrow(australia))
#' 
#' # Replacing values
#' australia$foo <- sample(
#'   LETTERS[1:5], 
#'   size = nrow(australia), 
#'   replace = TRUE
#' )
#' australia[['bar']] <- sample(
#'   c(TRUE, FALSE), 
#'   size = nrow(australia), 
#'   replace = TRUE
#' )
#' 
#' # Promote Spectra to SpectraDataFrame
#' s <- as(australia, 'Spectra')
#' class(s)
#' s$foo <- runif(nrow(s))
#' 
setMethod("[", c("Spectra", "ANY", "ANY", "missing"),
  function(x, i, j, ..., drop = FALSE) {

    missing.i <- missing(i)
    missing.j <- missing(j)

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

    # WAVELENGTHS
    if (missing.j) {
      j <- TRUE
    } else {
      # If the indices are all negative, cols are removed
      if (all(j < 0)) {
        j <- setdiff(as.numeric(wl(x)), abs(j))
      }
      j <- which(as.numeric(wl(x)) %in% j)

      # If some of teh WL that have been asked for are not in
      # the wl range, these are removed from j and a warning is produced
      if (any(!(j %in% wl(x)))) {
        j.in.wl <- which(j %in% wl(x))
        j <- j[j.in.wl]
      }
      
      if (length(j) == 0)
        stop("Wrong wavelengths selection. Your wavelength selection is not available.")
    }

    Spectra(wl = x@wl[j], nir=x@nir[i, j, drop = FALSE], id = x@id[i, , drop = FALSE]) 
  }
)


## Promote a Spectra object to a SpectraDataFrame

# if (!isGeneric('features<-'))
  setGeneric('features<-', function(object, ..., value)
    standardGeneric('features<-')
  )

setReplaceMethod("features", signature("Spectra", "ANY"),
  # safe enables id check
  # key gives the column name of the ids in the data.frame
  function(object, ..., value) {
     
    # hack to avoid the 'value' must be on the right side' thing at R CMD check
    dots <- list(...)
    ifelse('safe' %in% names(dots), safe <- dots$safe, safe <- TRUE)
    key <- dots$key # NULL if key is not in dots
    ifelse('exclude_id' %in% names(dots), exclude_id <- dots$exclude_id, exclude_id <- TRUE)

    if (!inherits(value, "data.frame"))
      stop('invalid initialization for SpectraDataFrame object')

    if (safe) {
      if (is.null(key))
        stop("In the safe mode, you need to provide either the column name of the sample ids to the key option.")
      if (length(key) != 1)
        stop("Please provide only ONE id column.")

      # Actual ID sanity check
      spectra_ids <- ids(object, as.vector = FALSE)
      if (is.numeric(key)) {
        key <- names(value)[key]
      }

      # Using the "key" name for ids
      names(spectra_ids) <- key
      
      # Put data together
      data <- join(spectra_ids, value,  by = key, type = "left", match = "first")
      
      # removing the id column      
      if (exclude_id)
        data <- data[, -1*which(names(data) == key)]
    }
    else {
      warning("Sample ID check has been disabled. This mode assumes you made sure the order of the rows in your data is consistent with the order in which these samples appear in the Spectra object.")
      data <- value
    }
    SpectraDataFrame(object, data = data)
  }
)

## Adding objects together
# Maybe to be moved into the Spectra() and SpectraDataFrame() method.

# if (!isGeneric("add"))
#   setGeneric("add", function(x, y, ...)
#     standardGeneric("add"))
# 
# #' Adds two Spectra objects together
# #'
# #' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
# .add.Spectra <- function(x, y){
#   tmp <- list()
# 
#   if (identical(x@wl, y@wl))
#     tmp$wl <- x@wl
#   else
#     stop('You can not add objects with different wavelength ranges')
# 
#   if (identical(ncol(x@wl), ncol(y@wl)))
#     tmp$nir <- rbind(x@nir, y@nir)
#   else
#     stop('You can not add objects with different wavelength ranges')
# 
#   if (!any(x@id %in% y@id))
#     tmp$id <- rbind(x@id, y@id)
#   else
#     stop('You can not add objects with overlapping IDs')
# 
#   if (x@units %in% y@units)
#     tmp$units <- x@units
#   else
#     stop('You can not add objects with different wavelength units')
# 
#   if (("data" %in% slotNames(x)) & ("data" %in% slotNames(y))) {
#     tmp$data <- join(x@data, y@data, type="full")
#     res <- SpectraDataFrame(wl=tmp$wl, nir=tmp$nir, id=tmp$id, units=tmp$units, data=tmp$data)
#   }
#   else
#     res <- Spectra(wl=tmp$wl, nir=tmp$nir, id=tmp$id, units=tmp$units)
# 
#   res
# }
# 
# add.Spectra <- function(...){
#   dotargs <- list(...)
#   if ( !all(sapply(dotargs, function(x) is(x,"Spectra") )) )
#     stop('the arguments must be Spectra objects')
# 
#   res <- dotargs[[1]]
#   if (nargs() >= 2) {
#     for (i in 2:length(dotargs))
#       res <- .add.Spectra(res, dotargs[[i]])
#   }
#   res
# }
# 
# setMethod("add", signature=c("Spectra", "Spectra"),
#   function(x,y,...) add.Spectra(x, y, ...))
# 
# setMethod("add", signature=c("SpectraDataFrame", "SpectraDataFrame"),
#   function(x,y,...) add.Spectra(x, y, ...))

#' @title Stacking \code{Spectra} objects together
#' @name rbind 
#' @description This method stacks two or more \code{Spectra*} objects together.
#' @aliases rbind.Spectra rbind.SpectraDataFrame
#' @param \dots The \code{Spectra} objects to be combined.
#' @param create_new_ids allows creation of new ids if the ids of the
#' \code{Spectra*} objects you are trying to stack are redundant
#' @param new_ids vector of new ids to be given to the new object
#' @return a \code{Spectra*} object.
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' s <- rbind(australia, australia, create_new_ids = TRUE)
#' 
#' l <- separate(australia, calibration = 0.6)
#' s <- rbind(l$validation, l$calibration)
rbind.Spectra <- function(..., create_new_ids = FALSE, new_ids = NULL) {
  dots <- list(...)
  names(dots) <- NULL

  # wl
  wls <- lapply(dots, wl)
  wl_test <- all(laply(wls, function(x) identical(wls[[1]], x)))
  if (!wl_test) {
    stop("To be added together, all Spectra objects must share the same wavelengths.")
  } else {
    wls <- apply(do.call('rbind', wls), 2, unique)
  }

  # units
  wl_uts <- unique(laply(dots, wl_units))
  if (length(wl_uts) > 1)
    stop("To be added together, all Spectra objects must share the same wavelength units.")

  # nir
  nir <- do.call("rbind", lapply(dots, spectra))
  
  # id
  if (is.null(new_ids)) { # We try to keep the old IDs
    ids <- do.call('c', lapply(dots, ids, as.vector = TRUE))
  } else {
    # New IDs are provided
    if (length(new_ids) == length(ids))
      ids <- new_ids
    else
      stop("new_ids must have the same length as the total number of spectra you are trying to collate together.")
  }
    
  # If ids are not unique
  if (length(unique(ids)) != length(ids)) {
    if (create_new_ids) {
      warning("Redundant IDs found. Creating new set of unique IDs.")
      ids <- 1:length(ids)
    }
    else {
      stop("Redundant IDs found. Either allow the creation of new ids using the create_new_ids option, or provide the function with a set of unique ids using the new_ids option.")
    }
  }
  
  res <- Spectra(wl = wls, nir = nir, id = ids, units = wl_uts)

  # data
  test_data <- laply(dots, function(x) "data" %in% slotNames(x))
  if (all(test_data)) {
    # Unify id colname for join
    data <- llply(dots, features, exclude_id = FALSE)
    data <- llply(data, function(x) {names(x)[1] <- 'id'; x})
    data <- do.call("rbind", data)
    res <- SpectraDataFrame(res, data = data)
  }
  
  res
}

#' @name rbind
rbind.SpectraDataFrame <- rbind.Spectra

#' @title Split a Spectra* object using factors
#' @name split
#' @aliases split split.Spectra split,Spectra-method
#' @docType methods
#' @description Splits a Spectra* object into groups using a factor, either a provided as a vector or as an attribute in the features of the object.
#' 
#' @details This is an adaptation of the \code{split} function in the base package.
#' @param x Spectra object
#' @param f either a vector of factors (for objects inheriting from
#' \code{Spectra}), or the name or indice of an attribute in the data slot of
#' an object inheriting from \code{SpectraDataFrame}
#' @param drop ignored
#' @param ... further potential arguments passed to methods.
#' @return A list of objects of same class as \code{x}.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # On a Spectra object, we need to provide a vector of factors
#' # to split the object
#' s <- as(australia, 'Spectra')
#' # We make up some kind of factor to split the data. 
#' idx <- sample(letters[1:5], replace = TRUE, size = nrow(s)) # This is a vector
#' r <- split(s, idx)
#' str(r)
#' 
#' # On a SpectradataFrame object, we can also provide the name or index 
#' # of an attribute
#' 
#' # Generate some kind of factor
#' australia$fact <- sample(LETTERS[1:3], size = nrow(australia), replace = TRUE) 
#' r <- split(australia, 'fact')
#' str(r)
#' 
#' # A list is returned, and is thus ready for use with lapply, or any
#' # of the l*ply functions from the plyr package
#' lapply(r, nrow)
setMethod("split", "Spectra", function(x, f, drop = FALSE, ...){
  
  # If length(f) <= 1, we consider f is giving the colname or index
  if (length(f) <= 1) {
    f <- features(x)[, f]
  }

  lapply(split(seq_len(nrow(x)), f, ...), function(ind) x[ind, ])
})

#' @title Mutate a Spectra* object by transforming the spectra values, and/or adding
#' new or replacing existing attributes.
#' 
#' @description This function is very similar to \code{transform} but it executes the
#' transformations iteratively so that later transformations can use the
#' columns created by earlier transformations. Like transform, unnamed
#' components are silently dropped.
#' 
#' Either the spectra, and/or the attributes (if the \code{.data} inherits from
#' the \code{SpectraDataFrame} class) can be affected: \itemize{ \item To
#' affect the spectra, one should use the \code{nir} placeholder, eg \code{nir
#' = log(1/nir)} \item To affect the attributes of the object, the definitions
#' of new columns are simply given using attributes names, \code{newAttribute =
#' 1/sqrt(attribute)} \item Both spectra and attrbutes can be transformed in
#' one command.}
#' 
#' @name mutate
#' @aliases mutate mutate.Spectra mutate,Spectra-method
#' @docType methods
#' @param .data an object inheriting from the \code{Spectra} class
#' @param ... named parameters giving definitions of new columns
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}, from code from
#' Hadley Wickham
#' @references Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data
#' Analysis. Journal of Statistical Software, 40(1), 1-29. URL
#' http://www.jstatsoft.org/v40/i01/.
#' @export
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Modifying spectra
#' m <- mutate(australia, nir = log1p(1/nir))
#' plot(m)
#' 
#' # Modifying and creating attributes
#' m <- mutate(
#'   australia, 
#'   sqrt_carbon = sqrt(carbon), 
#'   foo = clay + ph, 
#'   nir = log1p(1/nir)
#' )
#' plot(m)
setMethod("mutate", "Spectra", function (.data, ...){

  wls <- wl(.data)
  uns <- wl_units(.data)
  ids <- ids(.data, as.vector = FALSE)

  cols <- as.list(substitute(list(...))[-1])
  cols <- cols[names(cols) != ""]

  # transformations on the spectra
  if (any(names(cols) == 'nir')) {
    nir <- melt(spectra(.data), varnames=c('id', 'wl'), value.name = "nir")
    nir[["nir"]] <- eval(cols[["nir"]], nir, parent.frame())
    nir <- acast(nir, id ~ wl)
    # remove it from the cols list
    cols[['nir']] <- NULL
  }
  # no transformations on the spectra
  else
    nir <- spectra(.data)

  r <- Spectra(wl = wls, nir = nir, id = ids, units = uns)

  # transformations on the data - only for classes inheriting from SpectraDataFrame
  if (("data" %in% slotNames(.data)) & (length(cols) > 0)) { # testing if theres transformations left
    d <- sapply(cols, function(x) eval(x, features(.data), parent.frame()))
    r <- SpectraDataFrame(r, data = data.frame(features(.data), d))
  }

  r
})

## Separate calibration set vs validation set

# if (!isGeneric("separate"))
  setGeneric("separate", function(obj, calibration, ...)
    standardGeneric("separate"))

#' @title Separates a \code{Spectra*} object into a calibration and a validation set.
#' @name separate
#' @aliases separate separate.Spectra separate,Spectra-method
#' @docType methods
#' @description Separates a \code{Spectra*} object into a calibration and a validation set.
#' @param obj an object inheriting from class \code{SpectraDataFrame}
#' @param calibration The fraction of the dataset to be put in the calibration
#' set
#' @return An list with two \code{SpectraDataFrame} objects, one for the
#' calibration, and the other for the validation.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' l <- separate(australia, calibration=0.7)
#' # The result is a list of two Spectra* objects
#' str(l)
#' lapply(l, nrow)
setMethod("separate", "Spectra", function(obj, calibration){
  if (calibration < 1)
    calibration <- floor(calibration*nrow(obj))
  calib <- sample(x=seq_len(nrow(obj)), size=calibration, replace = FALSE)
  valid <- setdiff(seq_len(nrow(obj)), calib)
  list(calibration=obj[calib, ], validation=obj[valid, ])
})

# if (!isGeneric('melt_spectra'))
  setGeneric('melt_spectra', function(obj, attr = NULL,...)
    standardGeneric('melt_spectra')
)

#' @name melt_spectra
#' @aliases melt_spectra melt_spectra-methods melt_spectra,Spectra-method
#' melt_spectra,SpectraDataFrame-method
#' @docType methods
#' @title Melts the spectra data of a Spectra object and returns it as wide format.
#' 
#' @description This function is very useful when wanting to plot spectra using the lattice or ggplot2 packages
#' @usage 
#' melt_spectra(obj, attr=NULL, ...)
#' @param obj an object of class \code{"Spectra"} or inheriting from this class
#' @param attr vector of id variables against which the spectra will be melted (see \code{melt})
#' @param ... further arguments passed to or from other methods
#' @section Methods: \describe{
#' 
#' \bold{x=Spectra}
#' 
#' \code{melt_spectra(obj, ...)}
#' 
#' \tabular{rll}{ \tab \code{obj} \tab A \code{Spectra} object \cr \tab
#' \code{...} \tab Ignored \cr }
#' 
#' \bold{x=SpectraDataFrame}
#' 
#' \code{melt_spectra(obj, attr=NULL, ...)}
#' 
#' \tabular{rll}{ \tab \code{obj} \tab A \code{SpectraDataFrame} object \cr
#' \tab \code{attr} \tab Character, the name of an attribute in the object data
#' to split the spectra against. \cr \tab \code{...} \tab Ignored \cr }
#' 
#' }
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @references Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data
#' Analysis. Journal of Statistical Software, 40(1), 1-29. URL
#' http://www.jstatsoft.org/v40/i01/.
#' @import reshape2
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Simple melt
#' r <- melt_spectra(australia)
#' head(r)
#' 
#' \dontrun{
#' # Melt against some factor (or continuous data), and plot
#' # using ggplot2
#' 
#' # Create some factor
#' australia$fact <- sample(
#'   LETTERS[1:3], 
#'   size = nrow(australia), 
#'   replace = TRUE
#' )
#' r <- melt_spectra(australia, attr = 'fact')
#' 
#' # Create plot
#' library(ggplot2)
#' p <- ggplot(r) + 
#'   geom_line(aes(x=wl, y=nir, group=id, colour=fact)) + 
#'   theme_bw()
#' print(p)
#' }
#' 
setMethod("melt_spectra", "Spectra", function(obj, attr = NULL, ...){

  id.nm <- names(ids(obj, as.vector = FALSE))
  x <- data.frame(ids(obj, as.vector = FALSE), spectra(obj))
  names(x) <- c(id.nm, wl(obj))
  res <- melt(x, id.vars = id.nm, variable.name = 'wl', value.name="nir") 
  res$wl <- as.numeric(as.character(res$wl))
  res
})

## Selecting/cutting wavelengths

# Negative values will be to remove, positive to select
# if (!isGeneric("cut")) 
    setGeneric("cut", function(x, ...)
        standardGeneric("cut"))

#' Manipulating the wavelength range of a \code{Spectra} object
#' 
#' This methods allows to either select or remove a specific range of
#' wavelengths from a \code{Spectra} object.
#' 
#' The wavelengths are extracted if \code{wl > 0}, or removed if \code{wl < 0}.
#' You can't mix negative and positive index in \code{wl}.
#' 
#' @name cut
#' @aliases cut cut,Spectra-method
#' @docType methods
#' @param x an object inheriting from class \code{Spectra}
#' @param wl a vector of the wavelengths to either select or remove from
#' \code{x}
#' @param ... ignored
#' @return An object of the same class as \code{x}.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @examples
#' 
#' # Loading example data
#' data(australia)
#' spectra(australia) <- sr_no ~ ... ~ 350:2500
#' 
#' # Extracting a specific wavelengths range
#' s <- cut(australia, wl = 450:550)
#' plot(s)
#' 
#' s <- cut(australia, wl = c(450:550, 1850:2050))
#' plot(s)
#' 
#' # Removing a specific wavelengths range
#' s <- cut(australia, wl = -1*450:550)
#' plot(s)
#' 
#' s <- cut(australia, wl = -1*c(450:550, 1850:2050))
#' plot(s)
#' 
setMethod("cut", "Spectra", function(x, ..., wl) {
  
  # If wl is negative, we REMOVE these
  if (any(wl < 0) & any(wl > 0) | any(wl == 0))
    stop("You can't mix positive and negative wavelengths, or use zero.")

  if (all(wl < 0)) {
    wl <- abs(wl)
    wl <- setdiff(wl(x), wl)
  } 

  # Checking that wl in available wavelengths
  if (!all(wl %in% wl(x))) {
    stop("Selected wavelengths not present in the object")
  }

  # Getting indices of the wavelengths to select
  idx <- laply(wl, function(w) which(wl(x) == w))
  # Subsetting spectra matrix
  nir <- spectra(x)[, idx]
  
  res <- Spectra(wl = wl, nir = nir, id = ids(x, as.vector = FALSE), units = wl_units(x))

  if ("data" %in% slotNames(x)) {
    res <- SpectraDataFrame(res, data = features(x))
  }
  
  res
})
