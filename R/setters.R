#' Tries to guess the wavelengths from a stringr
#'
#' This simple wrapper around str_extract simply tries to pull out
#' digits contained in a chain of characters.
#' @param x a character or a vector of characters
#' @return ta character or a vector of characters containing he digits
#' that have been extracted from x.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @import plyr stringr
#' @noRd
.guessWl <- function(x){
  # simply returns the wavelengths from a string
  as.numeric(laply(.data=str_extract_all(x, "\\d"), .fun=paste, collapse=""))
}

#' Tries to find the columns of the spectra from the spectral range
#'
#' This function tries to guess the location of the columns corresponding
#' to each wavelength given by a numeric vector, based on the column
#' names. This is mainly a workaround the fact that R is usually putting
#' a "X" in front of column names that are numbers.
#'
#' @param data a \code{data.frame} object
#' @param wl a numeric or character vector containing the values of the
#' wavelengths to look for.
#' @return the index in data of the column corrresponding to the wavelengths
#' given in wl.
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @import plyr stringr
#' @noRd
.findSpectraCols <- function(data, wl, ...){

  # for each colname, we check if a part of it is in the spectral range
  # first implementation - simple test
  ind_col_spectra <- which(names(data) %in% as.character(wl))

  # if that did not succeed, we bite the bullet and
  # try to extract wl from the colnames
  if (length(ind_col_spectra) == 0) {
    # more tricky - looking for X350, etc:
    nm <- .guessWl(names(data))
    ind_col_spectra <- which(nm %in% wl)
  }
  if (length(ind_col_spectra) == 0) {
    stop('No columns found.')
  }
  ind_col_spectra
}

## setters for Spectra objects

if (!isGeneric('wl<-')) {
  setGeneric('wl<-', function(object, value)
    standardGeneric('wl<-'))
}

#' @rdname wl
setReplaceMethod("wl", "data.frame",
  function(object, value) {

    # value as to be a numeric vector
    if (is(value, 'numeric')) {
      # finding which cols contrain the spectra
      ind_nir <- .findSpectraCols(data=object, wl=value, .progress='text')
      nir <- object[, ind_nir, drop=FALSE]

      res <- Spectra(wl=as.numeric(value), nir=as.matrix(nir))

      # If there are some columns left, we use them to initiate a SpectraDataFrame object
      if (ncol(nir) < ncol(object)) {
        data <- object[, -ind_nir, drop=FALSE]
        res <- SpectraDataFrame(res, data = data)
      }
    }
    else {
      stop('Bad initialisation, please provide wavelengths as a numeric vector.')
    }
  res
  }
)

#` Replacing the wavelength range of a Spectra object. Handle with care!
#' @rdname wl
setReplaceMethod("wl", "Spectra",
  function(object, value) {

    # value as to be a numeric vector
    if (is(value, 'numeric')) {
      # if the same number of wl is provided, we simply change the content of
      # the @wl slot of the object
      if (length(value) == length(object)) {
        res <- object
        res@wl <- value
      }

      # if the number of wl provided is smaller than the number of wavelengths
      # in the @wl slot, it is interpreted as a reduction of the object to a certain
      # set of wavelengths
      else { 
        if (length(value) < length(object)) {
          ind.wl <- which(wl(object) %in% value)
          nir <- spectra(object)[, ind.wl, drop = FALSE]
          res <- Spectra(id = ids(object, as.vector = FALSE), wl = value, nir = nir, units = units(object))
          if ("data" %in% slotNames(object))
            res <- SpectraDataFrame(res, data = features(object))
          }

        else {
          stop('More wavelengths that the object originally has.')
        }
      }
    }
    else {
      stop('Please provide wavelengths as a numeric vector.')
    }

    res
  }
)


#' Parsing of the formula interface to spectra<-
#'
#' spectra(df) <- id ~ attr1 + attr2 ~ ...
#' spectra(df) <- id ~ ... ~ 350:2500
#' spectra(df) <- id ~ ...
#' spectra(df) <- ~ ...
#'
#' Inspired from Hadley Wickham's parse_formula
#' https://github.com/hadley/reshape/blob/master/R/formula.r
#'
#' @param formula a formula fot the spectra()<- setter
#' @param object a data.frame
#' @return returns a list of column names for the id slot,
#' the data slot and the nir slot of teh Spectra* object
#' @author Pierre Roudier \email{pierre.roudier@@gmail.com}
#' @import stringr plyr
#' @noRd
.parse_formula <- function(formula, object){
  formula <- str_c(deparse(formula, 500), collapse="")

  elements <- str_split(formula, fixed("~"))[[1]]
  length_elements <- aaply(elements, 1, str_length)
  elements <- elements[which(length_elements > 0)]

  formula <- lapply(str_split(elements, "[+*]"), str_trim)
  n_elements <- length(formula)

  # PLACEHOLDERS
  #
  # ... : all the columns that havent been used in the formula
  # : : sequence of integers, like 350:2500
  # :n: : sequence of numbers, emulates seq(x,y, by=n), like 350:n:2500

  # if used the ":"  placeholder
  if (any(str_detect(formula[[length(formula)]], ":"))) { # allowed only on the right hand element of the formula
    nir_seq <- aaply(unlist(str_split(formula[[length(formula)]], "[:]")), 1, as.numeric)

    # if this is a sequence of wl with by=1
    if (length(nir_seq) == 2) {
      nir_wl <- seq(nir_seq[1], nir_seq[2], by=1)
    }
    # if this is a sequence of wl with by!=1
    else {
      if (length(nir_seq) == 3) {
        nir_wl <- seq(nir_seq[1], nir_seq[3], by=nir_seq[2])
      }
      # if there's more than two ":" placeholders
      else {
        stop("Bad formula.")
      }
    }
    # finding the corresponding col names
    cols_nir  <- names(object)[.findSpectraCols(data = object, wl = nir_wl)]
    # replacing the placeholder by the actual col names
    formula[[length(formula)]] <- cols_nir
  }

  all_vars <- unlist(formula)

  # if used the "..." placeholder
  if (any(all_vars == "...")) {
    remainder <- setdiff(names(object), c(all_vars, 'id')) # setting id as a reserved name for id columns

    replace.remainder <- function(x) {
      if (any(x == "...")) {
        c(x[x != "..."], remainder)
      } else {
        x
      }
    }

    formula <- lapply(formula, replace.remainder)
  }

  if (n_elements == 1) { # case spectra(df) <- ~ ...
    cols_id <- NULL
    cols_data <- NULL
    cols_nir <- formula[[1]]
  }
  else if (n_elements == 2) {# case spectra(df) <- id ~ ...
    cols_id <- formula[[1]]
    cols_data <- NULL
    cols_nir <- formula[[2]]
  }
  else if (n_elements == 3) {# spectra(df) <- id ~ attr1 + attr2 ~ ...
    cols_id <- formula[[1]]
    cols_data <- formula[[2]]
    cols_nir <- formula[[3]]
  }
  else {
    stop('wrong formula.')
  }
  list(id = cols_id, data = cols_data, nir = cols_nir)
}

## setting the spectra of a Spectra* object
##
## - if applied to a data.frame --> we create a Spectra* object
## - if applied to a Spectra* --> we change its @nir slot
##
## IMPORTANT: The spectra() functionis for wide-formatted data. See spectra_long()
## method for long-formated data
##
if (!isGeneric('spectra<-')) {
  setGeneric('spectra<-', function(object, ..., value)
    standardGeneric('spectra<-'))
}

## for a data.frame

## samples by row
.set_spectra_by_row.data.frame <-  function(object, value) {

  # if given a formula
  if (is(value, 'formula')) {
    # parsing the formula to retrieve the different slots (id, data, nir)
    ind.vars <- lapply(.parse_formula(value, object), function(x) which(names(object) %in% x))

    if (length(ind.vars$nir) == 0) {
      ind.vars$nir <- .findSpectraCols(object, .parse_formula(value, object)$nir)
    }
    nir <- object[, ind.vars$nir, drop = FALSE]

    if (length(ind.vars$id) == 0) {
      ids <- as.character(NA)
    }
    else {
      ids <- object[, ind.vars$id, drop = FALSE]
    }

    wl <- .guessWl(names(nir))
    res <- Spectra(id = ids, wl = wl, nir = nir)

    cat("Wavelength range: ")
    cat(min(wl(res), na.rm=TRUE), " to ", max(wl(res), na.rm = TRUE)," ", wl_units(res), "\n", sep="")
    cat("Spectral resolution: ", res(wl(res)) , " ",  wl_units(res), "\n", sep="")

    if (length(ind.vars$data != 0)) {
      res <- SpectraDataFrame(res, data=object[, ind.vars$data, drop = FALSE])
    }
  }

  # if given a numeric vector (interpreted as the index of the cols)
  # eg spectra(df) <- 11:2161
  else if (is(value, 'numeric')) {
    nir <- object[, value, drop = FALSE]
    wl <- .guessWl(names(nir))
    res <- Spectra(wl=wl, nir=nir)

    # if there's some cols left, we create a SpectraDataFrame
    if (length(value) < ncol(object)) {
      data <- object[, setdiff(1:ncol(object), value), drop = FALSE]
      res <- SpectraDataFrame(res, data=data)
    }
  }

  # if given a character vector (interpreted as the names of the cols)
  # eg spectra(df) <- c('X450', 'X451', 'X452')
  else if (is(value, 'character')) {
    ind.nir <- which(names(object) %in% value)
    nir <- object[, ind.nir, drop = FALSE]

    wl <- .guessWl(names(nir))
    res <- Spectra(wl=wl, nir=nir)

    # if there's some cols left, we create a SpectraDataFrame
    if (length(value) < ncol(object)) {
      data <- object[, setdiff(1:ncol(object), ind.nir), drop = FALSE]
      res <- SpectraDataFrame(res, data=data)
    }
  }

  else {
    stop('Wrong Spectra initialisation.')
  }

  res
}

## samples by cols
.set_spectra_by_col.data.frame <- function(object, value) {
  if (is(value, 'formula')) {
    ind.vars <- lapply(.parse_formula(value, object), function(x) which(names(object) %in% x))
    
    # In that case, id = names(nir), wl = id
    nir_mat <- object[, ind.vars$nir, drop = FALSE]
    nir <- t(nir_mat)
    ids <- rownames(nir)
    wl <- as.numeric(as.character(object[, ind.vars$id]))
    # If wavelengths are not given, using 1:N.
    if (length(wl) == 0)
      wl <- 1:ncol(nir)
  }
  else {
    stop("Only the formula interface is supported by 'by_col' mode for the time being. Refer to the man page for help using it.")
  }
  
  r <- Spectra(id = ids, wl = wl, nir = nir)
  cat("Wavelength range: ")
  cat(min(wl(r), na.rm = TRUE), " to ", max(wl(r), na.rm = TRUE)," ", wl_units(r), "\n", sep = "")
  cat("Spectral resolution: ", res(wl(r)) , " ",  wl_units(r), "\n", sep = "")

  r
}

#' @rdname spectra-methods
setReplaceMethod("spectra", "data.frame", 
  function(object, ..., value) {
    
    # If no mode is being given, default to "by_row"
    dots <- list(...)
    ifelse('mode' %in% names(dots), mode <- dots$mode, mode <- "rowwise")

    if (mode == "rowwise") {
      r <- .set_spectra_by_row.data.frame(object, value)
    } 
    else {
      if (mode == "colwise") {
        r <- .set_spectra_by_col.data.frame(object, value)
      }
      else {
        stop("Wrong mode.")
      }
    }

    r
  }
)

## for a Spectra* object
#' @rdname spectra-methods
setReplaceMethod("spectra", "Spectra",
  function(object, value) {
    if (is(value, 'matrix')) {

      # this method should not allow to change the number of samples in the colection
      if (nrow(value) != nrow(spectra(object))) {
        stop("Dimensions of the matrix do not match the number of spectra in the Spectra object.")
      }
      # matrix of same dimensions is given
      if (ncol(value) == ncol(spectra(object))) {
        object@nir <- value
      }
      # matrix of different number of columns is given
      else {
        object@nir <- value
        object@wl <- as.numeric(colnames(value))
      }
    }
    else {
      stop(paste("You can't set the spectra of a Spectra* object by an object of class ", class(value), ". It has to be a matrix.", sep=""))
    }
    object
  }
)

## id

if (!isGeneric('ids<-')) {
  setGeneric('ids<-', function(object, value)
    standardGeneric('ids<-'))
}

#' @rdname ids
setReplaceMethod("ids", "Spectra",
  function(object, value) {
    if (length(value) != nrow(object)) {
      stop("length of the new ID does not match the length of the object")
    }
    if (!is.character(value)) {
      value <- as.character(value)
    }
    nm <- names(object@id)
    id <- data.frame(value)
    names(id) <- nm
    object@id <- id
    object
  }
)

#' @rdname ids
setReplaceMethod("ids", "SpectraDataFrame",
  function(object, value) {
    if (is(value, 'formula')) {
      # the id needs to be unique!
      if (length(all.vars(value)) == 1) {
        mf <- model.frame(formula=value, data=object)
        # assigning the id slot
        nm <- names(object@id)
        id <- data.frame(as.character(mf[, 1]))
        names(id) <- nm
        object@id <- id
        # removing the id col from the data slot
        object@data <- object@data[, -which(names(object@data) == names(mf))]
        # if nothing left in the data slot, back to a Spectra object!
        if (ncol(object@data) == 0) {
          object <- Spectra(wl=object@wl, nir=object@nir, id=object@id, units=object@units)
        }
      }
      else {
        stop('wrong id initialisation: id must be unique.')
      }
    }
    else{
      if (inherits(value, 'numeric')) {
        if (length(value) != length(object)) {
          stop("length of the new ID does not match the length of the object")
        }
        if (!is.character(value)) {
          value <- as.character(value)
        }
        nm <- names(object@id)
        id <- data.frame(value)
        names(id) <- nm
        object@id <- id
      }
      else {
        nm <- names(object@id)
        id <- data.frame(as.character(value))
        names(id) <- nm
        object@id <- id
      }
    }
    object
  }
)
 
