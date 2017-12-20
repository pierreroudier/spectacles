[![Travis-CI Build Status](https://travis-ci.org/spectacles/spectacles.svg?branch=master)](https://travis-ci.org/pierreroudier/spectacles)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spectacles)](https://CRAN.R-project.org/package=spectacles)

# Welcome to the spectacles project page

## Scope of the Package

The `spectacles` package is making it easy (or at least *easier*!) to handle spectroscopy data. It provides the user with dedicated classes (namely `Spectra` and `SpectraDataFrame`), so that most of the useful information about the spectral dataset is available in one R object:

* the spectral values 
* the wavelengths at which these have been recorded
* some kind of ID
* if available, some associated data (typically, some lab measurements)

## Installation

`spectacles` is not on CRAN (yet?), but you can install the latest version using the `devtools` package:

```
# Install devtools if you don't have it on your machine
# install.packages('devtools')
devtools::install_github("pierreroudier/spectacles")
```

## Graphical Capabilities

It also provides easy ways to plot a collection of spectra:

* simple line plots of the individual spectra
* offset plots of the individual spectra
* stacked plots of the individual spectra
* summary plots of a whole collection, or aggregated against a given factor
* tools to code more advanced visualisations yourself using eg `ggplot2` or `lattice`

It also gives overloads to the most common operators such as `$`, `[`, or `[[`, so that any user familiar with `data.frame` object would fell right at home.

## Processing

The philosophy of the package is really just to make it easier to work with quite complex data. There are a lot of tools already existing in R to do spectral preprocessing (`signal`, etc.). A few additional tools have been added in `spectacles`, such as the ASD splice correction. 

The idea is for the package to work quite well with the pipe (`%>%`) operator from the `magrittr` package, to create chains of pre-processing operators. The function `apply_spectra` makes it easy to work with any function whose input is either a `numeric` vector or a `matrix`:

```
# Example of splice correction, followed by
# a first derivative, followed by a SNV

my_spectra %>% 
  splice %>% 
  apply_spectra(diff, 1) %>%
  apply_spectra(snv)
  
# Another example using prospectr
my_spectra %>% 
  splice %>% 
  apply_spectra(prospectr::continuumRemoval, wav = wl(.)) %>% 
  plot
```

## Regression and Classification

Again, lots of existing methods available, so `spectacles` is not re-implementing any of these. There's various ways to use `spectacles` with the different methods available, but my favoured option is to use it in conjonction with the `caret` package, which gives a unique API to 160+ models in R:

```
fit <- train(
  y = s$carbon,
  x = spectra(s),
  method = "pls"
)
spectroSummary(fit)
```

## Hey, that sounds a lot like `inspectr`?!?

Yes, I once had a package called `inspectr` on Github, and `spectacles` is very much the continuation of `inspectr`. The only reason why `inspectr` changed name is that someone pushed a package called `inspectr` on CRAN (despite `inspectr` being quite visible on Github.... :-/). So, lesson learnt this time!
