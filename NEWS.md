# rasterdiv 0.3.6
## Minor changes
* Refixes \alias to `load_copNDVI()` (see warining raised by Brian Ripley)
* Sets Classical Rao's = 0 for window with only 1 non-NA value (see bug raised by Jakub Nowosad)

# rasterdiv 0.3.5p
## Minor changes
* Adds \alias to load_copNDVI()
* Removes TODO document

# rasterdiv 0.3.5
## Major changes
* Adds `twdtw` as a distance methods for `paRao()`
* Adds back `Rao()` as an alias for `paRao(alpha=2)`

## Minor changes
* Changes how `N` is intended in `mpaRaoAreaS`
* Other improvements in `paRaoS` and `paRaoP`
* Fixes bug in progBar for `paRao()`

# rasterdiv 0.3.4
## Major changes
* Adds back CRE function (leaving old function syntax)

## Minor changes
* Organises "accessory functions" in a separate file
* Documents almost all "accessory functions"
* Tries to improve progress bar
* Simplify heliPrep()
* Matches package with book chapter

# rasterdiv 0.3.3
## Major changes
* Adds heliPrep and heliPlot for helical graphs
* Temporarily removes CRE function
* Adds a vignettes for accRao (4)
* Adds a vignettes for helical plots (5)

## Minor changes
* Small changes in all vignettes
* Updated documentation to Roxygen style
* Uses terra in the place of raster
* Changes data type of copNDVI (SpatRaster now)
* Changes data type of world (SpatVector now)
* Improves documentation

# rasterdiv 0.3.2
## Minor changes
* Fixes minor "bug" for area based Rao that did not allow for polygons containing raster matrices with just one value.
* Adds an additional examples for Area based Rao.

# rasterdiv 0.3.1
## Minor changes
* Fixes repetition in paRao() manual

# rasterdiv 0.3.0
## Major changes
* Adds multidimension area-based Rao's index in *paRao()*.

## Minor changes
* Reorganises vignettes.
* Adds new GitHub site with *pkgdown()*.

# rasterdiv 0.2.5.2
## Minor changes
* Adds testthat tests for area-based Rao.

# rasterdiv 0.2.5.1
## Minor changes
* Made area-based *paRao* more efficient.
* Adapted vignettes to what changed.

# rasterdiv 0.2.5
## Major changes
* Starting to implement area-based *paRao*.

## Minor changes
* Vignettes reduced
* Added world vector dataset

# rasterdiv 0.2.4
## Major changes
* accRao now *RaoAUC*
* RaoAUC working with *paRao* in multidimensional mode.

## Minor changes
* Fixed some inconsistencies in the manual.

# rasterdiv 0.2.2p1
## Major changes
* Fixed issue due to high alphas which caused parametric Rao's to go to infinite (just added a warning for the multidimensional version).
* Rao still deprecated but connected in the background with *paRao(..., alpha=1)* for continuity with older versions.

# rasterdiv 0.2.2
## Major changes
* *Rao()* is now deprecated, replaced by *paRao(..., alpha=1)*.
* Multiple window sizes can now be indicated in paRao(..., window=c(3,5)).
* Multidimension paRao now working for np>1.

## Minor changes
* Vignette for multidimension Rao (rasterdiv_advanced_multimensionRao).
* Added *tests* folder for testthat.

# rasterdiv 0.2.1
## Major changes
* A new function **accRao** to derive the accumulation function of parametric Rao's index over a range of alphas. A vignette (rasterdiv_advanced_accRao) to show how to apply this function was added to the package.

# rasterdiv 0.2.0
## Major changes
* New parametrization of Rao's quadratic entropy. This parametrisation reconciles abundance-based indexes with distance-based indexes (Rao). The new function is called **paRao**, which stands for parametric Rao's entropy.

## Minor changes
* Minor improvements in the structure of all main functions.
* Option to save the output as a rasterstack for all main functions.
* pbapply and pbmcapply for progress bars.
* Main vignette (rasterdiv_basics) simplified.

# rasterdiv 0.1.0
First release which includes functions to compute on numerical matrices the following indexes:
* Shannon–Wiener index of diversity
* Rényi entropy
* Simpson index
* Berger-Parker index
* Hill numbers of order q or effective number of species
* Pielou evenness
* Rao's quadratic entropy
* Cumulative Residual entropy