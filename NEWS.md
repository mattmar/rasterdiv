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