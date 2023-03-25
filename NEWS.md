# NScluster 1.3.6

* Removed C wrapper functions and registered entry points for the routines accessed by the `.Fortran` interface to call Fortran subroutines.

* Added a ‘`NEWS.md`’ file.


# NScluster 1.3.5

* Added JSS registration information to the ‘`DESCRIPTION`’ file and the help page references.

* Added a ‘`CITATION`’ file.


# NScluster 1.3.4

* Removed unused or undefined variables in Fortran COMMON blocks.


# NScluster 1.3.3

* Fixed LTO (Link-Time Optimization) problems reported in CRAN package check results.

* Changed the specification of `boot.mple()`, i.e., by setting of random seed
 reproduction was enabled.

* Added new argument `trace = FALSE` to `boot.mple()` to enable/disable progress bar.


# NScluster 1.3.1

* Removed `Simulatexxxxx()`, `Estimatexxxxx()` and `Palmxxxxx()` functions, and replaced them
 with the main functions `sim.cppm()`, `mple.cppm()` and `palm.mple()` with `model` argument for simulation, MPLEs
 (the maximum Palm likelihood estimates), nonparametric and parametric estimation of Palm intensity.

* Added a new function `boot.mple()` to perform the bootstrap by repeatedly simulating and fitting point patterns under the fitted model.

* Added S3 methods `plot()`, `coef()` and `summary()` for returned class objects.

* Added AIC calculations to evaluate models in the `mple.cppm()` function.

* Added defaults for initial values of `mple.cppm()`.

* Revised the guide so that it can be referenced in `vignette("NScluster")`.


# NScluster 1.2.0

* Changed the implementation of the random number generator to Mersenne-Twister.

* Changed the names of simplex functions to `EstimateIP()`, `EstimateThomas()`, `EstimateTypeA()`, `EstimateTypeB()` and `EstimateTypeC()`.

* Changed the range of x and y labels in the output plots of the `Estimatexxxxx()` functions.

* Removed `ty` arguments from `PalmIP()`, `PalmThomas()`, `PalmTypeA()`, `PalmTypeB()` and `PalmTypeC()`.
  Fixed `ty = 1` in R programs.

* Renamed arguments `x2` to `uplimit = 0.3` of functions `EstimateIP()`, `EstimateTypeA()`, `PalmIP()` and `PalmTypeA()`.

* Renamed arguments `offspring` to `data` of Estimate functions and Palm functions.

* Added default values `0.001` for `delta` of Palm functions.

* Instead of removing the `pa` argument of the Palm function, the vector of true parameters `par1 = NULL` and the vector of MPLEs `par2 = NULL` have been added.

* Added new arguments `parents.distinct = FALSE` to `SimulateTypeB()` and `SimulateTypeC()` 
 to distinguish between two groups of plotted points.

* Added ‘`src/init.c`‘.


# NScluster 1.1.1

* Merged **NSclusterOMP** package (parallel versions of parameter estimation by the simplex method) into **NScluster** package.

* Fixed Fortran code according to the warning message when using gfortran with `-Wall -pedantic`.

* Fixed Simplex-IPfp.f, Simplex-TypeAfp.f, Plam-IPf.f and Palm-TypeAf.f for error handling.

* Corrected `True value` to `Initial value` in the expression of simplex function outputs.

* Removed `MPLE : simplex` and `normalized` outputs from simplex functions.

* Added legends to plots of simplex functions.


# NScluster 1.1.0

* Corrected `PKG_LIBS` in ‘`NScluster/src/Makevars`’ to `$(SHLIB_OPENMP_CFLAGS)`. (Reported by Brian Ripley.)

* Removed the `setOmpNumThread()` function.


# NScluster 1.0.2

* Corrected ‘`inst/doc/index.html`’ according to the error message for the W3C Markup Validation Service.

* Reduced example in `SimplexThomas()` to less than 5 sec.


#  NScluster 1.0.1

* Removed ‘`inst/LICENCE`’ and added ‘`inst/OUTHORS`’.
