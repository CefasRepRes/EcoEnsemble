## R CMD check results

0 errors | 0 warning | 1 notes 

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Resubmission 1

This is a resubmission. In this version I have:

* Inlcuded an additional author for the package.
* Removed the redundant "in R" from the package description.
* Ensured options() and par() are reset after being run in examples and vignettes.

## Resubmission 2

This is a resubmission. The package was archived on CRAN on 2022-10-24 as issues were not corrected in time.

In this version we have:

* Changed the maintainer to michael.spence@cefas.gov.uk for james.martindale@cefas.co.uk. This is because James has changed career path.
* Changed the EcoEnsemble_types.h file so that it passes CRAN checks on a Mac.

## R CMD check results
There were no ERRORs or WARNINGS

There were 3 NOTES

* checking for GNU extensions in Makefiles ... NOTE GNU make is a SystemRequirements.

This is a consequence of using Rstan which requires us to use GNU make to compile the Stan model as per their package usage instructions

* checking dependencies in R code ... NOTE Namespaces in Imports field not imported from: ‘RcppParallel’ ‘rstantools’ All declared Imports should be used.

This is a false positive. These packages are used by the makefile which compiles the Stan code when the package is installed. The files containing these libraries are generated on the fly when the package is installed and are thus not visible to R CMD CHECK. This process is controlled by the Stan developers.

* checking installed package size ... NOTE installed size is 119.8Mb sub-directories of 1Mb or more: libs 118.8Mb

This is a consequence of using Rstan which produces quite large binaries when compiled.

Thanks to the https://github.com/insightsengineering/rbmi/blob/main/cran-comments.md for helping describe these notes.

## Resubmission 3

* Changed http://mc-stan.org/ to https://mc-stan.org/ in vignette.

# Version 1.0.2

* Changed .Rbuildignore so that it will be compatible with a future 'rstan' release.

## R CMD check results
There were no ERRORs or WARNINGS

There were 3 NOTES

* checking for GNU extensions in Makefiles ... NOTE GNU make is a SystemRequirements.

This is a consequence of using Rstan which requires us to use GNU make to compile the Stan model as per their package usage instructions


* checking installed package size ... NOTE installed size is 120.3Mb sub-directories of 1Mb or more: libs 119.3Mb

This is a consequence of using Rstan which produces quite large binaries when compiled.

Thanks to the https://github.com/insightsengineering/rbmi/blob/main/cran-comments.md for helping describe these notes.

* checking CRAN incoming feasibility ... NOTE Maintainer: 'Michael A. Spence <michael.spence@cefas.gov.uk>'
Found the following (possibly) invalid URLs:
  URL: https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310
    From: DESCRIPTION
          man/EcoEnsemble-package.Rd
          README.md
    Status: 403
    Message: Forbidden
    
The URL remained from the previous CRAN submission. The URL links to the correct website.

# Version 1.0.3

* Changed `EcoEnsemble.Rmd` so that it will be compatible with a future `ggplot2` release.

## R CMD check results
There were no ERRORs or WARNINGS

There were 3 NOTES

* checking for GNU extensions in Makefiles ... NOTE GNU make is a SystemRequirements.

This is a consequence of using Rstan which requires us to use GNU make to compile the Stan model as per their package usage instructions


* checking installed package size ... NOTE installed size is 120.3Mb sub-directories of 1Mb or more: libs 119.3Mb

This is a consequence of using Rstan which produces quite large binaries when compiled.

Thanks to the https://github.com/insightsengineering/rbmi/blob/main/cran-comments.md for helping describe these notes.

* checking CRAN incoming feasibility ... NOTE Maintainer: 'Michael A. Spence <michael.spence@cefas.gov.uk>'
Found the following (possibly) invalid URLs:
  URL: https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310
    From: DESCRIPTION
          man/EcoEnsemble-package.Rd
          README.md
    Status: 403
    Message: Forbidden
    
The URL remained from the previous CRAN submission. The URL links to the correct website.
