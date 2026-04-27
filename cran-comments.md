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

# Version 1.0.4

* Changed `EcoEnsemble.Rmd` so that it will be compatible with a future `ggplot2` release.

## R CMD check results
There were no ERRORs or WARNINGS

There were 3 NOTES

* checking for GNU extensions in Makefiles ... NOTE GNU make is a SystemRequirements.

This is a consequence of using Rstan which requires us to use GNU make to compile the Stan model as per their package usage instructions


* installed size is 99.8Mb
  sub-directories of 1Mb or more:
    libs  98.7Mb

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

# Version 1.0.5

* Changed `EcoEnsemble` so that it will be compatible with a future the recent `rstan` release (2.26.0).

## R CMD check results

* For the win release there was one NOTE

checking CRAN incoming feasibility ... NOTE Maintainer: 'Michael A. Spence <michael.spence@cefas.gov.uk>'
Found the following (possibly) invalid URLs:
  URL: https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310
    From: DESCRIPTION
          man/EcoEnsemble-package.Rd
          README.md
    Status: 403
    Message: Forbidden
    
The URL remained from the previous CRAN submission. The URL links to the correct website.

* The Mac release gave an error:

checking package dependencies ... ERROR
Packages required and available but unsuitable versions: 'rstan', 'rstan'

I discussed this issue here https://github.com/CefasRepRes/EcoEnsemble/pull/2, where it was proposed that the problem was with the binaries on the Mac release check and it would sort itself out when they are updated. I presume that once updated there will be two Notes, like in all the other versions.

* From previous submission I have changed:

The Description field contains
     models"<https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12310>.
   Please use permanent DOI markup for linking to publications as in <doi:prefix/suffix>.

to:

models: models"<doi:10.1111/faf.12310>.

# Version 1.0.6

* Extended `get_mcmc_ensemble_model` to return all stanmodels in package

#Version 1.1.0

For the Mac release there were two NOTEs:

* NOTE GNU make is a SystemRequirements.

This is a consequence of using Rstan which requires us to use GNU make to compile the Stan model as per their package usage instructions. 

* NOTE installed size is 219.1Mb
  sub-directories of 1Mb or more:
    doc     3.5Mb
    libs  214.6Mb

This is a consequence of using Rstan which produces quite large binaries when compiled.

From the previous submission I have changed:

Added new prior option and the option to include drivers. Added vignettes 'BetaConjugatePrior' and 'IncludingDrivers'. Also added citation file.

#Version 1.1.1

For the Mac release there were two NOTEs:

* NOTE GNU make is a SystemRequirements.

This is a consequence of using Rstan which requires us to use GNU make to compile the Stan model as per their package usage instructions. 

* NOTE installed size is 219.1Mb
  sub-directories of 1Mb or more:
    doc     3.5Mb
    libs  214.6Mb

This is a consequence of using Rstan which produces quite large binaries when compiled.

From the previous submission I have changed:

Added new prior option and the option to include drivers. Added vignettes 'BetaConjugatePrior' and 'IncludingDrivers'. Also added citation file.

#Version 1.2.0

Added alternative sampler. Required new Rstan files, explanation vignettes and tests.

Mac check receives WARNING. See below.

Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/Volumes/Builds/templates/cran-build/sonoma-arm64/Rlib/4.6/QuickJSR/libs/QuickJSR.so':
  dlopen(/Volumes/Builds/templates/cran-build/sonoma-arm64/Rlib/4.6/QuickJSR/libs/QuickJSR.so, 0x0006): Symbol not found: _R_getRegisteredNamespace
  Referenced from: <70FA6501-FF39-348D-9485-5F13E9D957C9> /Volumes/Builds/templates/cran-build/sonoma-arm64/Rlib/4.6/QuickJSR/libs/QuickJSR.so
  Expected in:     <9A4FD863-B5A1-33CC-BDC9-62F15EEC56CE> /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libR.dylib
Calls:  ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
Execution halted
ERROR: configuration failed for package ‘EcoEnsemble’
* removing ‘/Volumes/PkgBuild/work/1777288809-78f894df3ff947c2/packages/sonoma-arm64/results/4.6/EcoEnsemble.Rcheck/EcoEnsemble’

The warning did not occur for check_mac_release() before R 4.6 release. Since, only minor vignette changes have occurred. Further, checking previous CRAN-accepted versions of EcoEnsemble receive the same warning. 
Therefore, I think this error does not show an issue with the code. 


