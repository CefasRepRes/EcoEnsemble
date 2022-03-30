## R CMD check results

0 errors | 1 warning | 3 notes 

* This is a new release.

The warning for this package is:
 
* Warning: class "stanfit" is defined (with package slot 'rstan') but no metadata object found to revise superClass information---not imported?  Making a copy in package 'EcoEnsemble'

This warning is caused because one of my classes has a slot containing a "stanfit" object from the "rstan" package. However this object is not exported by "rstan" so it has no metadata available. These "stanfit" objects are returned by many "rstan" functions so there should be no issue using these objects as slots.
