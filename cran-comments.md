## Test environments

* Arch linux, R 3.6.0
* Ubuntu 16.04.6 (on travis-ci), R 3.6.0
* Fedora Linux (on rhub), R-devel, clang, gfortran
* Ubuntu Linux 16.04 LTS (on rhub), R-release, GCC
* Windows Server 2008 R2 SP1 (on rhub), R-devel, 32/64 bit

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking dependencies in R code ... NOTE
Missing or unexported objects:
  ‘pomp::covariate_table’ ‘pomp::gillespie_hl’

  These objects will be available from pomp version 2 and up. There is
  an if-statement that checks the installed version of pomp to
  determine whether they should be used.

## Downstream dependencies

There are currently no downstream dependencies for this package.
