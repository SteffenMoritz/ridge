## Changes Version 3.3 -  Steffen Moritz

  * Fixed small bug with a unwanted print statement in linearRidge

## Changes Version 3.2 -  Steffen Moritz

  * Fixed issue with predict() function:
    Wrong results, when predicting the training data without newdata argument.
    When supplying them via newdata, everything seemed fine.
    (https://github.com/SteffenMoritz/ridge/issues/16)
    Thanks to Mekala Sundaram for reporting the issue.

## Changes Version 3.1 -  Steffen Moritz

  * Fixes to remain on CRAN
  
  * GitHub Actions as CI Tool 


## Changes Version 2.7 until 3.0 -  Steffen Moritz

  * Fixes related to autoconf (needed for passing new CRAN checks)


## Changes Version 2.6 -  Steffen Moritz

  * Readme update

  * Fixes to pass CRAN checks

## Changes Version 2.5 -  Dan Frankowski

  * Fix predict bug with factor variables

  * Add vcov.ridgeLinear

  * Add testthat tests
  
  * Changed documentation partially to roxygen2

## Changes Version 2.4 -  Steffen Moritz

  * bugfix in logisticRidge and linearRidge see: https://github.com/SteffenMoritz/ridge/issues/2
  
  * bugfix to pass cran checks 'warning: missing template: HAVE_GSL_HEADER autoheader'
  
  * Updated Readme


## Changes Version 2.3 -  Steffen Moritz

  * Some fixed to remain on CRAN and comly with CRAN policy
  
  * Improved Description file
  
  * Changed NEWS file to markup document
  

## Changes Version 2.2 -  Steffen Moritz

  * Made package CRAN ready again
  
  * Created github repository for the package
  
  * Fixed warning using Wpendantic gcc
  
  * Adapted the DESCRIPTION file to latest CRAN requirements
  
  * Renamed CHANGELOG to NEWS
  
  * NEWS template update
  
  * Changes to .C Method registration
  
  * NAMESPACE fixes
  

## Changes 2014-3-02 -	Erika Cule

	* Fixed layout of .Rd files
	
	* Added deletion of Makevars to cleanup script
	

## Changes 2012-9-27 -	Erika Cule 

	* Flat text (.txt) data files were moved from ridge/data to ridge/inst/extdata 
  	(in the source, which becomes ridge/extdata in the installed package). The .txt files 
  	should be in inst/extdata because they files are used by the package examples 
	  (albeit in not run sections), as described in Writing R Extensions 1.1.5 Data in packages.
	
	* Some users were reporting problems installing the package on some Linux OS. configure 
	  has been modified to fix this 
	  problem.


## Changes 2012-8-21 -   Erika Cule 

	* Bug fix in src/commonFunctions.c


## Changes 2012-8-21 -  Erika Cule

	* Added configure.ac so that package will install if GSL >= 1.14 is not available 
	  (with linearRidgeGenotypes, logisticRidgeGenotypes, linearRidgeGenotypesPredict and 
	  logisticRidgeGenotypesPredict disabled). 
	  
	* configure.ac detects whether openblas is available and if it is found, links to that. 
  	This speeds up computation. (http://xianyi.github.com/OpenBLAS/)
	
## Changes 2012-7-19   Erika Cule

  * fixed a bug in linearRidge when scaling = "none"
  
  * added functions linearRidgeGenotypes and logisticRidgeGenotypes and their predicting counterparts
    linearRidgeGenotypesPredict and logisticRidgeGenotypespredict. 
    These functions fit linear and logistic ridge regression models for genome-wide SNP data, 
    optionally automatically choosing the ridge parameter
   
	* minor bug fix: "Intercept" now prints as "(Intercept)"
		(as is the case for lm and glm models)
	
