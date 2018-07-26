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
	
Changes 2012-7-19   Erika Cule

  * fixed a bug in linearRidge when scaling = "none"
  
  * added functions linearRidgeGenotypes and logisticRidgeGenotypes and their predicting counterparts
    linearRidgeGenotypesPredict and logisticRidgeGenotypespredict. 
    These functions fit linear and logistic ridge regression models for genome-wide SNP data, 
    optionally automatically choosing the ridge parameter
   
	* minor bug fix: "Intercept" now prints as "(Intercept)"
		(as is the case for lm and glm models)
	
