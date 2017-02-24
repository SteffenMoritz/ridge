#Adopting Orphaned Package
I would like to adopt the ridge package.
(together with Lingbing Feng I want to revive his imputeR package,
which was archived, because ridge, which it imported was also archived)

I actually would have prefered, if the original author (Erica Cule) brought it back to CRAN.
Her university mail dropped, that is why the package was orphaned.
I did some research about her and tried to contact her via Facebook and ResearchGate some weeks ago.
(but thus far I got no answer)
If she maybe answers back in the future and is interested in maintaining the package again, 
I'll be more than happy to give maintainer rights back to her.

Most of the changes I made to the package were to improve useability + documentation and to comply to CRAN checks again.



## Test environments
* local OS X install, R 3.3.2
* ubuntu 12.04 (on travis-ci), R 3.3.2
* win-builder (devel and release)

## R CMD check results

I get from local OS X & ubuntu travis-ci:
0 errors | 0 warnings | 0 notes


I get from win-builter:
0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Steffen Moritz <steffen.moritz10@gmail.com>'

-----> As stated above I would like to adopt the orphaned package


Possibly mis-spelled words in DESCRIPTION:
  SNP (5:57)
  nucleotide (5:32)
  polymorphism (5:43)
  
-----> These words are technical, but spelled correctly
