## Update for publication

17/04/18
* removeOutliers has been made a optional process for survivALL calculations

## Update & bug-fix

01/03/18
* some computations were returning incorrect values inside suppressWarnings() - removed
* HRs and bootstrap results now calculate a bootstrap p-value in addition to the results of coxph
  * this has now been sped up ~100x
* function to calculate significance using a continuous variable added
* vignettes have been re-worked to reflect package updates

## Update

04/10/17
* further vignette improvements
* suspected spelling mistakes are unchanged from previous accepted version 

## Update

24/08/17
* updated plotALL() function to allow multivariate analysis
* updated vignettes


## Resubmission

11/08/17
* added references to DESCRIPTION
* enclosed package name (survivALL) in '' in DESCRIPTION
* formatted DESCRIPTION to 80 characters wide

10/08/17
* corrected title font
* tidied up vignettes
* sped up examples
* removed breastCancerNKI dependency that was failing to install

## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* win-builder (release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs :)

## Possible mis-spellings
IHC - initialism for immunohistochemistry, correctly spelled
Biomarker/biomarker - correctly spelled
survivALL - the name of the package
