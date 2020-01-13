
<img src="vignettes/figures/NR-logo_utvidet_r32g60b136_small.png" style="position:absolute;top:0px;right:0px;" align="right" height="30px"/>


# The rWHAP package
### Wave Height Analysis and Prediction
#### Authors: Hugo L. Hammer, Tor Arne Oigard, Thordis Thorarinsdottir and Hanne Rognebakke

***

## Overview

### Installation and loading the rWHAP package
The rWHAP package depends on the following R-packages

  - *glmnet*
  - *forecast*
  - *moments*

<!-- end list -->

These dependencies will be automatically installed when installing the rWHAP package.

The most recent version of the *rWHAP* package is hosted on a git repository at
<https://github.com/NorskRegnesentral/rWHAP.git>.

In order to install the rWHAP package run the following command:

```
devtools::install_git(url = "https://github.com/NorskRegnesentral/rWHAP.git")
```

### Example data
A small example data set is included. It contains Significant Wave Height (SWH), Sea Level Preassure (SLP) and gradients for the SLP for the years 2006 to 2015. The data is within the latitude range of 53.25 - 57.75 and longitude range of  54.00 - 58.50. The time resolution is 6 hours, i.e., 4 measurements available for each day.

THORDIS: Kanskje vi bør skrive mer om data? Nevne Era Interim? Sjekk med paperet!!!

The table below lists the various variables in the data set.

The data can be loaded using

```r
load(SWHdata)
```

To see summary statistics of the SWH in the data:
```r
summary(SWH)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4726  2.1296  3.0538  3.4142  4.3071 16.5310
```

For summary statistics of the SLP in the data:
```r
summary(SLP)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   93301  100127  101100  100975  101975  104539
```

Summary statistics for the gradient of the SLP:
```r
summary(SLP.grad)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.000000 0.004531 0.011082 0.017110 0.022871 0.538579
```


### Model fitting
