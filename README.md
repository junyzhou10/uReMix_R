# Unsupervised REgression MIXtures (uReMix) 

This is a branch from chamroukhi/uReMix_R, which does not provide an R package.

I adjusted the codes and put it in the R package framework so that user could directly install it by 

```r
devtools::install_github("https://github.com/junyzhou10/uReMix_R/tree/packageVersion")
```


# Usage

The main function is `emSRM(X, Y, ...)`, where X is the observation time which should be stored in long form, and Y is the observations stored in wide form. Currently, the X needs not to be regular, but the number of observations are required to be the same.
