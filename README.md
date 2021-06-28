# Unsupervised REgression MIXtures (uReMix) 

This is a branch from chamroukhi/uReMix_R, which does not provide an R package.

I adjusted the codes and put it in the R package framework so that user could directly install it by 

```r
devtools::install_github("https://github.com/junyzhou10/uReMix_R/tree/packageVersion")
```


# Usage

The main function is `emSRM(X, Y, ...)`, where X is the observation time which should be stored in long form, and Y is the observations stored in wide form. The observation time X needs to be regular, i.e. the total number of observations are the same and aligned for each subjects.

---

Now I created another function `emSRM(X, Y, id, ...)` which is able to handle irregular longitudinal data. 
