# Unsupervised REgression MIXtures (uReMix) 

This is a branch from chamroukhi/uReMix_R, which does not provide an R package.

I adjusted the codes and wrap it as an R package so that user could directly install it by 

```r
devtools::install_github("https://github.com/junyzhou10/uReMix_R/tree/packageVersion")
```
and apply it easily.

# Usage

The main function is `emSRM(X, Y, ...)`, where X is the observation time which should be stored in long form, and Y is the observations stored in wide form. The observation time X needs to be regular, i.e. the total number of observations are the same and aligned for each subjects.

---

Updates: I created another function `emSRMir(X, Y, id, ...)` which is able to handle irregular longitudinal data, i.e., the nubmer of observations could be different among subjects and the observation time points can be different across individuals. Here the `X` is long form, but `Y` is also long form as well. The observations from any individual are identified by input `id`. The rest function arguments are the same as `emSRM(X, Y, ...)`.


