# Unsupervised REgression MIXtures (uReMix) 

This is a branch from chamroukhi/uReMix_R, which does not provide an R package.

This longitudinal data clustering method relies on the assumption of normal error and ignore the within-subject correlations during the model set-up. The mixtures are all linear models, and can be solved straightforwardly by EM algorithm.

The most appealing feature of this method is it can determine the number of clusters automatically through an Robust EM framework. As for most model-based longitudinal data clustering algorithms, people need to pre-specify the number of clusters ahead of time to obtain the expression of likelihood function.

I adjusted the codes and wrap it as an R package so that user could directly install it by 

```r
devtools::install_github("https://github.com/junyzhou10/uReMix_R/tree/packageVersion")
```
and apply it easily.

# Usage

The main function is `emSRM(X, Y, ...)`, where X is the observation time which should be stored in long form, and Y is the observations stored in wide form. The observation time X needs to be regular, i.e. the total number of observations are the same and aligned for each subjects.

---

Updates: I created another function `emSRMir(X, Y, id, ...)` which is able to handle irregular longitudinal data, i.e., the nubmer of observations could be different among subjects and the observation time points can be different across individuals. Here the `X` is long form, but `Y` is also long form as well. The observations from any individual are identified by input `id`. The rest function arguments are the same as `emSRM(X, Y, ...)`.


