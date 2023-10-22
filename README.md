# BICpack

The function `bic_oc` can be used for computing the BIC for a model with order constraints and one-sided constraints using the order-constrained BIC (Mulder & Raftery, 2019). The function requires a fitted model (e.g., glm, survival) as input as well as a string that specifies a set of order constraints on the regression coefficients. The function uses some functionality of [BFpack](https://cran.r-project.org/web/packages/BFpack/index.html).

Last Modified 10/22/23

Licensed under the GNU General Public License version 2 (June, 1991)


Installation
------------

You can install BICpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BICpack")
```

Basic example
-------------

``` r
library(BICpack)

# Consider a main unconstrained model
mtcars.standardized <- as.data.frame(scale(mtcars))
salfit <- glm(mpg ~ cyl + disp + hp, family = gaussian, data = mtcars.standardized)

# and a model selection problem between:
# M1: cyl < 0 & disp < 0
# M2: (cyl < 0 & disp < 0) (complement space of M1)
# M3: cyl = disp = 0

# Compute constrained BIC of model M1
bic_oc1 <- bic_oc(salfit, constraints = "cyl < 0 & disp < 0")
# Note that the same order-constrained model can be written using brackets
bic_oc(salfit, constraints = "(cyl , disp ) < 0")

# Compute constrained BIC of model M2
bic_oc2 <- bic_oc(salfit, constraints = "cyl < 0 & disp < 0", complement = TRUE)

# Compute BIC of model M3
salfit0 <- glm(mpg ~ hp, family = gaussian, data = mtcars.standardized)
bic_oc3 <- bic_oc(salfit0) #when the 'constraints' are omitted the standard bic is given

# get posterior probabilities of the three above models assuming equal prior model probabilities
bicvec <- c(bic_oc1[[1]], bic_oc2[[1]], bic_oc3[[1]])
postprob(bicvec) #largest posterior probability for multivariate one-sided model

# model selection of an order-constrained model versus its complement
salfit2 <- glm(mpg ~ gear + qsec + vs + am, family = gaussian, data = mtcars.standardized)
bic_oc4 <- bic_oc(salfit2, constraints = "am > vs > qsec")
bic_oc5 <- bic_oc(salfit2, constraints = "am > vs > qsec", complement = TRUE)
#more evidence for the order-constrained model than for its complement
postprob(c(bic_oc4[[1]],bic_oc5[[1]]))


```


Citing BICpack
------------

You can cite the package 

Mulder, J., & Raftery, A. E. (2022). BIC Extensions for Order-constrained Model Selection.
Sociological Methods & Research, 51(2) 471–498. [doi:10.1177/0049124119882459](https://journals.sagepub.com/doi/full/10.1177/0049124119882459)

