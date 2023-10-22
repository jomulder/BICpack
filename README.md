# BICpack

The function `bic_oc` can be used for computing the BIC for a model with order constraints and one-sided constraints using the order-constrained BIC (Mulder & Raftery, 2019). The function requires a fitted model (e.g., glm, survival) as input as well as a string that specifies a set of order constraints on the regression coefficients. The function uses some functionality of **BFpack** (Mulder et al., 2021, https://cran.r-project.org/web/packages/BFpack/index.html).

BICpack was written by Joris Mulder <j.mulder3@tilburguniversity.edu>

Last Modified 10/22/23

Licensed under the GNU General Public License version 2 (June, 1991)

Basic example
-------------

``` r
library(BICpack)

# Testing a model assuming a positive gender effect versus a model assuming a negative
# gender effect versus a model assuming no gender effect.

# fit a model with all effects present
mtcars.standardized <- as.data.frame(scale(mtcars))
salfit <- glm(mpg ~ cyl + disp + hp, family = gaussian, data = mtcars.standardized)
# a model which assumes a negative effect of 'cyl' and 'disp'
bic_oc1 <- bic_oc(salfit, constraints = "cyl < 0 & disp < 0")
# same model with a different, equivalent notation of the constraints using brackets
bic_oc(salfit, constraints = "(cyl , disp ) < 0")

# fit a model which encompasses the complement of the constrained space of 'cyl' and 'disp'
bic_oc2 <- bic_oc(salfit, constraints = "(cyl , disp ) < 0", complement = TRUE)

# fit a model without 'cyl' and 'disp'
salfit0 <- glm(mpg ~ hp, family = gaussian, data = mtcars.standardized)
bic_oc3 <- bic_oc(salfit0) #when the 'constraints' are omitted the standard bic is given

# get posterior probabilities of the three above models assuming equal prior model probabilities
bicvec <- c(bic_oc1[[1]], bic_oc2[[1]], bic_oc3[[1]])
postprob(bicvec) #largest posterior probability for multivariate one-sided model

# compute bic of an ordered model
salfit2 <- glm(mpg ~ gear + qsec + vs + am, family = gaussian, data = mtcars.standardized)
bic_oc(salfit2, constraints = "am > vs > qsec")
bic_oc(salfit2, constraints = "am > vs > qsec", complement = TRUE)
#more evidence for the order-constrained model than for its complement
```

Installation
------------

You can install BICpack from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BICpack")
