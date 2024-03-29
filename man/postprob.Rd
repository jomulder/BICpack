% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIC_orderconstraints.R
\name{postprob}
\alias{postprob}
\title{Compute posterior model probabilities for given BIC values}
\usage{
postprob(bic, priorprob = 1)
}
\arguments{
\item{bic}{A vector of length greater than 1 containing BIC-values of the models of interest.}

\item{priorprob}{A vector containing the prior probabilities of the models. The vector does
not have to be normalized. The default setting gives equal prior probabilities.}
}
\value{
A vector of posterior model probabilities
}
\description{
\code{postprob} computes posterior model probabilities based on given \code{bic}-values
assuming prior probabilities specified as \code{priorprob}.
}
\examples{
\donttest{

# compute bic of an order-constrained model
salfit <- glm(mpg ~ gear + qsec + vs + am, family = gaussian, data = mtcars.standardized)
bic_oc1 <- bic_oc(salfit, constraints = "am > vs > qsec")
# compute bic of the complement of this order-constrained model
bic_oc2 <- bic_oc(salfit, constraints = "am > vs > qsec", complement = TRUE)
postprob(c(bic_oc1[[1]],bic_oc2[[1]]))
}
}
