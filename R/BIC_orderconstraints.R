
# The function 'bic_oc' can be used for computing the order-constrained BIC
# for a model-object (e.g., glm, coxph) with additional order constraints on
# certain parameters. The order-constrained BIC is based on the local unit-
# information prior. Please refer to Mulder & Raftery's "BIC extensions for
# order-constrained model selection."
#
# This package was written by Joris Mulder, MTO, Tilburg University, 2018, in
# collaboration with Anton Olsson Collentine as part of his internship at MTO.

#' Compute the BIC of an order-constrained model.
#'
#'\code{bic_oc} computes the order-constrained BIC for a fitted model object (e.g., lm-object,
#'a glm-object, or a coxph-object) with order constraints on certain effects.
#'
#' @param object A fitted model object, such as a glm-object
#' @param constraints A string specifying order constraints on certain effects of the modeling
#' object. The ampersant (\code{&}) can be used to specify separate constraints.
#' If the value is NULL the ordinary BIC is computed.
#' @param complement A logical scalar that specifies if the order-constrained subspace is
#' considered (FALSE) or its complement (TRUE). Default is FALSE.
#' @param N The sample size that was used to fit the model \code{object}. If it is left empty
#' (\code{NULL}), then the sample size will be determine based on the dimension of the
#' \code{fitted.values} element of \code{object}.
#' @return Return a list of which the order-constrained BIC of modeling \code{object} with
#' \code{constraints} based on the local unit-information prior as first element, the posterior
#' probability that the constraints hold as second element, and the prior probability that the
#' constraints hold as third element. If \code{complement} is TRUE then the complement of the
#' order-constrained subspace is considered.
#'
#' @references Mulder, J., and Raftery, A.E. (2022). BIC Extensions for Order-constrained Model
#' Selection. Sociological Methods & Research, 51 (2), 471-498. <DOI:10.1177/0049124119882459>
#'
#' @examples
#' \donttest{
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + .5 * x1 + 1 * x2 + rnorm(n)
#' df <- data.frame(y=y,x1=x1,x2=x2)
#' glm1 <- lm(y ~ 1 + x1 + x2, data=df)
#' # Compute the BIC of the fitted model `glm1' with order constraints that state that
#' # the effect of 'x2' on 'y' is larger than the effect of 'x1' on 'y', and both effects
#' # are assumed positive.
#' bic_oc(glm1,"x2 > x1 > 0")
#' # the same result would be obtained by separating the constraints with a '&'
#' bic_oc(glm1,"x2 > x1 & x1 > 0")
#'
#' # a model where both coefficients are assumed to be positive
#' bic_oc(glm1,"x2 > & x1 > 0")
#' # the same model where both coefficients are assumed to be positive using the brackets notation
#' bic_oc(glm1,"(x2 , x1 ) > 0")
#' }
#' @export
bic_oc <- function(object, constraints=NULL, complement=FALSE, N=NULL){

  # The function 'bic_oc' can be used for computing the order-constrained BIC
  # for a model-object (e.g., glm, coxph) with additional order constraints on
  # certain parameters. The order-constrained BIC is based on the local unit-
  # information prior.

  #compute unconstrained log likelihood
  logLike <- logLik(object)
  numpara <- attr(logLike,"df")
  logLike.unc <- logLik(object)[1]
  #
  if(is.null(N)){
    N <- length(object$fitted.values)
    if(class(object)[1]=="coxph"){
      N <- object$n
    }
    #message(paste0("The sample size was extracted from the object, resulting in N = ",
    #              as.character(N),", if this","\n",
    #              "is not correct specify the sample size manually via the 'N' argument.","\n"))
  }

  if(is.null(constraints)){ #compute regular BIC
    margLike <- margLike_unc <- logLike.unc - numpara/2*log(N)
    postprob <- 1
    priorprob <- 1
  }else{
    # check that constraints do not contain equalities
    if(grepl("=", constraints, fixed = TRUE)){
      stop("The 'constraints' argument should not contain equalities '='.")
    }

    # to get posterior and prior probabilities use BF function from BFpack
    estimates <- object$coefficients
    Sigma <- vcov(object)[1:length(estimates),1:length(estimates)]
    approxBF_order <- BFpack::BF(x=estimates,hypothesis=constraints,Sigma=Sigma,n=N)
    postprob <- approxBF_order$BFtable_confirmatory[1,4]
    priorprob <- approxBF_order$BFtable_confirmatory[1,2]

    #invert for complement if specified
    if(complement==T){
      postprob <- 1 - postprob
      priorprob <- 1 - priorprob
    }

    #regular BIC of a model excluding the inequality constraints
    margLike <- logLike.unc - numpara/2*log(N) + log(postprob) - log(priorprob)
    margLike_unc <- logLike.unc - numpara/2*log(N)
    names(margLike) <- NULL
  }

  BIC_OC <- -2 * margLike
  BIC_unc <- -2 * margLike_unc
  call1 <- match.call()
  out <- list(BIC_OC=BIC_OC,BIC_unc=BIC_unc,postprob=postprob,priorprob=priorprob,
              constraints=constraints,call=call1)
  class(out) <- "BIC_OC"
  return(out)
}

#' @method print BIC_OC
#' @export
print.BIC_OC <- function(x,
                         digits = 3,
                         na.print = "", ...){

  cat("Call:")
  cat("\n")
  print(x$call)

  cat("\n")
  cat("Order-constrained BIC:","\n", sep = "")
  cat(as.character(round(x$BIC_OC,digits)),"\n", sep = "")

  cat("\n")
  cat("BIC (ignoring the constraints):","\n", sep = "")
  cat(as.character(round(x$BIC_unc,digits)),"\n", sep = "")

  cat("\n")
  cat("Constraints:","\n", sep = "")
  cat(x$constraints,"\n", sep = "")

  cat("\n")
  cat("Posterior probability that the constraints hold under an unconstrained model:","\n", sep = "")
  cat(as.character(round(x$postprob,digits)),"\n", sep = "")

  cat("\n")
  cat("Prior probability that the constraints hold under an unconstrained model:","\n", sep = "")
  cat(as.character(round(x$priorprob,digits)),"\n", sep = "")

}

#' Compute posterior model probabilities for given BIC values
#'
#'\code{postprob} computes posterior model probabilities based on given \code{bic}-values
#'assuming prior probabilities specified as \code{priorprob}.
#'
#' @param bic A vector of length greater than 1 containing BIC-values of the models of interest.
#' @param priorprob A vector containing the prior probabilities of the models. The vector does
#' not have to be normalized. The default setting gives equal prior probabilities.
#' @return A vector of posterior model probabilities
#' @examples
#' \donttest{
#' BIC <- c(2343.23,2344.99,2349.12)
#' #Compute the posterior probabilities given these
#' #BIC-values and assuming equal prior model probabilities.
#' postprob(BIC)
#' }
#' @export
postprob <- function(bic,priorprob=1){
  #compute posterior model probabilities based on the BIC
  #assuming having prior probabilities 'priorprob'

  if(!is.numeric(bic)) stop("bic's have to be numeric values")
  if(length(bic)==1) stop("multiple bic's are needed")
  if(min(bic)<0) stop("bic's cannot be negative")
  if(!is.numeric(priorprob)) stop("prior probabilities are not numeric values")
  if(length(priorprob)!=length(bic) && length(priorprob)!=1) stop("number of prior probabilities does not match number of bic's")
  if(length(priorprob)==1) priorprob = 1
  priorprob <- priorprob/sum(priorprob) #normalize prior probabilities if necessary
  logmarglike <- -bic/2
  marglike_scaled <- exp(logmarglike - max(logmarglike))
  postprob <- marglike_scaled*priorprob/sum(marglike_scaled*priorprob)
  return(round(postprob,3))
}

#' @importFrom stats logLik vcov
#' @importFrom BFpack BF

create_matrices_oc <- function(object, constraints){
  # This is a slight modification of a function developed by Anton Olsson Collentine
  # as part of his internship at MTO, Tilburg University (April 2018).
  # The functions translates a string with order constraints to an augmented matrix
  # with coefficients for the order constraints [RI | rI], where RI %*% effects > rI.

  varnames <- names(object$coefficients) #provides the variable names of the linear model object, including intercept
  if(is.null(varnames)) stop("Please input proper linear model object")

  hyp2 <- gsub(" ", "", constraints) #removes all whitespace
  if(!grepl("^[0-9a-zA-Z><,().-]+$", hyp2)) stop("Impermissable characters in hypotheses") #Self-explanatory. NEW parentehese
  if(grepl("[><]{2,}", hyp2)) stop("Only use '>' or '<' for specifying order constraints")

  step1 <- unlist(strsplit(hyp2, split = "[<>,()]")) #split by special characters and unlist
  input_vars <- step1[grep("[a-zA-Z]+", step1)] #extract subunits that contain at least one letter
  if(!all(input_vars %in% varnames)) stop("Hypothesis variable(s) not in object, check spelling") #Checks if input variables exist in model-object

  pos_comparisons <- unlist(gregexpr("[<>=]", hyp2)) #Gives the positions of all comparison signs
  leftside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
  rightside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
  pos1 <- c(-1, pos_comparisons) #positions to extract data to the leftside of comparisons
  pos2 <- c(pos_comparisons, nchar(hyp2) + 1) #positions to extract data to the rightside of comparisons
  for(i in seq_along(pos1)){
    leftside[i] <- substring(hyp2, pos1[i] + 1, pos1[i+1] - 1) #Extract all variables or outcomes to the leftside of a comparison sign
    rightside[i] <- substring(hyp2, pos2[i] + 1, pos2[i+1] - 1) #Extract all variables or outcomes to the rightside of a comparison sign
  }
  leftside <- leftside[-length(leftside)] #remove last element which is a NA due to loop formatting
  rightside <- rightside[-length(rightside)] #remove last element which is a NA due to loop formatting
  comparisons <- substring(hyp2, pos_comparisons, pos_comparisons) #Extract comparison signs
  framed <- data.frame(left = leftside, comp = comparisons, right = rightside, stringsAsFactors = FALSE) #hypotheses as a dataframe


  if(any(grepl(",", framed$left)) || any(grepl(",", framed$right))){ #Larger loop that deals with commas if the specified hypothesis contains any
    if(nrow(framed) > 1){
      for(r in 1:(nrow(framed)-1)){ #If a hypothesis has been specified with commas e.g., "X1 > 0, X2 > 0" or "(X1, X2) > X3"
        if(all.equal(framed$right[r], framed$left[r+1])){ #The right hand side of the hypothesis df will be equal to the next row left side
          if(substring(framed$right[r], 1, 1) == "(") { #If the first row begins with a ( as when "X1 > (X2, X3)" and opposed to "(X2, X3) > X1"
            framed$right[r] <- sub("),.+", ")", framed$right[r])#If so, remove everything to the right of the parenthesis on the right hand side
            framed$left[r+1] <- sub(".+),", "", framed$left[r +1])#and everything to the left of the parenthesis on the left hand side to correct the df
          } else{
            framed$right[r] <- sub(",.+", "", framed$right[r]) #else, remove everything to the right of the comma on the right hand side
            framed$left[r+1] <- sub("[^,]+,", "", framed$left[r+1]) #and everything to the left of the comma on the left hand side to correct the df
          }
        }
      }
    }

    commas_left <- framed$left[grep(",", framed$left)] #At this point all remaining elements that contain commas should also have parentheses, check this
    commas_right <- framed$right[grep(",", framed$right)] #Necessary to use is isTRUE below in case one of these contains no commas, and 'any' for several rows
    if(isTRUE(any(!grepl("\\(.+)", commas_left))) || isTRUE(any(!grepl("\\(.+)", commas_right))) || #Check so rows contain parenthesis
       isTRUE(any(grepl(").+", commas_left))) || isTRUE(any(grepl(").+", commas_right))) || #Check so parentheses are not followed by anything
       isTRUE(any(grepl(".+\\(", commas_left))) || isTRUE(any(grepl(".+\\(", commas_right)))) { #chekc so parentheses are not preceded by anything
      stop("Incorrect hypothesis syntax or extra character, check specification")
    }


    framed$left <- gsub("[()]", "", framed$left) #drop remaining parentheses
    framed$right <- gsub("[()]", "", framed$right)
    commas <- unique(c(grep(",", framed$left), grep(",", framed$right))) #Gives us the unique rows that still contain commas (multiple comparisons) from left or right columns

    if(length(commas) > 0){ #If there are any multiple comparisons e.g., (X1, X2) below loop separates these in
      multiples <- vector("list", length = length(commas)) #Empty vector to store results for each row in loop below

      for(r in seq_along(commas)){ #for each row containing commas
        several <- framed[commas,][r, ] #select row r

        if(several$comp == "="){ #If e.g., (X1, X2) = X3, convert to X1 = X2 = X3

          several <- c(several$left, several$right)
          separate <- unlist(strsplit(several, split = ",")) #split by special characters and unlist
          if(any(grepl("^$", several))) stop("Misplaced comma in hypothesis") #if empty element
          hyp2 <- paste(separate, collapse = "=") #convert to X1 = X2 = X3 shape

          pos_comparisons <- unlist(gregexpr("=", hyp2)) #Gives the positions of all comparison signs
          leftside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
          rightside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
          pos1 <- c(-1, pos_comparisons) #positions to extract data to the leftside of comparisons
          pos2 <- c(pos_comparisons, nchar(hyp2) + 1) #positions to extract data to the rightside of comparisons
          for(i in seq_along(pos1)){
            leftside[i] <- substring(hyp2, pos1[i] + 1, pos1[i+1] - 1) #Extract all variables or outcomes to the leftside of a comparison sign
            rightside[i] <- substring(hyp2, pos2[i] + 1, pos2[i+1] - 1) #Extract all variables or outcomes to the rightside of a comparison sign
          }
          leftside <- leftside[-length(leftside)] #remove last element which is a NA due to loop formatting
          rightside <- rightside[-length(rightside)] #remove last element which is a NA due to loop formatting
          comp <- substring(hyp2, pos_comparisons, pos_comparisons) #Extract comparison signs
          multiples[[r]] <- data.frame(left = leftside, comp = comp, right = rightside, stringsAsFactors = FALSE) #hypotheses as a dataframe

        } else{ #If inequality comparison
          leftvars <- unlist(strsplit(several$left, split = ",")) #separate left hand var
          rightvars <- unlist(strsplit(several$right, split = ",")) #separate right hand vars
          if(any(grepl("^$", leftvars)) || any(grepl("^$", rightvars))) stop("Misplaced comma in hypothesis") #if empty element

          left <- rep(leftvars, length.out = length(rightvars)*length(leftvars)) #repeat each leftvars the number of rightvars
          right <- rep(rightvars, each = length(leftvars)) #complement for rightvars
          comp <- rep(several$comp, length(left)) #repeat the comparison a corresponding number of times

          multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE) #save as df to be able to combine with 'framed'
        }
      }

      framed <- framed[-commas,] #remove old unfixed rows with commas
      multiples <- do.call(rbind, multiples) #make list into dataframe
      framed <- rbind(multiples, framed) #recombine into one dataframe
    }
  } #end comma loop

  inequality <- framed[!framed$comp == "=",]

  #********Inequality
  if(nrow(inequality) == 0) { #If there are no '>' or '<' comparisons set to NULL
    list_inequality <- NULL
  } else{
    outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE) #Conversion to matrix in case there was only one row in outcomes
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
    specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
    specified <- specified[!is.na(specified)] #extract specified comparison values
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
    r_i <- matrix(r_i) #convert to matrix

    leq <- which(inequality$comp == "<") #gives the rows that contain '<' (lesser or equal) comparisons
    var_locations <- apply(inequality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0)) #convert non-variables to 0 and others are given their locations
    var_locations <- matrix(var_locations, ncol = 2) #Necessary if only one comparison row

    R_i <- matrix(rep(0, nrow(inequality)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix

    for(i in seq_along(r_i)){ # for each row i in R_i, replace the columns specified in var_locations row i
      if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)

        value <- if(i %in% leq) -1 else 1 #If comparison is 'lesser or equal' set to -1, if 'larger or equal' set to 1
        R_i[i, var_locations[i,]] <- value #Set this variable to 1 in R_i row i

      } else{ #If two variables specified
        value <- if(i %in% leq) c(-1, 1) else c(1, -1) #If comparison is 'leq' take var2 - var1, if 'larger or equal' take var1 - var2
        R_i[i, var_locations[i,]] <- value #Set one column to 1 and the other to -1 in R_i row i
      }
    }

    list_inequality<- list(R = R_i, r = r_i) #Note column 1 in R_i is for intercept
  }

  matrices <- list_inequality; matrices #final output

}


