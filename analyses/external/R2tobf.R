R2.to.bf <- function (n, nPred, r2) {
  # Computes Bayes factor of Models M to the Null Model and to all other models  
  # from R^2 obtained from hierarchical regression analysis
  # The output can be provided for all models conducted for the same sample                
  # Last update: 04-08-2014                                    
  # Args:
  #   n: Sample size
  #   nPred: an array containing mount of predictors
  #   r2: an array containing R^2 values
  # Returns:
  #   output: a matrix containing combinations of all Bayes factors
  #           first two columns correspond to the models, e.g. 2 1
  #           is a Bayes factor of model 2 to model 1
  # Additionally, the function displays the output in a descriptive form
  # in the console
  # Example:
  #   Hierarchical regression analysis was conduced for 500 participants,
  #   where the predictors are entered in 3 steps:
  #   Step 1) 1 predictors, R^2=.38
  #   Step 2) 2 predictors, R^2=.56
  #   Step 3) 10 predictors, R^2=.70
  #   To obtain Bayes factors type in the console:
  #   R2.to.bf(198, c(1, 2, 10), c(.38, .56, .70))
  
  #Load neccessary packages
  library(BayesFactor)
  library(caTools)
  library(coda)
  library(lattice)
  library(mvtnorm)
  library(pbapply)
  library(bitops)
  #source("linReg.R")
  
  nModels <- length(r2) #Aomunt of models
  
  #Initialise output Bayes Factors comparing Model M to the Null Model
  bfm0 <- c(rep(0, nModels))
  
  #Obtain Bayes Factors comparing Models M to the Null Model
  for(i in seq(1:nModels)) {
    temp <- linearReg.R2stat(n, nPred[i], r2[i],  rscale = 1)
    bfm0[i] <- exp(temp[1]$bf)
  }
  
  #Obtain Bayes factors comparing Models M with each other
  myCombs <- combs(seq(1:nModels), 2) #Create all combinations of BF
  bfmm <- c(rep(0, nrow(myCombs))) #Initialise vector containing BF among other models
  bfmmnames <- c(rep("", nrow(myCombs))) #Initialise vector containing names of BFs
  for(i in seq(1:nrow(myCombs))){
    bfmm[i] <- bfm0[myCombs[i, 2]] / bfm0[myCombs[i, 1]]
    bfmmnames[i] <- paste("BF", paste(as.character(myCombs[i, 2]), as.character(myCombs[i, 1]), collapse=NULL), collapse=NULL)
  }
  
  output <- cbind(myCombs, bfmm)
  return(output)
}