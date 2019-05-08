###########################
#This is a function to calculate
#the sample size needed to produce
#the desired power to test whether
#a coefficient of a logistic
#regression is 0.

#The method used is a binary search
#tree that calculates the power for 
#a lower and upper bound of n using
#simulation, and then splits the 
#interval in half, iterating until
#the desired power is found.

#Author: Scott Koeneman
###########################



#First, define true model. Must specify all coefficients, including intercept
truebetas <- c(.5,-.5,-.5,.5)

#We must then specify our upper and lower limits for sample size,
#as well as our desired power and alpha, and our predictor of interest
nupper <- 1000
nlower <- 50
des_power <- .90
alpha <- .05
predinterest <- 1



#Now, we specify our predictors and assumptions for each
predictors <- c("pred1","pred2","pred3")
vartype <- c("Cont","Cont","Cont")
predmean <- c(0,0,0)
predvar <- c(1,1,1)
Predictors <- data.frame(predictors,vartype,predmean,predvar)
print(Predictors)


#We now begin the simulation. We set a tolerance for finding the power
tolerance <- .01


#Create dummy programming indicators to save time
iuppower <- 1
ilowpower <- 1

uppower <- 100
lowpower <- 0

###We first find the initial power for the upper n
results <- numeric(10000)
for (k in 1:10000){
  #We generate some combinations of predictors
  pred1sample <- rnorm(nupper,mean=Predictors$predmean[1],sd=sqrt(Predictors$predvar[1]))
  pred2sample <- rnorm(nupper,mean=Predictors$predmean[2],sd=sqrt(Predictors$predvar[2]))
  pred3sample <- rnorm(nupper,mean=Predictors$predmean[3],sd=sqrt(Predictors$predvar[3]))
  #We simulate outcomes of these, by finding the probability of success
  #given by our model for each combination of predictors and then performing
  #a bernoulli trial
  logits <- truebetas[1] + truebetas[2]*pred1sample +truebetas[3]*pred2sample +truebetas[4]*pred3sample
  pis <- (exp(logits))/(1+exp(logits))
  outcomes <-numeric(length(pis))
  for (j in 1:length(outcomes)){
    outcomes[j] <- rbinom(1,1,prob=pis[j])
  }
  #We can now form our model from the simulated data
  the_model<- glm(outcomes ~ pred1sample + pred2sample + pred3sample, family = binomial)
  #We call summary and see the result of the Wald test
  the_summary <- summary(the_model)
  p_value <- the_summary$coefficients[[2,4]]
  #We store the result of the test
  results[k] <- (p_value<alpha)
}
uppower <- mean(results)
###

###We repeat the process to find the power for lower n

results <- numeric(10000)
for (k in 1:10000){
  #We generate some combinations of predictors
  pred1sample <- rnorm(nlower,mean=Predictors$predmean[1],sd=sqrt(Predictors$predvar[1]))
  pred2sample <- rnorm(nlower,mean=Predictors$predmean[2],sd=sqrt(Predictors$predvar[2]))
  pred3sample <- rnorm(nlower,mean=Predictors$predmean[3],sd=sqrt(Predictors$predvar[3]))
  #We simulate outcomes of these, by finding the probability of success
  #given by our model for each combination of predictors and then performing
  #a bernoulli trial
  logits <- truebetas[1] + truebetas[2]*pred1sample +truebetas[3]*pred2sample +truebetas[4]*pred3sample
  pis <- (exp(logits))/(1+exp(logits))
  outcomes <-numeric(length(pis))
  for (j in 1:length(outcomes)){
    outcomes[j] <- rbinom(1,1,prob=pis[j])
  }
  #We can now form our model from the simulated data
  the_model<- glm(outcomes ~ pred1sample + pred2sample + pred3sample, family = binomial)
  #We call summary and see the result of the Wald test
  the_summary <- summary(the_model)
  p_value <- the_summary$coefficients[[2,4]]
  #We store the result of the test
  results[k] <- (p_value<alpha)
}
lowpower <- mean(results)

#We see which power is closer to have a "working estimate" for n,
#and if the upper and lower limits contain the desired power
if(uppower < des_power){
  stop("Upper bound is too low")
}
if(lowpower > des_power){
  stop("Lower bound is too high")
}
if (abs(uppower-des_power)>abs(lowpower-des_power)){
  the_power <- lowpower
  n_est <- nlower
} else {
  the_power <- uppower
  n_est <- nupper
}

#We now start the binary search to converge to a sample size
#We find the correct sample s
while(abs(the_power-des_power)>=tolerance){
    #Find the midpoint n
    nmid <- round((nupper+nlower)/2)
    
    #We now find the power for this new n
    results <- numeric(10000)
    #Perform the test many times
    for (k in 1:10000){
      #We generate some combinations of predictors
      pred1sample <- rnorm(nmid,mean=Predictors$predmean[1],sd=sqrt(Predictors$predvar[1]))
      pred2sample <- rnorm(nmid,mean=Predictors$predmean[2],sd=sqrt(Predictors$predvar[2]))
      pred3sample <- rnorm(nmid,mean=Predictors$predmean[3],sd=sqrt(Predictors$predvar[3]))
      #We simulate outcomes of these, by finding the probability of success
      #given by our model for each combination of predictors and then performing
      #a bernoulli trial
      logits <- truebetas[1] + truebetas[2]*pred1sample +truebetas[3]*pred2sample +truebetas[4]*pred3sample
      pis <- (exp(logits))/(1+exp(logits))
      outcomes <-numeric(length(pis))
      for (j in 1:length(outcomes)){
        outcomes[j] <- rbinom(1,1,prob=pis[j])
      }
      #We can now form our model from the simulated data
      the_model<- glm(outcomes ~ pred1sample + pred2sample + pred3sample, family = binomial)
      #We call summary and see the result of the Wald test
      the_summary <- summary(the_model)
      p_value <- the_summary$coefficients[[2,4]]
      #We store the result of the test
      results[k] <- (p_value<alpha)
    }
    #We now find the power of the midpoint
    midpower <- mean(results)
    
    #We now select which region has our desired power, and iterate
    if (midpower < des_power){
      nlower <- nmid
      lowpower <- midpower
    } else{
      nupper <- nmid
      uppower <- midpower
    }
    
    #Find best estimate for n again
    if (abs(uppower-des_power)>abs(lowpower-des_power)){
      the_power <- lowpower
      n_est <- nlower
    } else {
      the_power <- uppower
      n_est <- nupper
    }
  
}
#After the loop, we should be converged to about our desired power and
#have the desired sample size
print(n_est)
print(the_power)




