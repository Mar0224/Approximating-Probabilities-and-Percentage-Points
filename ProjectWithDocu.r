# STAT 181 PROJECT - Group 2
# Topic: Approximating Probabilities and Percentage Points
# Description: An R package BLAH BLAH

# APPROXIMATING NORMAL DISTRIBUTION
# A. CUMULATIVE DISTRIBUTION FUNCTION
# Simpson's Rule
# parameters: x-> value at which the function is evaluated; mean->mean; sd-> standard deviation
ndsimpsons<-function(x, mean, sd){
  z <- (x-mean)/sd
  #function to be integrated
  function0<-function(t)(exp(-t^2/2))                    # where z > 0
  intFunction<-function(z)(integrate(function0,0,z)$val) # code for integration
  #function for simpson's rule
  Fz<-0.5+(1/sqrt(2*pi))*intFunction(z)
  return(Fz)
}

# Closed Form - Hastings
# parameters: x-> value at which the function is evaluated; mean->mean; sd-> standard deviation
ndhast<-function(x, mean, sd){
  z <- (x-mean)/sd
  #constant values
  a<-c(0.196854,0.115194,0.000344,0.019527)
  #solve for sum
  sum<-(1+((a[1]*z)+(a[2]*(z^2))+(a[3]*(z^3))+(a[4]*(z^4))))
  sum1 = 1/(sum^4)
  #solve for probability
  Fz<-1-(0.5*sum1) # plus error but we approximate it to 0 since the value is very small
  #print probability
  return(Fz)
}

# Closed Form - Abromowitz and Stegun
# parameters: x-> value at which the function is evaluated; mean->mean; sd-> standard deviation
ndabstegun<-function(x, mean, sd){
  z <- (x-mean)/sd
  #constant values
  b <- c(0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
  #Solve for t
  t <- 1/(1 + (0.2316419 * abs(z)))
  
  f2 <- function(z){                       #computes probability of z
    (exp(-0.5 * (z^2)))/sqrt(2*pi)
  }
  y <- 0
  for (i in 1:5) {
    y <- y + b[i] * (t^i)
  }
  ans <- 1-(f2(z) * y)
  
  #Decision Rule based on the sign of z
  if (z >= 0) { return(ans) }
  else{ return(1-ans)}
}

# Continued Fraction
# parameters: x-> value at which the function is evaluated; mean->mean; sd-> standard deviation; K->number of recurrence
ndcfraction<-function(x,mean,sd, K){
  #the parameter K is for the value of k during forward recurrence
  z<-(x-mean)/sd
  Fz<-0 
  
  a<-c() #initialize vector of a_k
  for(i in 1:K){ #perform the fibonacci sequence to get the values of a_k
    previous<-0 #lower value
    current<-1 #higher value
    
    a[i]<-current #value of a[i] is the higher value
    
    current<-previous+current #add the two values
    previous<-a[i] #change the lower value to the previous higher value
  }
  A<-c(1,0) #A(-1) and A(0)
  B<-c(0,1) #B(-1) and B(0)
  V<-0 #initialize V as 0
  k<-1 #initialize k as 1
  for(i in 3:(K+2)){ #start i in 3 to access -1 and 0 of A and B
    A[i]<-(z*A[i-1])+(a[k]*A[i-2]) #formula for A[i]
    B[i]<-(z*B[i-1])+(a[k]*B[i-2]) #formula for B[i]
    V[k]<-A[i]/B[i] #divide A[i] by B[i] to get V[k]
    k<-k+1
  }
  A<-(exp((-z^2)/2)*V[K])
  
  if(z>0){ #if z is greater than 0
    Fz<-1-((1/sqrt(2*pi))*A)      #normal distribution using ERF
  }
  else if(z<0){ #if z is less than 0
    Fz<-abs((1/sqrt(2*pi))*A)     #normal distribution using ERFC
  }
  else{ #if z is exactly 0
    Fz<-0.5 #assign 0.5
  }
  return(Fz) #return the probability
}

# ERF and ERFC Function (Error Function)
# parameters: x-> value at which the function is evaluated; mean->mean; sd-> standard deviation
nderfunction<-function(x,mean,sd){
  z<-(x-mean)/sd 
  Fz<-0
  
  ERF<-function(x){
    2*pnorm(x*sqrt(2))-1 #from documentation of normal distribution
  } 
  ERFC<-function(x){
    2*pnorm(x*sqrt(2),lower=FALSE) #from documentation of normal distribution
  } 
  
  if(z>=0){ #if z is greater than or equal to 0
    Fz<-(1+ERF(z/sqrt(2)))/2 #normal distribution using ERF
  }
  else{ #if z is less than 0
    Fz<-ERFC(-1*z/sqrt(2))/2 #normal distribution using ERFC
  }
  return(Fz) #return the probability
}

# B. PERCENTAGE POINTS
# General Series Expansion (Hastings)
# parameter: p-> probability of the unknown standard normal variate (z)
npGSE<-function(p){ 
  p_prime<-0 #declare p'
  c<-c(2.515517,0.802853,0.010328) #create a vector for values of c
  d<-c(1.432788,0.189269,0.001308) #create a vector for values of d
  numerator<-0 #declare numerator to compute for z
  denominator<-0 #declare denominator to compute for z
  
  if(p>0.5){ #if the given p is greater than 0.5
    p_prime<-(1-p) #get 1-p for p'
  }
  else{ #if the given p is lesser than or equal to 0.5
    p_prime<-p #get p for p'
  }
  w<-sqrt(-2*log(p_prime)) #get the value of w
  
  for(i in 1:length(c)){
    numerator<-numerator+(c[i]*w^(i-1)) #calculate the numerator
  }
  for(i in 1:length(d)){
    denominator<-denominator+(d[i]*w^i) #calculate the denominator
  }
  denominator<-denominator+1 #add 1 to the denominator
  
  z=w-(numerator/denominator) #get the value of the z
  return(z) #return z
}


# APPROXIMATING STUDENT T-DISTRIBUTION
# A. DISTRIBUTION FUNCTION
# Wallace Approximate Transformation
# parameters: x-> value at which the function is evaluated; df->degrees of freedom
tdwallace<-function(x,df){
  z<-((8*df+1)/(8*df+3))*sqrt(df*log(1+(x^2/df)))
  Fz <- pnorm(z)
  return(Fz)
}

# B. PERCENTAGE POINTS
# Wallace approximate transformation by a random variable Z
# parameters: alpha-> alpha; df->degrees of freedom
tpwallace <- function(alpha,df){
  z <- qnorm(1-alpha)
  t <- sqrt(df*(exp((z^2/df)*((8*df+3)/(8*df+1))^2)-1))
  return(t)
}


# APPROXIMATING CHI-SQUARE DISTRIBUTION
# A. DISTRIBUTION FUNCTION
# Wilson and Hilferty approximate transformation
# parameters: x-> value at which the function is evaluated; df->degrees of freedom
cdwilhil <- function(x,df){
  z <- ((((x/df)^(1/3))-(1-(2/(9*df))))/(sqrt(2/(9*df)))) #compute for the value of z
  Fz <- pnorm(z)
  return(Fz)
}

# B. PERCENTAGE POINTS
# Wallace approximate transformation by a random variable Z
# parameters: alpha-> alpha; df->degrees of freedom
cpwilhil <- function(alpha,df){
  z<-qnorm(1-alpha)
  chi<-df*(((z*sqrt(2/(9*df))) + (1-(2/(9*df))))^3)
  return(chi)
}


# APPROXIMATING F-DISTRIBUTION
# A. DISTRIBUTION FUNCTION
# Simpson's Rule
# parameters:x-> value at which the function is evaluated; df1->numerator's degree of freedom; df2->denominator's degree of freedom
fdsimpsons <- function(x,df1,df2){
  f <- function(t){
    ((1+(df1/df2*t))^(-(df1+df2)/2))*t^((df1/2)-1)
  }
  numerator <- gamma((df1+df2)/2)*(sqrt((df1/df2)^df1))
  denominator <- gamma((df1/2))*gamma((df2/2))
  int <- (integrate(f,0,x)$val)
  Fx <- (numerator/denominator)*int
  return(Fx)
}

# Approximate Transformation
# parameters: x-> value at which the function is evaluated; df1->numerator's degree of freedom; df2->denominator's degree of freedom
fdapproxtrans <- function(x,df1,df2){
  numerator <- ((1-(2/(9*df2)))*x^(1/3))-(1-(2/(9*df1)))
  denominator <- sqrt(((2/(9*df2))*x^(2/3))+(2/(9*df1)))
  z <- numerator/denominator
  Fx <- pnorm(z)
  return(Fx)
}

# Closed Form Approximation
# parameters: x-> value at which the function is evaluated; df1->numerator's degree of freedom; df2->denominator's degree of freedom
fdclosed <- function(x,df1,df2){ 
  a <- c(0.278393,0.230389,0.000972,0.078108)
  
  # formula of u:
  u <- (((1-(2/(9*df2)))*(x^(1/3))) - (1-(2/(9*df1)))) / (sqrt(((2/(9*df2))*(x^(2/3)))+(2/(9*df1))))
  # formula of y:
  y <- round(u/(sqrt(2)), digits = 4)
  # formula of F(x):
  Fx <- 1-((1/2)*(1+(a[1]*y)+(a[2]*(y^2))+(a[3]*(y^3))+(a[4]*(y^4)))^(-4))
  
  return(Fx)
}

# B. PERCENTAGE POINTS
# Wallace approximate transformation by a random variable Z
# parameters: alpha-> alpha; df1->numerator's degree of freedom; df2->denominator's degree of freedom
fppoint <- function(alpha,df1,df2){ # function in computing for percentage point in F distribution
  z<-qnorm(1-alpha)
  a <- ((1-(2/(9*df2)))^2) - ((2*(z^2))/(9*df2))
  b <- -(1-(2/(9*df1)))*(1-(2/(9*df2)))
  c <- ((1-(2/(9*df1)))^2) - ((2*(z^2))/(9*df1))
  
  f <- round(((-b+(sqrt((b^2)-(a*c))))/a)^3, digits = 2) # to compute for percentage point
  
  return(f)
}