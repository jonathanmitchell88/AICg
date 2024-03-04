# Set working directory.

library(MSCquartets)

#####################
##### Functions #####
#####################

# Read all functions in.

# Error function

erf <- function(x) 2*pnorm(x*sqrt(2))-1

# Function to integrate for the bias of the MLE of mu0 for model T3.

integrand <- function(y,mu0,alpha0,beta0) {
  (1/sqrt(2*pi))*y*(exp(-(1/2)*(y-mu0)^2)*erf(y*(1/tan(beta0))/sqrt(2))+exp(-(1/2)*(y+mu0*sin(alpha0))^2)*(erf((y*(1/tan(alpha0))+mu0*cos(alpha0))/sqrt(2))+erf((y*tan(alpha0+beta0)-mu0*cos(alpha0))/sqrt(2))))
}

# Bias of the MLE of mu0 for model T3.

T3bias <- function(mu02,mu0,n) {
  phi02 <- min((1/(4*(n+mu02^2)))*(4*n+3*mu02^2-mu02*sqrt(8*n+9*mu02^2)),1)
  alpha02 <- atan(1/sqrt(3*(3-2*phi02)))
  beta02 <- (1/2)*(pi/2-alpha02)
  bias <- integrate(integrand,lower=0,upper=Inf,mu0=mu02,alpha0=alpha02,beta0=beta02)$value-mu0
  return(bias)
}

# Function to integrate for T3 bias correction.

integrand1 <- function(y,mu0,beta0) (y-mu0)^2*exp(-(1/2)*(y-mu0)^2)*erf(y/(sqrt(2)*tan(beta0)))

# Function to integrate for T3 bias correction.

integrand2 <- function(phi,mu0,alpha0) 8*exp(-(1/2)*mu0^2)*(mu0^2*cos(phi)*(cos(phi)+sin(phi)*cos(alpha0+phi)*sin(alpha0+phi))+2*cos(alpha0+phi)^2)-sqrt(2*pi)*mu0*exp(-(1/2)*mu0^2*cos(phi)^2)*(1+erf(mu0*sin(phi)/sqrt(2)))*(-mu0^2*sin(2*phi)*(2*cos(phi)+sin(phi)*sin(2*(alpha0+phi)))+sin(2*alpha0+phi)-3*sin(2*alpha0+3*phi))

# T3 bias correction

T3biascorrection <- function(mu0,alpha0,beta0) sqrt(2/pi)*integrate(integrand1,0,Inf,mu0=mu0,beta0=beta0)$value+1/(4*pi)*integrate(integrand2,-pi/2,beta0,mu0=mu0,alpha0=alpha0)$value

# Star tree penalty, equal to the T3 bias correction minus some arbitrary quantity, diff.

starpenalty <- function(biascorrT3,diff) biascorrT3-diff

# Function to append AICg weights to a resolved quartet table.

quartetTableResolvedAICg <- function(RQT,alpha,diff) {
  p_T3 <- vector(length=nrow(RQT))
  p_star <- vector(length=nrow(RQT))
  
  for (i in 1:nrow(RQT)) {
    obs <- as.numeric(RQT[i,(ncol(RQT)-2):ncol(RQT)])
    n <- sum(obs)
    obs <- sort(obs,decreasing=TRUE)
    expd <- c(obs[1],(1/2)*(n-obs[1]),(1/2)*(n-obs[1]))
    
    phi0 <- 3*(n-expd[1])/(2*n)
    mu0 <- max(0,sqrt(2*n)*(1-phi0)/sqrt(phi0*(3-2*phi0)))
    
    if (mu0==Inf) {
      biascorrT3 <- 2
    } else {
      if (mu0 <= 3*sqrt(3)/(2*sqrt(2*pi))) {
        mu0unbiased <- 0
      } else if (sign(T3bias(0,mu0,n))==sign(T3bias(mu0,mu0,n))) {
        mu0unbiased <- mu0
      } else {
        mu0unbiased <- uniroot(T3bias,c(0,mu0),mu0=mu0,n=n)$root
      }
      phi0unbiased <- (1/(4*(n+mu0unbiased^2)))*(4*n+3*mu0unbiased^2-mu0unbiased*sqrt(8*n+9*mu0unbiased^2))
      alpha0unbiased <- atan(1/sqrt(3*(3-2*phi0unbiased)))
      beta0unbiased <- (1/2)*(pi/2-alpha0unbiased)
      biascorrT3 <- T3biascorrection(mu0unbiased,alpha0unbiased,beta0unbiased)
    }
    
    starp <- starpenalty(biascorrT3,diff)
    
    AICgstar <- -2*dmultinom(obs,prob=rep(1/3,3),log=TRUE)+starp
    AICgT3 <- -2*dmultinom(obs,prob=expd/n,log=TRUE)+biascorrT3
    AICgnetwork <- -2*dmultinom(obs,prob=obs/n,log=TRUE)+4
    
    minAICg <- min(AICgstar,AICgT3,AICgnetwork)
    
    AICgstar <- AICgstar-minAICg
    AICgT3 <- AICgT3-minAICg
    AICgnetwork <- AICgnetwork-minAICg
    
    p_star[i] <- ifelse(AICgstar==min(c(AICgstar,AICgT3,AICgnetwork)),exp(-1/2*AICgstar)/(exp(-1/2*AICgstar)+exp(-1/2*AICgT3)+exp(-1/2*AICgfull)),0)
    p_T3[i] <- ifelse(AICgstar==min(c(AICgstar,AICgT3,AICgnetwork)),NA,exp(-1/2*AICgT3)/(exp(-1/2*AICgT3)+AICgnetwork))
  }
  
  tests <- length(p_T3[!is.na(p_T3)])
  
  pTable <- cbind(cbind(RQT,p_T3),p_star)
  
  ps <- pTable[,(ncol(pTable)-1)]
  ordering <- order(ps)
  psordered <- ps[ordering]
  
  for (i in 1:tests) {
    psordered[i] <- min(psordered[i]*(tests+1-i),1)
  }
  
  psordered[which(psordered>alpha)[1]:tests] <- 1
  ps[ordering] <- psordered
  ps[is.na(ps)] <- 1
  
  pTable[,(ncol(pTable)-1)] <- ps
  
  return(pTable)
}

####################
##### Analysis #####
####################

# Parameters to be set by the user.

alpha <- 0.01 # Significance level for Holm-Bonferroni adjusted Akaike weights.
beta <- 1/3 # Select star tree model when its AICg is no more than that of model T3 and the network model and its AICg weight exceeds this value.
diff <- 0.1 # Equals the T3 bias correction minus the star tree penalty.

# Inference of a network.

gtrees <-  read.tree(file=system.file("extdata","dataYeastRokas",package="MSCquartets")) # Gene trees file in Newick format.
tnames <- taxonNames(gtrees)
QT <- quartetTable(gtrees)
RQT <- quartetTableResolved(QT)
AICgTable <- quartetTableResolvedAICg(RQT,alpha,diff)

NANUQ(AICgTable,alpha=alpha,beta=beta,outfile=file.path(tempdir(),"NANUQdist"))
dist=NANUQdist(AICgTable,alpha=alpha,beta=beta,outfile=file.path(tempdir(),"NANUQdist"))

nn=neighborNet(dist) # Inferred network
plot(nn,"2D")
