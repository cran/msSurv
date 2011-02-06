\name{RCdata}
\alias{RCdata}
\docType{data}
\title{Simulated Right Censored Data Set for a 5 State Model}
\description{
Simulated data from a five-state progressive model with both branching
and tracking between states.  Possible transitions are 1->2, 1->3, 2->4,
and 2->5, with  states 3, 4, and 5 as terminal states.  
}
\usage{data(RCdata)}
\format{
  A data frame with 1000 individuals with the following 4 variables.
  \describe{
    \item{\code{id}}{Identification number}
    \item{\code{stop}}{Transition time}
    \item{\code{st.stage}}{State transitioning FROM}
    \item{\code{stage}}{State transitioning TO}
  }
}

\details{
We simulated a data set of 1000 individuals subject to independent right
censoring with 60% of individuals starting in state 1 at time 0 and 40%
starting in state 2.  Those in state 1 remain there until they
transition to the transient state 2 or the terminal state 3.
Individuals in state 2 remain there until they transition to 
either terminal state 4 or 5.
}

\examples{
data(RCdata)

####################################################################
##               Code used to generate data                       ##
####################################################################

## Specifying the tree structure for the simulation

Nodes <- c("1","2","3","4","5") #states possible in MSM
Edges <- list("1"=list(edges=c("2","3")),"2"=list(edges=c("4","5")),
           "3"=list(edges=NULL),"4"=list(edges=NULL),"5"=list(edges=NULL))
treeobj <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

## Simulating the data
set.seed(123)
n <- 1000
n1 <- 0.6*n
n2<-0.4*n
ill <- round(rweibull(n1,2),digits=4)
death1 <- round(rweibull(n1,2),digits=4)
allcensor <- round(rlnorm(n, meanlog=-0.5, sdlog=2),digits=4)
censor<-allcensor[1:n1]
censor2<-allcensor[(1:n2)+length(n1)]

c1 <- pmin(ill,death1,censor)
num.ill <- sum(ill<death1  & ill<censor) #number who are "ill"
stage <- ifelse(censor<death1 & censor<ill, 0,ifelse(ill<death1,2,3))
data1 <- data.frame(id=1:length(c1), stop=c1,st.stage=1,stage=stage)

ind <- which(data1$stage==2) #those transitioning to state 2
#no.st2 <- sum(as.numeric(data1$stage==2)) #number that made transition to 2
death24 <- round(qweibull(pweibull(ill,shape=2)+runif(n1,0,1)*(1-pweibull(ill,shape=2)),shape=2),digits=4)
death25 <- round(qweibull(pweibull(ill,shape=2)+runif(n1,0,1)*(1-pweibull(ill,shape=2)),shape=2),digits=4)
c2 <- pmin(death24[ind],death25[ind],censor[ind])
stage[ind] <- ifelse(censor[ind]<death24[ind] & censor[ind]<death25[ind],0,ifelse(death24[ind]<death25[ind],4,5))
data2 <- data.frame(id=ind, stop=c2, st.stage=data1$stage[ind],stage=stage[ind])

death24.2 <- round(rweibull(n2,2),digits=4)
death25.2 <- round(rweibull(n2,2),digits=4)
tran<-pmin(death24.2,death25.2,censor2)
stage2 <- ifelse(censor2<death24.2 & censor2<death25.2,0,ifelse(death24.2<death25.2,4,5))
data2.2 <- data.frame(id=(1:n2)+length(c1), stop=tran,st.stage=2,stage=stage2)

RCdata <- rbind(data1,data2,data2.2)

}

\keyword{datasets}
