\name{LTRCdata}
\alias{LTRCdata}
\docType{data}
\title{Simulated Left Truncated, Right Censored Data Set for a Irreversible Illness-Death Model}
\description{
Simulated data from an irreversible three-state illness-death  model with data subject to
independent right censoring and left truncation. 
}
\usage{data(LTRCdata)}
\format{
  A data frame with 1000 individuals with the following 4 variables.
  \describe{
    \item{\code{id}}{Identification number}
    \item{\code{start}}{Left truncation time, i.e.: start time for the period of observation after the individual enters state j }
    \item{\code{stop}}{Transition time}
    \item{\code{st.stage}}{State transitioning FROM}
    \item{\code{stage}}{State transitioning TO}
  }
}

\details{
We consider an irreversible three-state illness-death  model with data subject to
independent right censoring and left truncation.  A data
set of 1000 individuals  starting in state 1 at time 0 was simulated.
Individuals remain in state 1 until they
transition to the transient state 2 (ill) or the terminal state 3
(death).  Individuals in state 2 remain there until they transition to
the terminal state 3 (death).  All times were rounded to the fourth
decimal place for clarity of presentation.

For further details concerning the simulation, please see the package vignette.
}

\examples{
data(LTRCdata)

####################################################################
##               Code used to generate data                       ##
####################################################################

## Specifying the tree structure for the simulation
Nodes <- c("1","2","3") #states possible in MSM
Edges <- list("1"=list(edges=c("2","3")),"2"=list(edges=c("3")),
           "3"=list(edges=NULL)) #transitions from each state

treeobj2 <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

## Simulating the data
set.seed(123)
n<-1000
censor <- round(rlnorm(n, meanlog=-0.5, sdlog=2),digits=4)

## Left truncation times are mixture of 80% lognormal, 20% LT time=0
p1 <- n*0.8
p2 <- n*0.2
left.p1 <- round(rlnorm(p1, meanlog=-0.5, sdlog=2),digits=4) #start time
left.p2 <- rep(0,p2) #start time
left<-c(left.p1,left.p2)
ill <- round(rweibull(n,2),digits=4)
death <- round(rweibull(n,2),digits=4)

max.tran <- pmax(death,censor)
hidden <- which(max.tran < left | death < left | censor < left) #those who are never visible to investigator

present <- which(left < censor & left < death & left<ill)  #present in state 1
illstart <- which(ill < left & left < censor & left < death) #present in study but starts in "ill"

###All present @ any time is which(left < censor & left < death)

x <- cbind(left,death,censor,ill)

#those who enter in stage 1
c1 <-  pmin(ill[present],death[present],censor[present])
num.ill <- sum(ill[present]<death[present]  & ill[present]<censor[present]) #number who are "ill"
stage <- ifelse(censor[present]<death[present] & censor[present]<ill[present],
                0,ifelse(ill[present]<death[present],2,3))
data1 <- data.frame(id=1:length(c1), start=left[present], stop=c1,st.stage=1, stage=stage)


ind <- which(data1$stage==2) #those transitioning to stage 2
s.cens <- censor[present]
n2 <- sum(as.numeric(data1$stage==2)) #number that made transition to 2
death2 <- round(qweibull(pweibull(ill,shape=2)+runif(n,0,1)*(1-pweibull(ill,shape=2)),shape=2),digits=4)
c2 <- pmin(death2[ind],s.cens[ind])
stage[ind] <- ifelse(censor[ind]<c2,0,3)
data2 <- data.frame(id=ind, start=data1$stop[ind], stop=c2, st.stage=data1$stage[ind], stage=stage[ind])

if(length(illstart)>0){
#those who enter in stage 2 after truncation
	c2a <-  pmin(death2[illstart],censor[illstart])
	stage[illstart] <- ifelse(censor[illstart]<death2[illstart],0,3)
	data2a <- data.frame(id=(1:length(c2a))+max(data1$id), start=left[illstart], stop=c2a,st.stage=2, stage=stage[illstart])
	data <- rbind(data1,data2,data2a)}

if(length(illstart)==0)	data <- rbind(data1,data2)

LTRCdata<-with(data, data[order(stop), ])

}


\keyword{datasets}
