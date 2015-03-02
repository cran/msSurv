\name{LTRCdata}
\alias{LTRCdata}
\docType{data}
\title{Simulated Left Truncated, Right Censored Data Set for an
  Irreversible Illness-Death Model}
\description{
Simulated data from an irreversible three-state illness-death  model
with data subject to independent right censoring and left truncation.
}
\usage{data(LTRCdata)}
\format{
  A data frame with 1000 individuals with the following 4 variables.
  \describe{
    \item{\code{id}}{Identification number}
    \item{\code{start}}{Left truncation time, i.e.: start time for the
      period of observation after the individual enters state j }
    \item{\code{stop}}{Transition time}
    \item{\code{start.stage}}{State transitioning FROM}
    \item{\code{end.stage}}{State transitioning TO}
  }
}

\details{
A data set of 1000 individuals was simulated from an irreversible
three-state illness-death  model with data subject to
independent right censoring and left truncation.  All individuals were
assumed to starte in state 1 at time 0.
Individuals remained in state 1 until they
transitioned to the transient state 2 (ill) or the terminal state 3
(death).  Individuals in state 2 remained there until they transition to
the terminal state 3 (death).  All times were rounded to the fourth
decimal place for clarity of presentation.
}

\references{
  Nicole Ferguson, Somnath Datta, Guy Brock (2012).
  msSurv: An R Package for Nonparametric Estimation of Multistate Models.
  Journal of Statistical Software, 50(14), 1-24.
  URL http://www.jstatsoft.org/v50/i14/.
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

RCLTtree <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

## Simulating the data
set.seed(123)
n <- 1000
censor <- round(rlnorm(n, meanlog=0, sdlog=2),digits=4)

## 80% of sample is LT, rest has start time of 0
p1 <- n*0.8
left.p1 <- round(rlnorm(p1, meanlog=-1, sdlog=2),digits=4) ## start time
left.p2 <- rep(0, n-p1) ## start time
left <- c(left.p1, left.p2)
ill <- round(rweibull(n,2),digits=4)
death1 <- round(rweibull(n,2),digits=4)
## 2nd transition time for those entering illness state first
death2 <- round(qweibull(pweibull(ill, shape=2) + runif(n, 0, 1)*(1 -
                         pweibull(ill, shape=2)), shape=2), digits=4)
## use death2 for indiv w ill < death1, death1 for those w/death1 < ill
death <- ifelse(ill < death1, death2, death1)
last.tran <- pmin(death, censor)

## those who are never visible to investigator
hidden <- which(left > last.tran)
length(hidden) ## 407

## indiv starting in state 2 (illness)
ill.start <- which(left > ill & left < last.tran & ill < last.tran)
length(ill.start)  ## 28

## remainder should be indiv starting in state 1 (wellness)
## can double-check w/the following ..
first.tran <- pmin(death, pmin(ill, censor))
well.start <- which(left < first.tran)
length(c(hidden, well.start, ill.start))    ## 1000

x <- cbind(left, death, censor, ill)

## those who enter in state 1
t1 <-  pmin(ill[well.start], death[well.start], censor[well.start])
end.stage <- ifelse(censor[well.start] < death[well.start] &
                        censor[well.start] < ill[well.start], 0,
                    ifelse(ill[well.start] < death[well.start], 2, 3))
data1 <- data.frame(id=1:length(t1), start=left[well.start], stop=t1,
                    start.stage=1, end.stage=end.stage)

## those transitioning to stage 2
ind <- which(data1$end.stage==2)
t2 <- pmin(death[well.start][ind], censor[well.start][ind])
end.stage2 <- ifelse(censor[well.start][ind] < t2, 0, 3)
data2 <- data.frame(id=ind, start=data1$stop[ind], stop=t2,
                    start.stage=data1$end.stage[ind], end.stage=end.stage2)

## those who enter in stage 2 after truncation
if(length(ill.start) > 0) {
  t3 <-  pmin(death[ill.start], censor[ill.start]) ## same as last.tran[ill.start]
  end.stage3 <- ifelse(censor[ill.start] < death[ill.start], 0, 3)
  data2a <- data.frame(id=(1:length(t3)) + max(data1$id),
                       start=left[ill.start],
                       stop=t3, start.stage=2, end.stage=end.stage3)
  data <- rbind(data1, data2, data2a)
}

if(length(ill.start)==0)  data <- rbind(data1, data2)

LTRCdata <- with(data, data[order(id, stop), ])

}


\keyword{datasets}
