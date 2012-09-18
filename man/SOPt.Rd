\name{SOPt}
\alias{SOPt}

\title{
State Occupation Probability at Time t
}
\description{
  This function displays the state occupation probability at time t and
  then prints it.  The function also displays the corresponding
  variance if the user specifies.
}
\usage{
SOPt(object, t="last", deci=4, covar=FALSE)
}

\arguments{
  \item{object}{
    A msSurv object.
}
  \item{t}{
     The time at which we find the  state occupation probabilities.  Default is "last" which is the highest (or "last") event time.
}
  \item{deci}{
     Numeric argument specifying number of decimal places for estimates.  Default is 4.
}
  \item{covar}{
     Logical argument to determine if variance is reported.  Default is FALSE.
}
}
\details{

Displays the state occupation probabilities and corresponding variance for all states in the multistate model.  
   
}

\value{
  Returned invisibly:
  \describe{
    \item{SOPt}{The state occupation probability at time t.}
    \item{var.sop}{The variance of the state occupation probabilities at
      time t (if \code{covar==TRUE}).}
    }
}

\references{
  Nicole Ferguson, Somnath Datta, Guy Brock (2012).
  msSurv: An R Package for Nonparametric Estimation of Multistate Models.
  Journal of Statistical Software, 50(14), 1-24.
  URL http://www.jstatsoft.org/v50/i14/.

  Datta, S. and Satten G.A. (2001). Validity of the Aalen-Johansen
  estimators of stage occupation probabilities and Nelson-Aalen
  estimators of integrated transition hazards for non-Markov models. 
  Statistics and Probability Letters, 55(4): 403-411.
 
  Datta S, Satten GA (2002). Estimation of Integrated Transition Hazards
  and Stage Occupation Probabilities for Non-Markov Systems under
  Dependent Censoring. Biometrics, 58(4), 792-802.
}
\author{
   Nicole Ferguson <nicole.ferguson@kennesaw.edu>,
   Guy Brock <guy.brock@louisville.edu>,
   Somnath Datta <somnath.datta@louisville.edu>
}

\seealso{
   \code{\link{msSurv}}
}
\examples{

## Row elements of data 
p1 <- c(1,0,0.21,1,3)
p2 <- c(2,0,0.799,1,2)
p22 <- c(2,0.799,1.577,2,3)
p3 <- c(3,0,0.199,1,0)

## Combining data into a matrix
ex1 <- rbind(p1,p2,p22,p3)
colnames(ex1) <- c("id", "start", "stop", "start.stage", "end.stage")
ex1 <- data.frame(id=ex1[,1], start=ex1[,2], stop=ex1[,3],
                  start.stage=ex1[,4], end.stage=ex1[,5])


## Inputting nodes & edges of the tree structure
Nodes <- c("1","2","3") # states in MSM
Edges <- list("1"=list(edges=c("2","3")),"2"=list(edges=c("3")),
           "3"=list(edges=NULL)) ## allowed transitions between states
                                 ## edges=NULL implies terminal node

## Specifying tree
treeobj <- new("graphNEL", nodes=Nodes, edgeL=Edges,
                edgemode="directed")

## Running msSurv
ans1 <- msSurv(ex1,treeobj)
SOPt(ans1, t=2.5)

}

\keyword{survival}
