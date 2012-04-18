\name{Pst}
\alias{Pst}

\title{
Computation of P(s,t)
}
\description{
  This function calculates the transition probability matrix between any two values s and t and
  then prints it.  The function also
  calculates and prints the cov(P(s,t)) matrix if the user specifies.  Values are returned invisibly.
}
\usage{
Pst(object, s=0, t="last", deci=4, covar=FALSE)
}

\arguments{
  \item{object}{
    An msSurv object.
}
  \item{s}{
    The lower time.  Default is 0.
}
  \item{t}{
     The highter time.  Default is "last" which is the highest (or "last") event time.
}
  \item{deci}{
     Numeric argument specifying number of decimal places for estimates.  Default is 4.
}
  \item{covar}{
     Logical argument to determine if var(P(s,t)) is computed.  Default is FALSE.
}
}
\details{

   Computation of P(s,t) and var(P(s,t)) for multistate models are described in Andersen et al. (1993).  
   
}
\value{
  Returned invisibly:
  \describe{  
    \item{Pst}{The transition probability matrix between two times s and t.}
    \item{cov.Pst}{Covariance matrix for the transition probability matrix
      between two times s and t (if \code{covar==TRUE}).}
    }
}
\references{
  Andersen, P.K., Borgan, O., Gill, R.D. and Keiding,
  N. (1993). Statistical models based on counting processes. 
  Springer Series in Statistics. New York, NY: Springer. 
}
\author{
   Nicole Ferguson <nicole.ferguson@kennesaw.edu>,
   Guy Brock <guy.brock@louisville.edu>,
   Somnath Datta <somnath.datta@louisville.edu>
}
\note{
   If s = 0 and all subjects begin in the initial state at time
   0, then the top row of P(0,t) yields the state occupation
   probabilities at time t.  Data are assumed to follow a Markovian model.
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
Pst(ans1, s=0.25, t=2.5)

}

\keyword{survival}
