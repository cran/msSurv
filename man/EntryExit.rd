\name{EntryExit}
\alias{EntryExit}

\title{
State Entry and Exit Time Distributions at Time t
}
\description{
  This function displays the state entry and exit time distributions at
  time t.  The function also displays the corresponding 
  variance estimates if the user requests them.
}
\usage{
EntryExit(object, t="last", deci=4, covar=FALSE, norm=TRUE)
}

\arguments{
  \item{object}{
    A msSurv object.
  }
  \item{t}{
    The time to find state entry/exit distributions.  Default is "last" which is the highest (or "last") event time.
  }
  \item{deci}{
    Numeric argument specifying number of decimal places for estimates.  Default is 4.
  }
  \item{covar}{
    Logical argument to determine if the variance is displayed.  Default is FALSE.  If variance estimates are NULL, an error message will print.
  }
  \item{norm}{
    Logical argument to determine whether normalized or
  non-normalized (subdistribution) functions are displayed.  Default is
  normalized distributions.  
  }
}

\details{
   Display of the state entry and exit time distributions and
   corresponding variance for multistate models at each state in the
   system where computation makes sense.   
}
\value{
   Returns (invisibly) a list consisting of the state entry / exit
   distributions (either \code{entry.norm} and \code{exit.norm} or
   \code{entry.sub} and \code{exit.sub}), and (optionally) the variance
   estimates (either \code{entry.var.norm} and \code{exit.var.norm} or
   \code{entry.var.sub} and \code{exit.var.sub}).
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
\note{
   State entry distributions (and corresponding variance estimates) are
   displayed for all states where entry into the state occurs. 
   State exit distributions (and corresponding variance estimates) are
   displayed for all states where exit from the state occurs. 
   
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
EntryExit(ans1,t=0.8)

}



\keyword{survival}
