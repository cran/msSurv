\name{msSurv}
\alias{msSurv}
\title{
Nonparametric Estimation for Multistate Models
}
\description{
The function uses the counting process and at risk set of event times to calculate the state occupation
probabilities, as well as the state entry and exit time distributions,
for a general, possibly non-Markov, model.  It also calculates the
transition probability matrices and covariance matrices for the state occupation probabilities
and transition probabilities.
}
\usage{
msSurv(Data, tree, cens.type="ind", LT=FALSE, bs=FALSE, B=200)
}

\arguments{
  \item{Data}{
    Data with counting-process style of input. Columns should be named
    "id", "start" (needed if \code{LT=TRUE}, optional
    otherwise), "stop",  "start.stage", and "end.stage".
}
  \item{tree}{
    A \code{\link[graph:graphNEL-class]{graphNEL}} graph with the nodes
    and edges of the multistate model.
}
  \item{cens.type}{
    A character string specifying whether censoring is independent
    ("ind") or state dependent ("dep"). Default is "ind".
}
  \item{LT}{
    Logical argument specifying whether data are subject to left truncation.
    Default is \code{FALSE}.
}
  \item{bs}{
    Logical argument specifying whether to use the bootstrap to
  calculate the variance.  The bootstrap is needed for the variance of
  the state entry/exit distributions, for state-dependent censoring, and
  for state occupation probabilities when there is more than one
  possible initial starting state.  Default is \code{FALSE}.
}
  \item{B}{
    The number of bootstrap iterations, when \code{bs=TRUE}.
}
}
\details{
   Data are from a possibly non-Markovian multistate model with a
   directed tree structure and subject to right censoring and possibly
   left truncation.

   State 0 is reserved as the censoring state.

   State occupation probabilities are calculated according to the
   formula given in Datta & Satten (2001).


   State entry and exit time distributions are calculated using state
   occupation probabilities.


   The transition probability matrices are estimates using
   the Aalen-Johansen estimator described in Andersen et al. (1993).
   The transition probability is the (i,j)th entry of the transition matrix.
   Datta and Satten (2001) showed that the Aalen-Johansen estimator remained
   valid for non-Markov systems.

   The covariance matrix for transition probability is computed componentwise
   following formula (4.4.20) in Andersen et al. (1993, p. 295) for
   independent censoring.  The bootstrap is needed for state-dependent
   censoring, and also for the variance of
   the state entry/exit distributions,  and
   for state occupation probabilities when there is more than one
   possible initial starting state.

}
\value{

  An object of S4 class  \code{\linkS4class{msSurv}} with the following slots:
  \item{tree}{A \code{\link[graph:graphNEL-class]{graphNEL}} object with
  the nodes and edges of the multistate model.}
  \item{ns}{Number of states.}
  \item{et}{Event times.}
  \item{pos.trans}{Possible transitions between states.}
  \item{nt.states}{Non-terminal states.}
  \item{dNs}{The counting process.}
  \item{Ys}{The "at-risk" set.}
  \item{sum_dNs}{Counting process for total transitions out of each state.}
  \item{dNs.K}{The weighted counting process.}
  \item{Ys.K}{The weighted "at-risk" set.}
  \item{sum_dNs.K}{Weighted counting process for total transitions out of each state.}
  \item{ps}{State occupation probabilities.}
  \item{AJs}{An array containing matrices of Aalen-Johansen estimates.}
  \item{I.dA}{Array with all the I+dA transition matrices for
    Aalen-Johansen computation.}
  \item{cov.AJs}{Variance-covariance matrices of the transition
    probabilities (A-J estimates).}
  \item{var.sop}{Variance of state occupation probability.}
  \item{cov.dA}{A matrix containing the covariance of dA matrices used
    for computation of cov(P(s,t)).}
  \item{Fnorm}{Normalized state entry time distributions.}
  \item{Fsub}{State entry time sub-distributions.}
  \item{Gnorm}{Normalized state exit time distributions.}
  \item{Gsub}{State exit time sub-distributions.}
  \item{Fnorm.var}{Variance of the normalized state entry time
  distributions.}
  \item{Fsub.var}{Variance of the state entry time sub-distributions.}
  \item{Gnorm.var}{Variance of normalized state exit time distributions.}
  \item{Gsub.var}{Variance of state exit time sub-distributions.}
}
\references{

  Nicole Ferguson, Somnath Datta, Guy Brock (2012).
  msSurv: An R Package for Nonparametric Estimation of Multistate Models.
  Journal of Statistical Software, 50(14), 1-24.
  URL http://www.jstatsoft.org/v50/i14/.

  Andersen, P.K., Borgan, O., Gill, R.D. and Keiding,
  N. (1993). Statistical models based on counting processes.
  Springer Series in Statistics. New York, NY: Springer.


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
   Somnath Datta <somnath.datta@louisville.edu>,
   Guy Brock <guy.brock@louisville.edu>
}
\note{
   State 0 is reserved as the censoring state.
}


\seealso{
See the description of the \code{plot}, \code{print}, and \code{summary}
methods in the help page for S4 class \code{\linkS4class{msSurv}}.
}
\examples{
## 3-state illness-death multistate model (no left-truncation)

## Row data for 3 individuals
## Data in the form "id", "start", "stop", "start.stage", "end.stage"
p1  <- c(1,0,0.21,1,3)
p2  <- c(2,0,0.799,1,2)
p22 <- c(2,0.799,1.577,2,3)
p3  <- c(3,0,0.199,1,0)

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

ans1 <- msSurv(ex1, treeobj)

## printing msSurv object 'ans1'
print(ans1)

## summary for msSurv object
summary(ans1)

## plotting msSurv object 'ans1'
plot(ans1, plot.type="stateocc")
plot(ans1, plot.type="stateocc", states=c("1", "2"))
plot(ans1, plot.type="transprob")
plot(ans1, plot.type="entry.norm")
plot(ans1, plot.type="exit.norm")

## 3-state illness-death multistate model WITH left-truncation

## Row data for 3 individuals
## Data in the form "id", "start", "stop", "start.stage", "end.stage"

p1  <- c(1,  0.383,  1.400,  1,  0)
p2  <- c(2,  0.698,  0.999,  1,  2)
p22 <- c(2,  0.999,  1.180,  2,  0)
p3  <- c(3,  0.249,  0.391,  1,  2)
p32 <- c(3,  0.391,  0.513,  2,  3)

ex2 <- rbind(p1, p2, p22, p3, p32)
colnames(ex2) <- c("id",  "start",  "stop",  "start.stage",  "end.stage")
ex2 <- data.frame(ex2)

## inputting nodes & edgest of the tree structure
Nodes <- c("1", "2", "3")
Edges <- list("1"=list(edges=c("2", "3")), "2"=list(edges=c("3")), "3"=list(edges=NULL))
treeobj <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")

ans2 <- msSurv(ex2, treeobj, LT=TRUE)

## Summary for msSurv object 'ans2'
summary(ans2) ## gives estimates for IQR event times
summary(ans2,  all=TRUE) ## gives estimates for all event times
## Summary for state occupation probability only
summary(ans2,  trans.pr=FALSE, dist=FALSE)

}