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
msSurv(Data, tree,cens.type="ind",LT=FALSE,d.var=FALSE,B=200,start.states)
}

\arguments{
  \item{Data}{
    Data with columns id, start, stop, start state, end state
}
  \item{tree}{
    A graphNEL graph with the nodes and edges of the multistate model.
}
  \item{cens.type}{
    A character string specifying whether censoring is independent ("ind") or state dependent ("dep").
    Default is "ind".
}
  \item{LT}{
    Logical argument.  Whether data is subject to left truncation.  Default is FALSE.
}
  \item{d.var}{
    Logical argument.  Whether variance is estimated for the state entry and exit time distributions.  Default is FALSE.
}
  \item{B}{
    The number of bootstraps for dependent censored data.
}
  \item{start.states}{
    A vector of starting states for individuals in the data set at time
    0.  The vector is automatically calculated for right censored data 
    if not specified.  For left truncated data, all individuals are
    assumed to be in the initial state at time 0 unless the user specifies
    differently.   
}
}
\details{
   Data are from a possibly non-Markovian multistate model with a directed tree
   structure and subject to right censoring and possibly left truncation.

   State 0 is reserved as the censoring state.

   State occupation probabilities are calculated according to the formula given
   in Datta & Satten (2001).

   State entry and exit time distributions are calculated using state occupation
   probabilities.

   The transition probability matrices are estimates using
   the Aalen-Johansen estimator described in Andersen et al. (1993).  The transition
   probability is the (i,j)th entry of the transition matrix.
   Datta and Satten (2001) showed that the Aalen-Johansen estimator remained
   valid for non-Markov systems.

   The covariance matrix for transition probability is computed componentwise
   following formula (4.4.20) in Andersen et al. (1993, p. 295) for independent
   censoring.  Bootstrapping is used to estimate the variance for state entry &
   exit time distributions.  Bootstrapping is also used to estimate the variance
   for data with dependent censoring.

}
\value{

  An object of S4 class  \code{\linkS4class{msSurv}} with the following slots:
  \item{tree}{A graphNEL graph with the nodes and edges of the multistate model.}
  \item{ns}{Number of states.}
  \item{et}{Event times.}
  \item{pos.trans}{Possible transitions between states.}
  \item{nt.states}{Non-terminal states.}
  \item{dNs}{The counting process.}
  \item{Ys}{The "at risk" set.}
  \item{ps}{State occupation probabilities.}
  \item{all.ajs}{An array containing matrices of Aalen-Johansen estimates.}
  \item{Fs}{State entry time distributions.}
  \item{Gs}{State exit time distributions.}
  \item{out}{Variance-covariance matrices of transition probability.}
  \item{cov.p}{Variance of state occupation probability.}
  \item{sum.dNs}{Counting process for total transitions out of each state.}
  \item{dNs.K}{The weighted counting process.}
  \item{Ys.K}{The weighted "at risk" set.}
  \item{sum.dNs.K}{Weighted counting process for total transitions out of each state.}
  \item{all.I_dA}{Array with all the I+dA matrices.}
  \item{Fs.var}{Variance of state entry time distribution.}
  \item{Gs.var}{Variance of state exit time distribution.}
}
\references{
  Andersen, P.K., Borgan, O., Gill, R.D. and Keiding, N. (1993). Statistical models based on counting
  processes. Springer Series in Statistics. New York, NY: Springer.

  Datta, S. and Satten G.A. (2001). Validity of the Aalen-Johansen estimators of stage occupation
  probabilities and Nelson-Aalen estimators of integrated transition hazards for non-Markov models.
  Statistics and Probability Letters, 55(4): 403-411.

  Datta S, Satten GA (2002). Estimation of Integrated Transition Hazards and Stage Occu-
  pation Probabilities for Non-Markov Systems under Dependent Censoring." Biometrics,
  58(4), 792-802.
}
\author{
   Nicole Ferguson <nicole.ferguson@louisville.edu>,
   Guy Brock <guy.brock@louisville.edu>,
   Somnath Datta <somnath.datta@louisville.edu>
}
\note{
   State 0 is reserved as the censoring state.
}


\seealso{
See the description of the plot, print, and summary methods in the help
page for S4 class \code{\linkS4class{msSurv}}.
}
\examples{
## 3-state illness-death multistate model (no left-truncation)

## data for 3 individuals
## data in the form "id","start","stop","st.stage","stage"
p1 <- c(1,0,0.21,1,3)
p2 <- c(2,0,0.799,1,2)
p22 <- c(2,0.799,1.577,2,3)
p3 <- c(3,0,0.199,1,0)

##combining data into a matrix
ex1 <- rbind(p1,p2,p22,p3)
colnames(ex1) <- c("id","start","stop","st.stage","stage")
ex1 <- data.frame(ex1)


## inputting nodes & edges of the tree structure
Nodes <- c("1","2","3") # states in MSM
Edges <- list("1"=list(edges=c("2","3")),"2"=list(edges=c("3")),
           "3"=list(edges=NULL)) #allowed transitions between states
                                 #edges=NULL implies terminal node

##specifying tree
treeobj <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

ans1 <- msSurv(ex1,treeobj)


## printing msSurv object 'ans1'
print(ans1)
print(ans1,covar=TRUE) #to print with covariance

#plotting msSurv object 'ans1'
plot(ans1,plot.type="stateocc")
plot(ans1,plot.type="stateocc",states=c("1","2"))
plot(ans1,plot.type="transprob")
plot(ans1,plot.type="entry.d")
plot(ans1,plot.type="exit.d")

## 3-state illness-death multistate model WITH left-truncation

## data for 3 individuals
## data in the form "id","start","stop","st.stage","stage"

p1 <-c(1, 0.383, 1.400, 1, 0)
p2 <-c(2, 0.698, 0.999, 1, 2)
p22<-c(2, 0.999, 1.180, 2, 0)
p3 <-c(3, 0.249, 0.391, 1, 2)
p32<-c(3, 0.391, 0.513, 2, 3)

ex2 <- rbind(p1,p2,p22,p3,p32)
colnames(ex2) <- c("id", "start", "stop", "st.stage", "stage")
ex2 <- data.frame(ex2)

## inputting nodes & edgest of the tree structure
Nodes <- c("1","2","3")
Edges <- list("1"=list(edges=c("2","3")),"2"=list(edges=c("3")),"3"=list(edges=NULL))
treeobj <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

ans2 <- msSurv(ex2,treeobj,LT=TRUE)

## Summary for msSurv object 'ans2'
summary(ans2) #gives estimates for IQR event times
summary(ans2,all=TRUE) #gives estimates for all eent times
 #summary for state occupation probability only
summary(ans2,trans.pr=FALSE)

}