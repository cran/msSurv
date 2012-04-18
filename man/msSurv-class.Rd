\name{msSurv-class}
\Rdversion{1.1}
\docType{class}
\alias{msSurv-class}
\alias{plot,msSurv,missing-method}
\alias{print,msSurv-method}
\alias{summary,msSurv-method}
\alias{tree}
\alias{tree,msSurv-method}
\alias{ns}
\alias{ns,msSurv-method}
\alias{et}
\alias{et,msSurv-method}
\alias{pos.trans}
\alias{pos.trans,msSurv-method}
\alias{nt.states}
\alias{nt.states,msSurv-method}
\alias{dNs}
\alias{dNs,msSurv-method}
\alias{Ys}
\alias{Ys,msSurv-method}
\alias{sum.dNs}
\alias{sum.dNs,msSurv-method}
\alias{dNs.K}
\alias{dNs.K,msSurv-method}
\alias{Ys.K}
\alias{Ys.K,msSurv-method}
\alias{sum.dNs.K}
\alias{sum.dNs.K,msSurv-method}
\alias{ps}
\alias{ps,msSurv-method}
\alias{AJs}
\alias{AJs,msSurv-method}
\alias{cov.ps}
\alias{cov.ps,msSurv-method}
\alias{cov.AJs}
\alias{cov.AJs,msSurv-method}
\alias{I.dA}
\alias{I.dA,msSurv-method}
\alias{cov.dA}
\alias{cov.dA,msSurv-method}
\alias{Fnorm}
\alias{Fnorm,msSurv-method}
\alias{Fsub}
\alias{Fsub,msSurv-method}
\alias{Gnorm}
\alias{Gnorm,msSurv-method}
\alias{Gsub}
\alias{Gsub,msSurv-method}
\alias{Fnorm.var}
\alias{Fnorm.var,msSurv-method}
\alias{Fsub.var}
\alias{Fsub.var,msSurv-method}
\alias{Gnorm.var}
\alias{Gnorm.var,msSurv-method}
\alias{Gsub.var}
\alias{Gsub.var,msSurv-method}


\title{Class "msSurv"}
\description{The class \code{"msSurv"} contains nonparametric estimates
  for multistate models subject to right censoring and possibly left
  truncation by calling \code{\link{msSurv}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("msSurv", ...)}.
"msSurv" objects are also returned from function \code{\link{msSurv}}
}
\section{Slots}{
  \describe{
    \item{\code{tree}:}{Object of class \code{"graphNEL"}. A
      \code{\link[graph:graphNEL-class]{graphNEL}} object with nodes
      corresponding to the states in 
      the multistate model and the edges corresponding to the allowed
      transitions. }
    \item{\code{ns}:}{Object of class \code{"numeric"}. The number of
      unique states in the multistate model.}
    \item{\code{et}:}{Object of class \code{"numeric"}. The event times.}
    \item{\code{pos.trans}:}{Object of class \code{"character"}.
      Possible transtitions between states. }
    \item{\code{nt.states}:}{Object of class \code{"character"}.  The
      non-terminal states in the multistate model. }
    \item{\code{dNs}:}{Object of class \code{"array"}.  A matrix
      containing the differential of the counting processes
      for the event times. }
    \item{\code{Ys}:}{Object of class \code{"array"}.  A matrix
      containing the at risk sets for the event times. }
    \item{\code{sum.dNs}:}{Object of class \code{"array"}.  A matrix
      containing the differential for the counting process for total
      transitions out of each state, at each event time. }
    \item{\code{dNs.K}:}{Object of class \code{"array"}. A matrix
      containing the differential of the weighted counting process
      described in Datta and Satten (2001).}
    \item{\code{Ys.K}:}{Object of class \code{"array"}. A matrix
      containing the weighted at risk sets described in Datta and
      Satten (2001).  }
    \item{\code{sum.dNs.K}:}{Object of class \code{"array"}. A matrix
      containing the differential of the weighted counting process
      for total transitions out of each state, at each event time.}
    \item{\code{ps}:}{Object of class \code{"array"}.  A matrix with
      state occupation probabilities for each state at each event time.}
    \item{\code{AJs}:}{Object of class \code{"array"}.  An array
      containing matrices of Aalen-Johansen estimates (transition
      probabilities) at each event time. }
    \item{\code{I.dA}:}{Object of class \code{"array"}. A matrix
      containing the I+dA transition matrices for Aalen-Johansen computation.}
    \item{\code{cov.AJs}:}{Object of class \code{"array"}.  An array
      containing the variance-covariance matrices for transition
      probabilities at each event time.}
    \item{\code{var.sop}:}{Object of class \code{"array"}.  A matrix
      containing covariance estimates for the state occupation
      probabilities.} 
    \item{\code{cov.dA}:}{Object of class \code{"array"}. A matrix
      containing the covariance of dA matrices used for computation
      of cov(P(s,t)).}
    \item{\code{Fnorm}:}{Object of class \code{"array"}.  A matrix
      containing normalized state entry distributions.  Note: "NA" is
      recorded for \code{Fnorm} at the initial state(s) (node(s)).}
    \item{\code{Fsub}:}{Object of class \code{"array"}.  A matrix
      containing unnormalized (subdistribution) state entry
      distributions.  Note: "NA" is recorded for \code{Fsub} at the
      initial state(s) (node(s)).}
    \item{\code{Gnorm}:}{Object of class \code{"array"}.  A matrix
      containing normalized state exit distributions.  Note: "NA" is
      recorded for \code{Gnorm} at the terminal state(s) (node(s)).}
    \item{\code{Gsub}:}{Object of class \code{"array"}.  A matrix
      containing unnormalized (subdistribution) state exit
      distributions.  Note: "NA" is recorded for \code{Gsub} at the
      terminal state(s) (node(s)).}
    \item{\code{Fnorm.var}:}{Object of class \code{"array or NULL"}.  A matrix
      containing variance estimates for the normalized state entry
      distributions.  Will be \code{NULL} if the user does not specify
      \code{bs=TRUE}.}
    \item{\code{Fsub.var}:}{Object of class \code{"array or NULL"}.  A matrix
      containing variance estimates for the unnormalized
      (subdistribution) state entry distributions.  Will be \code{NULL}
      if the user does not specify \code{bs=TRUE}.}
    \item{\code{Gnorm.var}:}{Object of class \code{"array or NULL"}.  A matrix
      containing variance estimates for the normalized state exit
      distributions.  Will be \code{NULL} if the user does not specify
      \code{bs=TRUE}.}
    \item{\code{Gsub.var}:}{Object of class \code{"array or NULL"}.  A matrix
      containing variance estimates for the unnormalized
      (subdistribution) state exit distributions.  Will be \code{NULL}
      if the user does not specify \code{bs=TRUE}.}      
  }
}
\section{Methods}{
  \describe{
    \item{ACCESSOR FUNCTIONS}{\code{signature(object = "msSurv")}:
      Accessor functions are defined for each of the slots in an
      \code{msSurv} object, e.g. \code{tree}, \code{ns}, \code{et}, etc.
      The accessor functions all have the same name as the corresponding
      slot name, and all have the same signature.}
    
     \item{print}{\code{signature(x = "msSurv")}: Print method for
       "msSurv" objects.}

     \item{show}{\code{signature(object = "msSurv")}: Show method for
       "msSurv" objects.}
    
    \item{summary}{\code{signature(object = "msSurv")}: Summary function
      for "msSurv" objects.

      Additional arguements:
      \describe{
      \item{\code{digits=3}}{The number of significant digits to use for
	estimates.  Defalt is 3.}
      \item{\code{all=FALSE}}{Logical argument to determine whether
	summary information should be displayed for all event times or
	only for the key percentile time points (IQR).  Default is FALSE
	where \code{all=FALSE} corresponds to only the IQR of event
	times being displayed in the summary output.}
      \item{\code{times=NULL}}{Numeric vector of time-points at which to
	present summary information.  Overrides \code{all} if supplied.}
       \item{\code{ci.fun="linear"}}{Transformation applied to
	confidence intervals. Possible choices are "linear", "log",
	"log-log", and "cloglog".  Default is "linear".}
      \item{\code{ci.level=0.95}}{Confidence level.  Default is 0.95.}
      \item{\code{stateocc=TRUE}}{Logical argument specifying whether
	state occupation probabilities should be displayed.  Default is
	\code{TRUE}.}
      \item{\code{trans.pr=TRUE}}{Logical argument specifying whether
	state transition probabilities should be displayed.  Default is
	\code{TRUE}.} 
      \item{\code{dist=TRUE}}{Logical argument specifying whether
	state entry / exit distributions should be displayed.  Default is
	\code{TRUE}.} 
      \item{\code{DS=FALSE}}{Logical argument specifying whether
	Datta-Satten weighted counting processes.
	Default is \code{FALSE}.} 
    }}

    \item{plot}{\code{signature(x = "msSurv", y = "missing")}: Plotting
      method for "msSurv" objects.

      Additional arguments:
      \describe{
	\item{\code{states="ALL"}}{States in the multistate
	  model to be plotted.  Default is all states in the system.
          User may specify individual states or multiple states to plot. }
	\item{\code{trans="ALL"}}{Transitions in the
	  multistate model to be plotted.  Default is all transitions.  
          Transitions should be entered with a space between the two states,
          e.g.: "1 1".}
	 \item{\code{CI=TRUE}}{A logical argument to specify whether
	   pointwise confidence intervals should be plotted.  If the
	   user specifies \code{CI=FALSE}, only the estimates are
	   plotted.  If the user specifies \code{CI=TRUE}, plots of each
	   estimate and its corresponding   
	   confidence intervals are created (if appropriate variances
	   are available). The default is \code{TRUE}. }
	 \item{\code{ci.level=0.95}}{Confidence level.  Default is 0.95.}
	 \item{\code{ci.trans="linear"}}{Transformation applied to
	   confidence intervals. Possible choices are "linear", "log",
	   "log-log", and "cloglog".  Default is "linear".}
	 \item{\code{plot.type="stateocc"}}{Determines the type of
	   estimate to be 
	   plotted.  User may specify "transprob" for transition probability
	   plots, "stateocc" for state occupation probability plots,
	   "entry.norm"  / "entry.sub"  for normalized / unnormalized
	   state entry time distributions, 
	   or "exit.norm"  / "exit.sub"  for normalized / unnormalized
	   state exit time distributions.  "stateocc" is the default.
	 }
         \item{\code{\dots}}{Further arguments passed to
	   \code{\link[lattice:xyplot]{xyplot}}} 
     }}
   }
}


\references{
  Andersen, P.K., Borgan, O., Gill, R.D. and Keiding,
  N. (1993). Statistical models based on counting processes. Springer
  Series in Statistics. New York, NY: Springer. 

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

    The \code{summary} function prints estimates of the transition
    probabilities P(0,t), state occupation probabilities, and state entry
    and exit time distributions for each state. Users can explicitly 
    select which information they want to display.

    Summary information for the transition probabilities include the
    event time, estimate of transition probability, 
    variance stimate, lower and upper confidence intervals and the
    number at risk using methods in Andersen et al. (1993) (\code{"n.risk"})
    and, optionally, the weighted Datta & Satten (2001) estimates
    (\code{"n.risk.K"}). 
    For transitions into a different state, the number of events is also
    provided according to both methods previously described.  For
    transitions into the same state, the number remaining after the
    event time is also provided.  
    
    Summary information for the state occupation probabilities include
    the event time, estimate of state occupation probability, variance
    estimate, and lower and upper confidence intervals.

    Estimates for both normalized and non-normalized (sub-distribution)
    state entry and exit distributions are also displayed.
    
}


\seealso{
  For a description of the function \code{'msSurv'}  see \code{\link{msSurv}}.
  
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

print(ans1)  ## same as 'show(ans1)'

summary(ans1) ## prints IQR for ans1
summary(ans1, all=TRUE) ## prints all event times for ans1

## prints only state occupation probability info for all event times
summary(ans1, all=TRUE, trans.pr=FALSE, dist=FALSE)

plot(ans1) ## plots state occupation probability
plot(ans1, states="1") 
plot(ans1, states=c("1", "2"))
plot(ans1, plot.type="transprob") ## plots for transition probability 
plot(ans1, plot.type="transprob", trans=c("1 2", "1 3"))

}
\keyword{classes}
\keyword{survival}
