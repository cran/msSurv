\name{msSurv-class}
\Rdversion{1.1}
\docType{class}
\alias{msSurv-class}
\alias{plot,msSurv,missing-method}
\alias{print,msSurv-method}
\alias{summary,msSurv-method}
\alias{tree,msSurv-method}
\alias{tree}

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
      \code{graphNEL} object with nodes corresponding to the states in
      the multistate model and the edges corresponding to the allowed
      transitions. }
    \item{\code{ns}:}{Object of class \code{"numeric"}. The number of
      unique states in the multistate model.}
    \item{\code{et}:}{Object of class \code{"numeric"}. The event times.}
    \item{\code{pos.trans}:}{Object of class \code{"character"}.
      Possible transtitions between states. }
    \item{\code{nt.states}:}{Object of class \code{"numeric"}.  The
      non-terminal states in the multistate model. }
    \item{\code{dNs}:}{Object of class \code{"array"}.  A matrix
      containing the counting processes for the event times. }
    \item{\code{Ys}:}{Object of class \code{"array"}.  A matrix
      containing the at risk sets for the event times. }
    \item{\code{ps}:}{Object of class \code{"array"}.  A matrix with
      state occupation probabilities for each state at each event time.}
    \item{\code{all.ajs}:}{Object of class \code{"array"}.  An array
      containing matrices of Aalen-Johansen estimates at each event time. }
    \item{\code{Fs}:}{Object of class \code{"array"}.  A matrix
      containing state entry distributions.  Note: "NA" is recorded for
      Fs at the initial state (node).}
    \item{\code{Gs}:}{Object of class \code{"array"}.  A matrix
      containing state exit distributions.  Note: "NA" is recorded for
      Fs at terminal states (nodes).}
    \item{\code{out}:}{Object of class \code{"array"}.  An array
      containing the variance-covariance matrices for transition
      probabilities at each event time.}
    \item{\code{cov.p}:}{Object of class \code{"array"}.  A matrix
      containing variance estimates for state occupation probability. }
    \item{\code{sum.dNs}:}{Object of class \code{"array"}.  A matrix
      containing the counting process for total transitions out of each
      state. }
    \item{\code{dNs.K}:}{Object of class \code{"array"}. A matrix
      containing the weighted counting process described in Datta and
      Satten (2001).}
    \item{\code{Ys.K}:}{Object of class \code{"array"}. A matrix
      containing the weighted at risk sets described in Datta and
      Satten (2001).  }
    \item{\code{sum.dNs.K}:}{Object of class \code{"array"}. A matrix
      containing the weighted counting process for total transitions out
      of each state.}
    \item{\code{all.I_dA}:}{Object of class \code{"array"}. A matrix
      containing the I+dA matrices for Aalen-Johansen computation.}
    \item{\code{cov.dA}:}{Object of class \code{"array"}. A matrix
      containing the covariance of dA matrices used for computation
      of cov(P(s,t)).}
    \item{\code{Fs.var}:}{Object of class \code{"array or NULL"}.  A matrix
      containing variance estimates for state entry distribution. May be NULL
      if the user does not specify \code{d.var=TRUE}.}
    \item{\code{Gs.var}:}{Object of class \code{"array or NULL"}.  A matrix
      containing variance estimates for state exit distribution. May be NULL
      if the user does not specify \code{d.var=TRUE}.}


  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "msSurv", y = "missing")}: Plotting
      method for "msSurv" objects.

      Additional arguments:
      \describe{
	\item{\code{states="ALL"}}{States in the multistate
	  model to be plotted.  Default is all states in the system.
          User may specify individual states or multiple states to plot. }
	\item{\code{trans="ALL"}}{Transitions in the
	  multistate model to be plotted.  Default is all transitions.  
          Transitions should be entered without a space between the two states,
          ie: "11".}
	 \item{\code{CI=TRUE}}{A logical argument to specify whether
	   pointwise confidence intervals should be plotted.  If the
      user specifies CI=FALSE, only the estimates are plotted.  If the user
      specifies CI=TRUE, plots of each estimate and its corresponding 
      confidence intervals are created. The default is TRUE. }
	 \item{\code{ci.level=0.95}}{Confidence level.  Default is 0.95.}
	 \item{\code{ci.trans="linear"}}{Transformation applied to
	     confidence intervals. Possible choices are "linear", "log",
	     "log-log", and "cloglog".  Default is "linear".}
	 \item{\code{plot.type="stateocc"}}{Determines the type of estimate to be
      plotted.  User may specify "transprob" for transition probability
      plots, "stateocc" for state occupation probability plots, "entry.d" for state
      entry time distributions, or "exit.d" for state exit time distributions.  "stateocc" is the default.
      Note: CIs do not plot for state entry/exit time distributions.}

         \item{\code{\dots}}{Further arguments}

     }}

    \item{print}{\code{signature(x = "msSurv")}: Print function for
      "msSurv" objects.

      Additional arguments:
      \describe{
	\item{\code{covar=FALSE}}{A logical argument to specify whether
	  to print the variance-covariance estimates.  Default is FALSE.}
	\item{\code{ee.distn=TRUE}}{A logical argument to specify
      whether to print the entry and exit time distributions.  Default
      is TRUE.}

	}
    }

    \item{summary}{\code{signature(object = "msSurv")}: Summary function
      for "msSurv" objects.

      \describe{
      Additional arguements:
      \item{\code{digits=3}}{The number of significant digits to use for
	estimates.  Defalt is 3.}
      \item{\code{all=FALSE}}{Logical argument to determine whether
	summary information should be displayed for all event times or
	only for the key percentile time points (IQR).  Default is FALSE
	where all=FALSE corresponds to only the IQR of event times being
	displayed in the summary output.}
      \item{\code{ci.fun="linear"}}{Transformation applied to
	     confidence intervals. Possible choices are "linear", "log",
	     "log-log", and "cloglog".  Default is "linear".}
      \item{\code{ci.level=0.95}}{Confidence level.  Default is 0.95.}

    }}

    \item{tree}{\code{signature(object = "msSurv")}: Accesses the
      "tree" slot in "msSurv" objects}
	 }
}

\references{
  Andersen, P.K., Borgan, O., Gill, R.D. and Keiding, N. (1993). Statistical models based on counting   processes. Springer Series in Statistics. New York, NY: Springer.

  Datta, S. and Satten G.A. (2001). Validity of the Aalen-Johansen
  estimators of stage occupation probabilities and Nelson-Aalen estimators of integrated transition hazards for non-Markov models.  Statistics and Probability Letters, 55(4): 403-411.

}

\author{
    Nicole Ferguson <nicole.ferguson@louisville.edu>,
    Guy Brock <guy.brock@louisville.edu>,
    Somnath Datta <somnath.datta@louisville.edu>
}
\note{
    This function prints estimates of state occupation probabilities for each event time and a matrix giving the estimates of transition probability P(0,t).

    Summary information for the state occupation probabilities include the event time, estimate of state occupation probability, variance estimate, and lower and upper confidence intervals, and state entry and exit time distributions for each state.

    Summary information for the transition probabilities include the event time, estimate of transition probability,
    variance stimate, lower and upper confidence intervals and the number at risk using methods in Andersen et al. (1993) ("n.risk") and Datta & Satten (2001).  For transitions into a different state,
    the number of events is also provided according to both methods previously described.  For transitions into the same state, the number remaining after the event time is
    also provided.
}


\seealso{
  For a description of the function 'msSurv'  see \code{\link{msSurv}}.
  
}
\examples{

p1 <- c(1,0,0.21,1,3)
p2 <- c(2,0,0.799,1,2)
p22 <- c(2,0.799,1.577,2,3)
p3 <- c(3,0,0.199,1,0)

##combining data into a matrix
ex1 <- rbind(p1,p2,p22,p3)
colnames(ex1) <- c("id","start","stop","st.stage","stage")
ex1 <- data.frame(id=ex1[,1],start=ex1[,2],stop=ex1[,3],st.stage=ex1[,4],stage=ex1[,5])


## inputting nodes & edgest of the tree structure
Nodes <- c("1","2","3") # states in MSM
Edges <- list("1"=list(edges=c("2","3")),"2"=list(edges=c("3")),
           "3"=list(edges=NULL)) #allowed transitions between states
                                 #edges=NULL implies terminal node

##specifying tree
treeobj <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

ans1 <- msSurv(ex1,treeobj)

summary(ans1) #prints IQR for ans1
summary(ans1,all=TRUE) #prints all event times for ans1

#prints only state occupation probability info for all event times
summary(ans1,all=TRUE,trans.pr=FALSE)

plot(ans1) #plots state occupation probability
plot(ans1,states="1") 
plot(ans1,states=c("1","2"))
plot(ans1,plot.type="transprob") #plots for transition probability 
plot(ans1,plot.type="transprob",trans="12")
plot(ans1,plot.type="transprob",trans=c("12","13"))

}
\keyword{classes}
\keyword{survival}
