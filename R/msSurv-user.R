#########################################################
############# Transition Probability P(s,t) #############
#########################################################

## Function for P(s,t)
## Takes msSurv object, 1st (s), & last time (t)

Pst <- function(object,s=0,t="last",deci=4,covar=FALSE){

    if (!(0 <= s & s < t))
        stop("'s' and 't' must be positive, and s < t")
    if (t <= object@et[1] | s >= object@et[length(object@et)])
        stop("Either 's' or 't' is an invalid time")

	if(t=="last") t <- object@et[length(object@et)]
	
	#first <- length(object@et[object@et<= s]) + 1
	#last <- length(object@et[object@et<= t])

	idx <- which(s<=object@et & object@et<=t) #location of those [s,t]
	l.idx <- length(idx)

	cum.prod <- diag(object@ns)
	rownames(cum.prod) <- nodes(object@tree)
	red.all.ajs <- array(dim=c(object@ns,object@ns,nrow(object@dNs)),dimnames=list(rows=nodes(object@tree),cols=nodes(object@tree),dim=rownames(object@dNs)))

	for(i in idx){
		cum.prod <- cum.prod %*% object@all.I_dA[,,i]
	     	red.all.ajs[,,i] <- cum.prod
	}

    if(covar == TRUE){

	bl.Id <- diag(1,(object@ns)^2) #Ident matrix for Kronecker product
 	var.Pst <- array(0, dim=c(dim(object@all.I_dA[,,idx])[c(1, 2)]^2,nrow(object@dNs)))
      		colnames(var.Pst) <- rownames(var.Pst) <- paste(rep(nodes(object@tree),object@ns),sort(rep(nodes(object@tree),object@ns)))
 	Id <- diag(1,object@ns)

	for(i in idx){
	  if(i==idx[1]) var.Pst[, , i] <- bl.Id%*% object@cov.dA[,,i] %*% bl.Id
	  else var.Pst[, , i] <- (t(object@all.I_dA[, , i]) %x% Id) %*% var.Pst[, , i-1] %*%((object@all.I_dA[, , i]) %x% Id) +
	      (Id %x% red.all.ajs[, , i-1]) %*% object@cov.dA[,,i]  %*% (Id%x% t(red.all.ajs[, , i-1]))
	  #Note:  red.all.ajs[,,i-1] = P(s,t-), object@cov.dA[,,i] = the varcov (cov.dA) for the 1st time (no products), object@all.I_dA is I+dA matrix from original program
	} #end of for idx
    } #end of if var 

#	list(P.s.t=cum.prod)

	cat(paste("Estimate of P(",s,",",t,")\n", sep = ""))
       	print(round(cum.prod,digits=deci))
#	print(red.all.ajs[,,idx])
       	cat("\n")

       if (!is.null(object@out) & covar == TRUE) {
           cat(paste("Estimate of cov(P(",s,",",t,"))\n", sep = ""))
           print(round(var.Pst[,,max(idx)],digits=deci))
       }


}

########################################################
############# State Occ for Specific Time t ############
########################################################


st.t<- function(object,t="last",deci=4,covar=FALSE){

	if(t=="last") t <- object@et[length(object@et)]
	t.loc<- length(object@et[object@et<= t])

	cat(paste("The state occupation probabilities at time ",t," are:\n", sep = ""))
	for(i in nodes(object@tree)){
	   cat(paste("State ",i,": ",round(object@ps[t.loc,as.numeric(i)],deci),"\n",sep = ""))
	}
	cat("\n")

       if (!is.null(object@out) & covar == TRUE) {
           cat(paste("Covariance Estimates for State Occupation Probability: \n", sep = ""))

	   for(i in nodes(object@tree)){
	      cat(paste("State ",i,": ",round(object@cov.p[t.loc,as.numeric(i)],deci),"\n",sep = ""))
	   }
       }


}

#1/17/11 - Adding state entry/exit distribution function
#########################################################
########## State Entry/Exit Time Distribution ###########
######### Calculates BS Variance Estimate & CIs #########
##### Also gives est of Fs & Gs at specified time t #####
#########################################################

EntryExit <- function(object,t="last",deci=4,covar=FALSE){

## Add entry & exit arguments to display different information?

	if(covar==TRUE & is.null(object@Fs.var)){
		stop(paste("msSurv object does not have variance estimates for entry/exit time distributions.
Please re-run the msSurv object with the argument 'd.var=FALSE' and then try again. \n", sep=""))
	}

      entry.st <- which(!(sapply(inEdges(object@tree), function(x) length(x) == 0)))
	initial <- which(!(nodes(object@tree)%in%entry.st)) #initial states, no Fs
	exit.st <- which(sapply(edges(object@tree), function(x) length(x) > 0))
	terminal <- which(!(nodes(object@tree)%in%exit.st))

	if(t=="last") t <- object@et[length(object@et)]
	t.loc<- length(object@et[object@et<= t])

	cat(paste("The state entry distributions at time ",t," are:\n", sep = ""))
	for(i in entry.st){
	   cat(paste("State ",i,": ",round(object@Fs[t.loc,][[i]],deci),"\n",sep = ""))

	}
	   cat(paste("State entry distributions for state ",as.character(initial),"is omitted
since there are no transitions into that state."))
	cat("\n","\n")

       if (covar==TRUE) {

           cat(paste("Variance Estimates for State Entry Distributions: \n", sep = ""))

	   for(i in entry.st){
	      cat(paste("State ",i,": ",round(object@Fs.var[t.loc,][[i]],deci),"\n",sep = ""))

	   } #end of entry.st loop

	   cat("Variance estimates of state entry distributions for state ",as.character(initial),"is omitted
since there are no transitions into that state.")
 	   cat("\n")

       } #end of if covar loop


       cat("\n","\n")

	cat(paste("The state exit distributions at time ",t," are:\n", sep = ""))
	for(i in exit.st){
	   cat(paste("State ",i,": ",round(object@Gs[t.loc,][[i]],deci),"\n",sep = ""))

	}
	cat("State exit distributions for state(s) ",as.character(terminal),"is (are) omitted
since there are no transitions into that (those) state(s).")
	cat("\n","\n")


       if (covar==TRUE) {

           cat(paste("Variance Estimates for State Exit Distributions: \n", sep = ""))

	   for(i in exit.st){
	      cat(paste("State ",i,": ",round(object@Gs.var[t.loc,][[i]],deci),"\n",sep = ""))
	   } #end of entry.st loop

	   cat("Variance estimates of state exit distributions for state(s) ",as.character(terminal),"is (are) omitted
since there are no transitions into that (those) state(s).")

	   cat("\n")

       } #end of if covar loop



} #end of loop


############################################################
####################### Main Function ######################
############################################################


msSurv <- function(Data,tree,cens.type="ind",LT=FALSE,d.var=FALSE,B=200,start.states){
## cens.type = "ind" or "dep" - used for D-S calculation

   if (any(!(c("id", "stop", "st.stage", "stage")%in%colnames(Data))))
       stop("'Incorrect column names for 'Data'.  Column names should be 'id','stop','st.stage', or 'stage'.")

   if(!("start" %in% colnames(Data)) & LT==TRUE)
	stop("The 'start' times must be specified for left truncated data.")

   if(!("start" %in% colnames(Data)) & LT==FALSE)  Data=Add.start(Data)
 

## NOTE: If LT=TRUE and missing start.states, need to assume all start.states = initial state in tree
  if (missing(start.states) & LT == TRUE) {	
	start.probs <- numeric(length(nodes(tree)))
	names(start.probs) <- nodes(tree)
	start.probs[names(start.probs)== nodes(tree)[which(sapply(inEdges(tree), function(x) !length(x)>0))]] <- which(sapply(inEdges(tree), function(x) !length(x)>0)) 
      warning("'start.states' not specified.  Assuming all individuals start in the initial state at time 0.")
  }

##   if(missing(start.states)){
##	if(LT==TRUE)
##	 warning("'start.states' not specified.  Assuming all individuals start in the initial state at time 0.")
##   }

   if(!missing(start.states) & LT==TRUE){
	start.probs <- numeric(length(nodes(tree)))
	names(start.probs) <- nodes(tree)
	start.probs[names(start.probs)== names(table(start.states))] <- table(start.states)/length(start.states)
   }

   n <- length(unique(Data$id)) ## number of individuals in sample
   ns <- length(nodes(tree)) ## number of states

   Cens <- Add.States(tree)
   if(LT) {
	Data = LT.Data(Data)
	cp <- CP(tree,Cens$treeLT,Data,Cens$nt.states.LT)
   }

   if(!LT) cp <- CP(tree,Cens$tree0,Data,Cens$nt.states)

    ds.est<-DS(LT="LT",Cens$nt.states,cp$dNs,cp$sum.dNs,cp$Ys,Cens="0",cens.type)
   ## May want to include LT & Cens values in the data conversion function we plan to do
   ## Here we want to make sure the nt.states is for Cens and not LT!!
    cp.red <- Red(tree,cp$dNs,cp$Ys,cp$sum.dNs,ds.est$dNs.K,ds.est$Ys.K,ds.est$sum.dNs.K)

    if(missing(start.states)){
	 if(!LT)  start.probs=cp$start.probs
    }

    et <- as.numeric(rownames(cp.red$dNs))

   ##made a list of all possible transitions excluding thansition for terminal state into terminal state ...
    res.ci2 <- strsplit(colnames(cp.red$dNs), " ") #string splits names
    a <- sapply(res.ci2, function(x) x[2]) #pull 1st number
    b <- sapply(res.ci2, function(x) x[3]) #pull 2nd number
    pos.trans <- paste(a,b) #combine 1st & 2nd numbers
    stay <- paste(Cens$nt.states,Cens$nt.states) #creating names for 11,22, etc. (except staying in terminal states)
    pos.trans <- sort(c(stay,pos.trans)) #sorting all possible "transitions" from P(s,t)

   stateoccfn <- stocc(ns,tree,cp.red$dNs.K,cp.red$Ys.K,start.probs)
##11/10 added .K to the dNs & Ys to follow formula in DS ... also applicable to ind/dep cens
#   stateoccfn <- stocc(ns,tree,cp.red$dNs,cp.red$Ys,n)
   ent.exit <- Dist(stateoccfn$ps,ns,tree)

   variances <- var.fn(tree,ns,Cens$nt.states,cp.red$dNs,cp.red$Ys,cp.red$sum.dNs,stateoccfn$all.ajs, stateoccfn$all.I_dA,stateoccfn$ps)

## length(which(table(cp$res)>0)) ## use this to determine when to BS the "ind" cens var(ps) if>1 then need to Bs or if the # is not in the initial state need to BS
##1/26/11 - Add ps variances if different starting states
	no.start.st <-  length(which(start.probs>0))


   if(cens.type=="ind" & no.start.st==1){
	## will need to update tree location when it is no longer global
	if(d.var==TRUE){	
		ee.vars <- BS.var(Data,tree,ns,et,cp.red$dNs,cens.type,B,LT,start.states)
		var.Fs <- ee.vars$Fs
	 	var.Gs <- ee.vars$Gs
	} else {
	    var.Fs=NULL
	    var.Gs=NULL
	}
	   res <- new("msSurv", tree=tree,ns=ns,et=et,pos.trans=pos.trans,nt.states=Cens$nt.states,dNs=cp.red$dNs,
				Ys=cp.red$Ys,ps=stateoccfn$ps,all.ajs=stateoccfn$all.ajs,Fs=ent.exit$Fs,Gs=ent.exit$Gs,
				out=variances$out,cov.p=variances$cov.p,sum.dNs=cp.red$sum.dNs, dNs.K=cp.red$dNs.K,Ys.K=cp.red$Ys.K,
				sum.dNs.K=cp.red$sum.dNs.K,cov.dA=variances$varcov,all.I_dA=stateoccfn$all.I_dA,
                        	Fs.var=var.Fs,Gs.var=var.Gs)
	}

   if(cens.type=="ind" & no.start.st>1){
	   bsvar <- BS.var(Data,tree,ns,et,cens.type,B,LT,start.states)
	   res <- new("msSurv", tree=tree,ns=ns,et=et,pos.trans=pos.trans,nt.states=Cens$nt.states,dNs=cp.red$dNs,
				Ys=cp.red$Ys,ps=stateoccfn$ps,all.ajs=stateoccfn$all.ajs,Fs=ent.exit$Fs,Gs=ent.exit$Gs,
				out=bsvar$out,cov.p=bsvar$cov.p,sum.dNs=cp.red$sum.dNs, dNs.K=cp.red$dNs.K, Ys.K=cp.red$Ys.K,
				sum.dNs.K=cp.red$sum.dNs.K,cov.dA=variances$varcov,all.I_dA=stateoccfn$all.I_dA,
				Fs.var=bsvar$var.Fs, Gs.var=bsvar$var.Gs)}
   

  if(cens.type=="dep"){
	   bsvar <- BS.var(Data,tree,ns,et,cens.type,B,LT,start.states)
	   res <- new("msSurv", tree=tree,ns=ns,et=et,pos.trans=pos.trans,nt.states=Cens$nt.states,dNs=cp.red$dNs,
				Ys=cp.red$Ys,ps=stateoccfn$ps,all.ajs=stateoccfn$all.ajs,Fs=ent.exit$Fs,Gs=ent.exit$Gs,
				out=bsvar$out,cov.p=bsvar$cov.p,sum.dNs=cp.red$sum.dNs, dNs.K=cp.red$dNs.K, Ys.K=cp.red$Ys.K,
				sum.dNs.K=cp.red$sum.dNs.K,cov.dA=variances$varcov,all.I_dA=stateoccfn$all.I_dA,
				Fs.var=bsvar$var.Fs,Gs.var=bsvar$var.Gs)
	}
   

	

   return(res)
}





