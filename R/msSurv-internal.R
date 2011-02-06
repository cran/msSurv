############################################################
#################### Adding Start Times ####################
############################################################

Add.start <- function(Data){

	Data$start <- 0
	idx <- which(table(Data$id)>1)

	for(i in idx){
	   ab <-Data[which(Data$id==i),]
	   ab<-with(ab,ab[order(ab$stop),])
	   ab2<-which(Data$id==i) #row numbers in Data
	   start2<-vector(length=length(ab2))
	   start2[1]<-0
	   start2[2:length(ab2)]<-ab$stop[1:length(ab2)-1]
	   Data$start[ab2]<-start2
	} #end of for loop

	new.data <- data.frame(id=Data$id,start=Data$start,stop=Data$stop,st.stage=Data$st.stage,stage=Data$stage)
	res<-new.data
}

############################################################
################# Converting for Censoring #################
############################################################

Add.States <- function(tree){

   ##Adding censoring state to Nodes & Edges
   Nodes <- c("0",nodes(tree))
   Edges <- edgeL(tree)
   Edges[["0"]] <- list(edges=numeric(0))

   nt.states <- which(sapply(Edges, function(x) length(x$edges)>0)) #nonterminal states

  for(stage in nt.states) {
    Edges[[stage]]$edges <- c("0",Edges[[stage]]$edges)
  }

 ##tree for censored data
 tree0 <- new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

## Adding "Left Truncated" State
Nodes<- c("LT",nodes(tree0))
Edges[["LT"]] <- list(edges=nodes(tree)[nodes(tree)%in%names(nt.states)])
  nt.states.LT <- which(sapply(Edges, function(x) length(x$edges)>0)) #nonterminal states


treeLT <-new("graphNEL",nodes=Nodes,edgeL=Edges,edgemode="directed")

 list(tree0=tree0,nt.states=nt.states,nt.states.LT=nt.states.LT,treeLT=treeLT)

}
############################################################
############# Adding Dummy "LT" obs to Data set ############
############################################################

LT.Data <- function(Data){
 
  ## NOTE: Below assumes all the variables in Data have the names 'id', 'start', 'stop', etc., so make sure you set that up IN YOUR FUNCTION
  ## Best way to check that is to make an example with different column names, and pass that to your function
 
   Data <- Data[order(Data$id), ]  ## make sure id's line up below
   ids <- unique(Data$id)
   stop.time <- with(Data, tapply(start, id,min))
#   dummy <- data.frame(id = ids, start = -1, stop = stop.time, st.stage="LT", stage=1) #dummy initial stage
##1/9/11 -- changed "stage" argument in dummy so that it is the starting stage, added "enter.st" to do that
   enter.st<- with(Data, tapply(st.stage, id,min))
   dummy <- data.frame(id = ids, start = -1, stop = stop.time, st.stage="LT", stage=enter.st) #dummy initial stage

   Data <- rbind(Data, dummy) 
   Data <- with(Data, Data[order(id,stop), ])

  return(Data=Data)

}



############################################################
################ Counting Process & At Risk ################
####### Also, components for D-S estimator if needed #######
############################################################

## Assuming stage 0 is censored and stage 1 is the initial state

CP <- function(tree,tree0,Data,nt.states){
   #tree0=tree for uncens, tree0=tree0 for cens
   #default is Datta-Satten estimator FALSE

#   times <- sort(unique(Data$stop[!Data$stop==0])) #unique non-0 times in entire data set
##Change to work with the KM Dataset bmt since one stop time =0
   times <- sort(unique(Data$stop)) #unique non-0 times in entire data set

   lng <- sapply(edges(tree0)[nodes(tree0)%in%names(nt.states)], length)
   ds <- paste("dN", rep(nodes(tree0)[nodes(tree0)%in%names(nt.states)], lng),
          unlist(edges(tree0)[nodes(tree0)%in%names(nt.states)])) #used for dNs
   ys <- paste("y",unlist(nodes(tree0)))  #used for Ys

   ## index of obs in each stage/state/node
   indy <- vector(length=length(ys),mode="list")
   names(indy) <- ys

   #indicator for counting set, includes trans for cens, cens occurs in nonterminal states
   indD <- vector(length=length(ds),mode="list")

   # matrix of # of transitions, initialize to zeros
   dNs <- matrix(0, nrow=length(times), ncol=length(ds))

   # matrix of total # of transitions from a state, initialize to zeros
   sum.dNs <- matrix(0, nrow=length(times), ncol=length(nt.states))

   # matrix of at-risk sets for each stage at each time
   Ys <- matrix(NA, nrow=length(times), ncol=length(ys))

   #names of rows/columns for vectors/matrices
   rownames(dNs) <- rownames(sum.dNs) <- rownames(Ys) <- times
   names(indD) <- colnames(dNs) <- ds
   colnames(Ys) <- ys
   colnames(sum.dNs) <- paste("dN",names(nt.states),".")

n.vec<-vector(length=length(nodes(tree0)))

   for(i in nodes(tree0)){ #loop through nodes

     nam <- strsplit(names(indy)," ")
     idx <- which(sapply(nam, function(x) x[2]==i))
      indy[[ys[idx]]] <- which(Data$stage==i)  #list of those in the different stages/nodes

     if(length(inEdges(tree0)[[i]])==0) next #skip computations for those in initial state or have no length

     ld <- length(inEdges(tree0)[[i]]) ### length of the number of transitions per stage

     for(j in 1:ld){ #Fill-in no. transitioning between stages at each time

        nam2 <- paste("dN", inEdges(tree0)[[i]][j], i) # assigns name for different transitions
	## use the following if looking at previous time's ending state
##        indD[[nam2]] <- indy[[idx]][Data$stage[(indy[[idx]]-1)]==inEdges(tree0)[[i]][j]]
	##use the following if looking at the starting state for current time
        indD[[nam2]] <- indy[[idx]][Data$st.stage[indy[[idx]]]==inEdges(tree0)[[i]][j]]
        # create table with times, converts time to character & NAs to 0s
        tmp.tab <- table(Data$stop[indD[[nam2]]])  # no. trans at each trans time
        dNs[names(tmp.tab),nam2] <- tmp.tab
     }
   } #end of outer loop for dNs

## to figure out distribution of states for individuals at their 
## initial observation time

res <- by(Data, Data$id, function(x) x$st.stage[which.min(x$stop)])
res <- factor(res, levels=nodes(tree), labels=nodes(tree))
start.probs <- table(res)/length(res)

## Starting stage at TIME 0
## Two ways to have user explicitly give these

## 1. Give vector of length N, where each element in the vector
##    has the corresponding starting state for that individual
##    (equivalent to 'res' above)

## 2. Give a vector of length equal to the number of states, where
##    each element in the vector gives the probability of starting 
##    in that state (equivalent to 'start.probs' above)
##    IF want to be really flexible, can allow user to choose

   ### starting at risk computations ###
   for(i in nodes(tree0)){ #loop through nodes to find Ys


     n <- length(which(res==i))

     nam <- strsplit(names(indy)," ")
     idx <- which(sapply(nam, function(x) x[2]==i))

      if (length(inEdges(tree0)[[i]])>0)
        into.node <- paste("dN", inEdges(tree0)[[i]], i) else into.node <- NULL
      if (length(edges(tree0)[[i]])>0)
        from.node <- paste("dN", i, edges(tree0)[[i]]) else from.node <- NULL
#1/25/11 - remove initial reference and update Ys, correcting so initial state other than 1

  initial <- which(sapply(inEdges(tree0), function(x) !length(x)>0)) #nonterminal states
  transient <- which(sapply(edges(tree0),function(x) length(x)>0) & sapply(inEdges(tree0),function(x) length(x)>0))

	#to see initial state:  nodes(tree0)[initial]
	#to see transient states:  nodes(tree0)[transient]

	if (i==names(initial)){
        Ys[,idx] <- c(n, n + cumsum(rowSums(dNs[,into.node, drop=FALSE])) - cumsum(rowSums(dNs[,from.node, drop=FALSE])))[-(nrow(Ys)+1)]
	} else if(i==names(transient) && !n==0){
        Ys[,idx] <- c(n, n + cumsum(rowSums(dNs[,into.node, drop=FALSE])) - cumsum(rowSums(dNs[,from.node, drop=FALSE])))[-(nrow(Ys)+1)]
	} else Ys[,idx] <- c(0, cumsum(rowSums(dNs[,into.node, drop=FALSE])) - cumsum(rowSums(dNs[,from.node, drop=FALSE])))[-(nrow(Ys)+1)]

   } #end of loop for Ys

    ## Counting transitions from different stages (ie: stage sums)
    sum.dNs <- matrix(nrow=nrow(dNs),ncol=length(nt.states))
    rownames(sum.dNs) <- rownames(dNs) #
    colnames(sum.dNs) <- paste("dN",names(nt.states),".")
    a <- strsplit(colnames(sum.dNs), " ")
    a2 <- strsplit(colnames(dNs), " ")
    uni <- unique(sapply(a,function(x) x[2]))# gives the unique stages exiting

    for(i in uni){ #calculating the dNi.s
      b <- which(sapply(a,function(x) x[2]==i))
      b2 <- which(sapply(a2,function(x) x[2]==i))
	 sum.dNs[,b] <- rowSums(dNs[,b2])
    } #end of for loop for calculating dNi.s

    list(dNs=dNs,Ys=Ys,sum.dNs=sum.dNs,res=res,start.probs=start.probs)

} #end of function

############################################################
################ Datta-Satten Estimation ###################
############################################################

DS <- function(LT="LT",nt.states,dNs,sum.dNs,Ys,Cens="0",cens.type){
   ## Calculating dNs, sum.dNs, and Y from D-S(2001) paper
   ## Dividing dNs*, sum.dNs*, & Y* by K to get dNs^, sum.dNs^, & Ys^
   ## Make sure nt.states is from the non-LT

    ## Removing "LT" state and non-event times
       res <- strsplit(colnames(dNs), " ") #string splits names
       res2 <- strsplit(colnames(Ys)," ") #string split names of Ys
       res3 <- strsplit(colnames(sum.dNs)," ") #string splits names of dNs

       DS.col.idx <- which(sapply(res, function(x) x[3]==Cens)) # looks at censored columns,needed for D-S est
       DS2.col.idx <- which(sapply(res2, function(x) x[2]%in%names(nt.states))) # looks at censored columns,needed for D-S est
       DS3.col.idx <- which(sapply(res3, function(x) x[2]%in%names(nt.states))) # looks at censored columns,needed for D-S est
## This needs to be the nt.states from tree0 if we want to leave x[2]!=LT in the statement
## If we don't want that, we need nt.states from treeLT
   if(cens.type=="ind"){ ## for INDEPENDENT censoring

	K <- vector(length=nrow(dNs))
	dN0 <- rowSums(dNs[,DS.col.idx])
	Y0 <- rowSums(Ys[,DS2.col.idx]) #those at risk of being censored
	N.Y <- ifelse(dN0/Y0=="NaN",0,dN0/Y0)
	colnames(N.Y) <- NULL
	H.t <- cumsum(N.Y) #calculating the hazard
	k <- exp(-H.t)
	K <- c(1, k[-length(k)])

	dNs.K <- dNs/K  #D-S dNs
	Ys.K <- Ys/K  #D-S Ys
        sum.dNs.K <- sum.dNs/K
    } #end of ind censoring if

## Dependent censoring

    if(cens.type=="dep"){

      dN0 <- dNs[,DS.col.idx]
      Y0 <- Ys[,DS2.col.idx] #those at risk of being censored

      N.Y <- ifelse(dN0/Y0=="NaN",0,dN0/Y0)
	colnames(N.Y) <- paste(colnames(dN0),"/",colnames(Y0))

      H.t <- apply(N.Y, 2, function(x) cumsum(x))
      K <- exp(-H.t)
      ## K <- apply(k, 2, function(x) c(1, x[-length(x)]))  ## maybe don't need

	## This gives the dNs that need to be divided by Ys (ie: transitions b/t states in tree)

      ab <- which(sapply(res,function(x) x[2]%in%nt.states))
      ac <- which(sapply(res3,function(x) x[2]%in%nt.states)) #those in sum.dNs
	dNs.K <-dNs; Ys.K <- Ys; sum.dNs.K <- sum.dNs
	for(i in names(nt.states)){
		K.idx <- which(sapply(strsplit(colnames(N.Y)," "),function(x) x[2]==i))
		dN.idx <- which(sapply(res,function(x) x[2]==i))
                sum.dNs.idx <- which(sapply(res3,function(x) x[2]==i))
		Ys.idx <- which(sapply(res2,function(x) x[2]==i))
		dNs.K[,dN.idx] <- dNs[,dN.idx]/K[,K.idx]
                sum.dNs.K[,sum.dNs.idx] <- sum.dNs[,sum.dNs.idx]/K[,K.idx]
		Ys.K[,Ys.idx] <- Ys[,Ys.idx]/K[,K.idx]
	}
    } #end of if dependent censoring

      res <- list(dNs.K=dNs.K,Ys.K=Ys.K,sum.dNs.K=sum.dNs.K)
      return(res)

} ## end of D-S function

############################################################
############ Reducing dNs & Ys to event times ##############
############################################################

## By creating this function, I can simplify the CP function
## will need to update CP, remove the DS arguement in fn defn OK.
## will need to update call of CP in msSurv
## will need to add call of Red in msSurv function
## will need to update reference to dNs, Ys, sum.dNs in other functions
## mainly the call of functions stateoccfn, variances, et, & res
## will also need to add a separate function for D-S estimator ...
## remove DS arguement from CP

Red <- function(tree,dNs,Ys,sum.dNs,dNs.K,Ys.K,sum.dNs.K){

  ## tree is original tree currently inputted by user
  ## dNs, sum.dNS, & Ys come from CP function
  ## K comes from DS

   ##reducing dNs & Ys to just event times & noncens/nontruncated states
    res <- strsplit(colnames(dNs), " ") #string splits names
    col.idx <- which(sapply(res, function(x) x[2]%in%nodes(tree) & x[3]%in%nodes(tree))) # looks at noncensored columns
    row.idx <- which(apply(dNs[,col.idx], 1, function(x) any(x>0))) #identifies times where transitions occur
    dNs.et <- dNs[row.idx,col.idx] ## reduces dNs

    res2 <- strsplit(colnames(Ys)," ") #string split names of Ys
    nt.states.f <- which(sapply(edges(tree), function(x) length(x)>0)) #nonterminal states
    col2.idx <- which(sapply(res2,function(x) x[2]%in%names(nt.states.f))) #ids nonterminal columns
    Ys.et <- Ys[row.idx,col2.idx] ## reduces Ys

    col3.idx <- which(sapply(strsplit(colnames(sum.dNs)," "),function(x) x[2]%in%nodes(tree)))
    sum.dNs.et <- sum.dNs[row.idx,col3.idx]

    dNs.K.et <- dNs.K[row.idx,col.idx]
    Ys.K.et <- Ys.K[row.idx,col2.idx]
    sum.dNs.K.et <- sum.dNs.K[row.idx,col3.idx]

    ans <- list(dNs=dNs.et,Ys=Ys.et,sum.dNs=sum.dNs.et,dNs.K=dNs.K.et,Ys.K=Ys.K.et,sum.dNs.K=sum.dNs.K.et)
    return(ans)

}


############################################################
############# State Occupation Probabilities ###############
############################################################

stocc <- function(ns,tree,dNs.et,Ys.et,start.probs){
 ##currently ns is defined in main function
 ##tree needs to be uncensored tree

   cum.tm <- diag(ns)
   colnames(cum.tm) <- rownames(cum.tm) <- nodes(tree)

   ps <- matrix(NA, nrow=nrow(dNs.et), ncol=length(nodes(tree)))
      rownames(ps) <- rownames(dNs.et); colnames(ps) <- paste("p",nodes(tree))
all.dA <- all.I_dA <- all.ajs <- array(dim=c(ns,ns,nrow(dNs.et)),dimnames=list(rows=nodes(tree),cols=nodes(tree),dim=rownames(dNs.et)))
## The above would rename dim names for all.ajs (dim of array=event times)

#   all.dA <- all.I_dA <- all.ajs <- array(dim=c(ns,ns,nrow(dNs.et)))
#      colnames(all.ajs) <- colnames(ps)

   for(i in 1:nrow(dNs.et)){ ##loop through times

	I_dA <- diag(ns) #creates trans matrix for current time
	dA <- matrix(0,nrow=ns,ncol=ns)
	colnames(I_dA) <- rownames(I_dA) <- colnames(dA) <- rownames(dA) <- nodes(tree)

	idx <- which(dNs.et[i,]>0)  ## transition time i
	t.nam <- colnames(dNs.et)[idx] ## gets names of transitions (ie:  dN##)
	tmp <- strsplit(t.nam," ") ## splits title of dN##
	start <- sapply(tmp, function(x) x[2])
	end <- sapply(tmp, function(x) x[3])  ## pulls start & stop states as character strings
	idxs <- matrix(as.numeric(c(start, end)), ncol=2)
	idxs2 <- matrix(as.numeric(c(start, start)), ncol=2)

	dA[idxs] <- dNs.et[i,idx]/Ys.et[i,paste("y",start)]
      if(length(idx)==1) dA[start,start] <- -dNs.et[i,idx]/Ys.et[i,paste("y",start)]
		else dA[idxs2] <- -rowSums(dA[start, ])
	## 11/10/10 - updated this part for will compute if more than 1 type of transition at different event times

	I_dA <- I_dA + dA #I+dA matrix

	all.dA[,,i] <- dA  #stores all dA matrices
	all.I_dA[,,i] <- I_dA ## array for storing all tran matrices

	cum.tm <- cum.tm %*% I_dA # Multiply cur.tm and cum.tm to get matrix for current time
	all.ajs[,,i] <- cum.tm #A-J estimates, stored in array

#	ps[i,] <- all.ajs[1,,i] #just the state occupation probabilities
##1/26/11--added multiplication to calculate state occ for models with starting states other than 1
###start.probs=pk0
	ps[i,] <- start.probs%*%all.ajs[,,i] #just the state occupation probabilities

    } #end of loop

    list(ps=ps,all.ajs=all.ajs,all.I_dA=all.I_dA)
} #end of function

############################################################
############## State Entry/Exit Distributions ##############
############################################################

Dist <- function(ps,ns,tree){
 #ps from stocc function
 #tree needs to be uncensored tree

    initial <- which(sapply(inEdges(tree), function(x) !length(x)>0)) #initial states, no Fs
    terminal <- which(sapply(edges(tree), function(x) !length(x)>0)) #terminal states, no Gs

    Fs <- matrix(0, nrow=nrow(ps), ncol=ns) #entry distn
      rownames(Fs) <- rownames(ps)
      colnames(Fs) <- paste("F",nodes(tree))

    Gs <- matrix(0, nrow=nrow(ps), ncol=ns) #exit distn
      rownames(Gs) <- rownames(ps)
      colnames(Gs) <- paste("G",nodes(tree))

    for(i in 1:ns){#looping through nodes
	 node <- nodes(tree)[i]
	 later.stages <- names(acc(tree, node)[[1]])
	 stages <- c(node, later.stages)

       f.numer <- rowSums(ps[,paste("p", stages),drop=FALSE])
       Fs[,i] <- f.numer/f.numer[length(f.numer)]

	 if(length(stages)==1) next

	 g.numer <- rowSums(ps[,paste("p", later.stages),drop=FALSE])
	 Gs[,i] <- g.numer/g.numer[length(g.numer)]

     } #end of for loop

	Fr <- strsplit(colnames(Fs)," ")
	Fs.idx <- which(sapply(Fr,function(x) x[2]%in%names(initial)))
	Fs[,Fs.idx]<-NA

	Gr <- strsplit(colnames(Gs)," ")
	Gs.idx <- which(sapply(Gr,function(x) x[2]%in%names(terminal)))
	Gs[,Gs.idx]<-NA

    list(Fs=Fs,Gs=Gs)
} #end of function

############################################################
######################### Variance #########################
############################################################

var.fn <- function(tree,ns,nt.states,dNs.et,Ys.et,sum.dNs,all.ajs,all.I_dA,ps){

    #elements needed for computation
    varcov <- array(0, dim = c(ns^2,ns^2,nrow(dNs.et)))
      colnames(varcov) <- rownames(varcov) <- paste(rep(nodes(tree),ns),sort(rep(nodes(tree),ns)))
    bl.Id <- diag(1,(ns)^2) #Ident matrix for Kronecker product
    tm <- matrix(0,nrow=ns,ncol=ns) #tmp matrix to col var est
    res.array <- array(0,dim(tm)^2)
      colnames(res.array) <- rownames(res.array) <- paste(rep(nodes(tree),ns),sort(rep(nodes(tree),ns)))
    out <- array(0, dim=c(dim(all.I_dA)[c(1, 2)]^2,nrow(dNs.et)))
      colnames(out) <- rownames(out) <- paste(rep(nodes(tree),ns),sort(rep(nodes(tree),ns)))
    Id <- diag(1,ns)
    cov.p <- matrix(0,nrow=nrow(dNs.et),ncol=ns) #matrix for var matrix of state occup prob
      colnames(cov.p) <- paste("Var", "p",nodes(tree))
      rownames(cov.p) <- rownames(ps)
    v.p <- matrix(0,ns,ns) #variance of p-hat(0), needed to alter if start time not 0

    for(i in 1:nrow(dNs.et)){ ##loop through times


        #VARIANCE OF A-J (TRANS PROB MATRIX P(0,t))
	for(outer in 1:length(nt.states)){ #loop on the blocks (g)

	  tm <- matrix(0,nrow=ns,ncol=ns) #resets tm to 0 for next loop through outer
	  for(j in 1:ns){ #loop in the blocks

	     for(k in j:ns){  ## this just fills in upper diagonal matrix
                           ## use symmetry to fill in rest (or not bother)

	        if(Ys.et[i,outer]==0){  ## if Y_g = 0 the covariance = 0
        	  	tm[j,k] <- 0
	        	next
       	  } #end of if

       	  if (j == outer & k == outer) {  ## 3rd formula
			tm[j,k] <- (Ys.et[i,outer]-sum.dNs[i,outer])*sum.dNs[i,outer]/Ys.et[i,outer]^3
	        }  else if (j == outer & k != outer) {  ## 2nd formula
	            name <- paste("dN", outer, k)
			if (!name%in%colnames(dNs.et)) next
			tm[j,k] <- -(Ys.et[i,outer]-sum.dNs[i,outer])*dNs.et[i,name]/Ys.et[i,outer]^3
		  } else if (j != outer & k == outer) {  ## 2nd formula pt 2, for recurrent
	            name <- paste("dN", outer, j)
			if (!name%in%colnames(dNs.et)) next
			tm[j,k] <- -(Ys.et[i,outer]-sum.dNs[i,outer])*dNs.et[i,name]/Ys.et[i,outer]^3
		  } else { ## 1st formula
			namek <- paste("dN", outer, k)
			namej <- paste("dN", outer, j)
			if (!(namej%in%colnames(dNs.et) & namek%in%colnames(dNs.et))) next
			tm[j,k] <- (ifelse(j==k, 1, 0)*Ys.et[i,outer]-dNs.et[i,namej])*dNs.et[i,namek]/Ys.et[i,outer]^3
              } #end of if/else statements
	     } ## end of k loop
	  } ## end of j loop

      tm[lower.tri(tm)] <- t(tm)[lower.tri(tm)]

      res.array[(seq(1, ns*(ns-1)+1, by=ns)+outer-1), (seq(1, ns*(ns-1)+1, by=ns)+outer-1)] <- tm

	}#end of outer loop

	varcov[,,i] <- res.array #array holding var-cov matrix for trans prob

	if(i==1) out[, , i] <- bl.Id%*% varcov[,,i] %*% bl.Id
	  else out[, , i] <- (t(all.I_dA[, , i]) %x% Id) %*% out[, , i-1] %*%((all.I_dA[, , i]) %x% Id) +
	      (Id %x% all.ajs[, , i-1]) %*% varcov[,,i]  %*% (Id%x% t(all.ajs[, , i-1]))
	      #note:  all.ajs[,,i-1] corresponds to P(0,t-)
     		#out is the var/cov est of P(0,t)

## p.t <- as.vector(c(1,rep(0,ns-1)))
## can be used for time s=0 with v1,vv1,var1.mat, & var.pkj01 <- var1.mat
## if s!=0, everything in place to say if (i=1) p.t as above
## else p.t <- ps[i-1,] and then everything should be generalized

      ## calculating the variance of state occupation prob
	for (j in nodes(tree)){ #loop through states

         st.nam <- paste("1",j)
##         v1 <- strsplit(colnames(out), "")
##         vv1 <- which(sapply(v1,function(x) x[3]==j))
##         var1.mat <- out[vv1,vv1,i]  #creates a matrix of var for state j (2nd # in title is j)
	   part1 <- var.pkj0t <- out[st.nam,st.nam,i]


##	   part1 <- t(p.t)%*%var.pkj0t%*%p.t

	   res3 <- strsplit(colnames(ps)," ")
	   col.idx3 <- which(sapply(res3, function(x) x[2]== j)) # looks at state transitioned to
	   b.t <- all.ajs[,col.idx3,i] #creating vector of col for current state from trans prob

	   part2 <- t(b.t)%*%v.p%*%b.t #should be 0 when P(0,t)
	    ##right now forced to be 0 by the way v.p defined outside of time loop

	   res.varp <- part1+part2 #calculating var/cov matrix for time i, state j

	   cov.p[i,as.numeric(j)] <- res.varp #storing var/cov calc for time i, state j

	} #closes states loop
    } ## end of time loop

    list(out=out,varcov=varcov,cov.p=cov.p)

}#end of function


#########################################################
################ BS Variance for Dep Cens ###############
#########################################################

BS.var <- function(Data,tree,ns,et,cens.type,B,LT,start.states){

	n <- length(unique(Data$id)) # sample size
	ids <- unique(Data$id)

	bs.est <- array(dim=c(length(nodes(tree)),length(nodes(tree)),length(et),B),
				dimnames=list(rows=nodes(tree),cols=nodes(tree),dim=et)) # storage for bootstrap estimates of transition probability matrices 
	bs.ps <- array(dim=c(length(et),ns,B))
        rownames(bs.ps) <- et
        colnames(bs.ps) <- paste("p",nodes(tree))

	## For entry / exit distributions
	## Return these as well
	bs.Fs <- bs.ps; bs.Gs <- bs.ps #storage for BS Fs/Gs
		colnames(bs.Fs) <- paste("F",nodes(tree))
		colnames(bs.Gs) <- paste("G",nodes(tree))
      initial <- which(sapply(inEdges(tree), function(x) !length(x)>0)) #initial states, no Fs
      terminal <- which(sapply(edges(tree), function(x) !length(x)>0)) #terminal states, no Gs

      bs.cov.p <- matrix(0,nrow=length(et),ncol=ns) #matrix for var matrix of state occup prob
      colnames(bs.cov.p) <- paste("Var", "p",nodes(tree))
      rownames(bs.cov.p) <- et

      res.array <- array(0,dim=c(ns^2,ns^2,length(et)),dimnames=list(rows=paste(rep(nodes(tree),ns),sort(rep(nodes(tree),ns))),cols=paste(rep(nodes(tree),ns),sort(rep(nodes(tree),ns))),dim=et))

	for(b in 1:B){ #randomly selects the indices

	## Find the bootstrap sample, pull bs sample info from original data & put into data set
		bs=sample(ids, n, replace=TRUE)
		bs=factor(bs, levels=ids)
		bs.tab=data.frame(table(bs)) ##table with the frequencies
		Data.bs=merge(Data, bs.tab, by.x="id", by.y="bs") ## merging original data with bs freq table
		bs.id=unlist(apply(Data.bs[Data.bs$Freq>0,], 1, function(x) paste(x["id"], 1:x["Freq"], sep="."))) ## creating bs id
		idx=rep(1:nrow(Data.bs), Data.bs$Freq) ##indexing the bs sample
		Data.bs=Data.bs[idx,] ##creating a bs dataset
		Data.bs.originalID=Data.bs$id
		Data.bs$id=bs.id #changing id column to bs.id to use functions
		Data.bs=Data.bs[order(Data.bs$stop),] #ordered bs dataset

                ## Calling functions using bs dataset
		Cens <- Add.States(tree)
                ##1/26/11 - Added the LT case
		  if(LT) {
			Data.bs = LT.Data(Data.bs)
			cp <- CP(tree,Cens$treeLT,Data.bs,Cens$nt.states.LT)
			## here calculate start.probs based on vector of starting states
			## start.probs <- table(XXX)
			res <- factor(start.states, levels=nodes(tree), labels=nodes(tree))
			start.probs <- table(res)/length(res)
		   }

		   if(!LT) {
			cp <- CP(tree,Cens$tree0,Data.bs,Cens$nt.states)
			start.probs <- cp$start.probs
		  }


#		cp <- CP(tree,Cens$tree0,Data.bs,Cens$nt.states)
		ds.est<-DS(LT="LT",Cens$nt.states,cp$dNs,cp$sum.dNs,cp$Ys,Cens="0",cens.type)
		cp.red <- Red(tree,cp$dNs,cp$Ys,cp$sum.dNs,ds.est$dNs.K,ds.est$Ys.K,ds.est$sum.dNs.K)
		stateoccfn <- stocc(ns,tree,cp.red$dNs.K,cp.red$Ys.K,start.probs) ##ns is output of msSurv object

		idx <- which(dimnames(bs.est)[[3]] %in% dimnames(stateoccfn$all.I_dA)[[3]])
		idx2 <- which(!(dimnames(bs.est)[[3]] %in% dimnames(stateoccfn$all.I_dA)[[3]]))
		bs.IA <- bs.est		

		bs.IA[,,idx,b] <- stateoccfn$all.I_dA  ## or wherever you have stored the transition probabilities ... 
		bs.IA[,,idx2,b] <- diag(ns)

		bs.est[,,1,b] <- bs.IA[,,1,b]
		bs.ps[1,,b] <- start.probs%*%bs.est[,,1,b]

		for(j in 2:length(et)){
			bs.est[,,j,b] <- bs.est[,,j-1,b] %*% bs.IA[,,j,b]
			bs.ps[j,,b] <- start.probs%*%bs.est[,,j,b]
		} ## end of j for loop 

	  ## Entry / Exit variance as well 
	  for(i in 1:ns){#looping through nodes
		 node <- nodes(tree)[i]
		 later.stages <- names(acc(tree, node)[[1]])
		 stages <- c(node, later.stages)

	       bs.f.numer <- rowSums(bs.ps[,paste("p", stages),b,drop=FALSE])
		if(sum(bs.f.numer)==0)  bs.Fs[,i,b]<-0
	          else bs.Fs[,i,b] <- bs.f.numer/bs.f.numer[length(bs.f.numer)]

	
		 if(length(stages)==1) next

		 bs.g.numer <- rowSums(bs.ps[,paste("p", later.stages),b,drop=FALSE])
		if(sum(bs.g.numer)==0)  bs.Gs[,i,b]<-0
		 else bs.Gs[,i,b] <- bs.g.numer/bs.g.numer[length(bs.g.numer)]

	  } #end of for loop

	} ## end of bs loop

	Fs.var <- apply(bs.Fs,c(1,2),var)
	Fs.var[,initial]<-NA #setting the initial state variances = NA
	Gs.var <- apply(bs.Gs,c(1,2),var)
	Gs.var[,terminal] <- NA #setting distn for terminal states to NA since don't exist

	bs.var <- apply(bs.est, c(1,2,3), var) ## calculate var of each transition matrix at each time point
#	bs.cov.p <- t(bs.var[1,,]) #if all start in state 1
#	bs.p <- start.probs%*%bs.est[,,,b]
##1/26/11 - replaced BS cov(p) computation to account for starting state other than 1
	bs.cov.p <- apply(bs.ps,c(1,2),var)
        colnames(bs.cov.p) <- paste("Var", "p",nodes(tree))
        rownames(bs.cov.p) <- et

	bs.cov <- array(dim=c(ns^2,ns^2,length(et)),dimnames=list(rows=paste(rep(1:ns,ns), rep(1:ns, each=ns)),cols=paste(rep(1:ns,ns), rep(1:ns, each=ns)),dim=et))
	for(i in 1:length(et)){
		bs.est2 <- matrix(bs.est[,,i,],nrow=B, ncol=ns^2, byrow=TRUE)
		bs.cov[,,i] <- cov(bs.est2)
	}  ##this for loop creates a B x (# of states)^2 x (# of event times)

	list(out=bs.cov,cov.p=bs.cov.p, Fs.var=Fs.var,Gs.var=Gs.var)

} ## end of function


##1/14/11 -- Adding Var for Fs/Gs
#########################################################
############ BS Variance for Entry/Exit Distn ###########
#########################################################

Dist.BS.var <- function(Data,tree,ns,et,dNs.K,cens.type,B,LT,start.probs){

#	B <- 200 # number of replicates
	n <- length(unique(Data$id)) # sample size
	ids <- unique(Data$id)
        initial <- which(sapply(inEdges(tree), function(x) !length(x)>0)) #initial states, no Fs

        terminal <- which(sapply(edges(tree), function(x) !length(x)>0)) #terminal states, no Gs


	bs.est <- array(dim=c(length(nodes(tree)),length(nodes(tree)),length(et),B),
				dimnames=list(rows=nodes(tree),cols=nodes(tree),dim=rownames(dNs.K))) # storage for bootstrap estimates of transition probability matrices 

	bs.ps <- array(dim=c(length(et),ns,B))
        rownames(bs.ps) <- et
        colnames(bs.ps) <- paste("p",nodes(tree))

	bs.Fs <- bs.ps; bs.Gs <- bs.ps #storage for BS Fs/Gs
		colnames(bs.Fs) <- paste("F",nodes(tree))
		colnames(bs.Gs) <- paste("G",nodes(tree))


	for(b in 1:B){ #randomly selects the indices

	## Find the bootstrap sample, pull bs sample info from original data & put into data set
		bs=sample(ids, n, replace=TRUE)
		bs=factor(bs, levels=ids)
		bs.tab=data.frame(table(bs)) ##table with the frequencies
		Data.bs=merge(Data, bs.tab, by.x="id", by.y="bs") ## merging original data with bs freq table
		bs.id=unlist(apply(Data.bs[Data.bs$Freq>0,], 1, function(x) paste(x["id"], 1:x["Freq"], sep="."))) ## creating bs id
		idx=rep(1:nrow(Data.bs), Data.bs$Freq) ##indexing the bs sample
		Data.bs=Data.bs[idx,] ##creating a bs dataset
		Data.bs.originalID=Data.bs$id
		Data.bs$id=bs.id #changing id column to bs.id to use functions
		Data.bs=Data.bs[order(Data.bs$stop),] #ordered bs dataset

	## Calling functions using bs dataset
		Cens <- Add.States(tree)
#1/26/11 - Added the LT case
		  if(LT) {
			Data.bs = LT.Data(Data.bs)
			cp <- CP(tree,Cens$treeLT,Data.bs,Cens$nt.states.LT)
		   }

		   if(!LT) cp <- CP(tree,Cens$tree0,Data.bs,Cens$nt.states)

#		cp <- CP(Cens$tree0,Data.bs,Cens$nt.states)
		ds.est<-DS(LT="LT",Cens$nt.states,cp$dNs,cp$sum.dNs,cp$Ys,Cens="0",cens.type)
		cp.red <- Red(tree,cp$dNs,cp$Ys,cp$sum.dNs,ds.est$dNs.K,ds.est$Ys.K,ds.est$sum.dNs.K)
		stateoccfn <- stocc(ns,tree,cp.red$dNs.K,cp.red$Ys.K) ##ns is output of msSurv object

		idx <- which(dimnames(bs.est)[[3]] %in% dimnames(stateoccfn$all.I_dA)[[3]])
		idx2 <- which(!(dimnames(bs.est)[[3]] %in% dimnames(stateoccfn$all.I_dA)[[3]]))
		bs.IA <- bs.est		

		bs.IA[,,idx,b] <- stateoccfn$all.I_dA  ## or wherever you have stored the transition probabilities ... 
		bs.IA[,,idx2,b] <- diag(ns)

		bs.est[,,1,b] <- bs.IA[,,1,b]
		bs.ps[1,,b]<-start.probs%*%bs.est[,,1,b]

		for(j in 2:length(et)){
			bs.est[,,j,b] <- bs.est[,,j-1,b] %*% bs.IA[,,j,b]
			bs.ps[j,,b]<-start.probs%*%bs.est[,,j,b]
		} ## end of j for loop 

	  for(i in 1:ns){#looping through nodes
		 node <- nodes(tree)[i]
		 later.stages <- names(acc(tree, node)[[1]])
		 stages <- c(node, later.stages)

	       bs.f.numer <- rowSums(bs.ps[,paste("p", stages),b,drop=FALSE])
		if(sum(bs.f.numer)==0)  bs.Fs[,i,b]<-0
	          else bs.Fs[,i,b] <- bs.f.numer/bs.f.numer[length(bs.f.numer)]

	
		 if(length(stages)==1) next

		 bs.g.numer <- rowSums(bs.ps[,paste("p", later.stages),b,drop=FALSE])
		if(sum(bs.g.numer)==0)  bs.Gs[,i,b]<-0
		 else bs.Gs[,i,b] <- bs.g.numer/bs.g.numer[length(bs.g.numer)]

	  } #end of for loop

	} ## end of bs loop

	Fs.var <- apply(bs.Fs,c(1,2),var)
	Fs.var[,initial]<-NA #setting the initial state variances = NA
	Gs.var <- apply(bs.Gs,c(1,2),var)
	Gs.var[,terminal] <- NA #setting distn for terminal states to NA since don't exist

##Fs.var[,which(!(nodes(treeobj)%in%initial))]  #Fs.var for only the non-inital states

	list(Fs.var=Fs.var,Gs.var=Gs.var)

} ## end of function

#########################################################
######## CONFIDENCE INTERVALS for p(t) & P(s,t) #########
#########################################################

MSM.CIs <- function(x,ci.level=0.95,ci.trans="linear"){

#uses uncensored tree (tree)
# if state occup prob, cov.p
# if trans prob matrix, out
##we will also need to specify which type of CI wanted

#default ci.level is 0.95, default CI type (ie: ci.trans) is linear
    if(ci.level < 0 || ci.level > 1)
      stop("confidence level must be between 0 and 1")

    z.alpha <- qnorm(ci.level + (1 - ci.level) / 2)

    ci.trans <- match.arg(ci.trans,c("linear","log","cloglog","log-log"))

#CIs on state occup probability

#creating storage for CIs and labeling columns/rows
#each array corresponds to a different possible transition (each pos.trans)
#each row in the array corresponds to a time point

    CI.p <- array(0,dim=c(nrow(x@dNs),3,length(nodes(x@tree))),dimnames=list(rows=rownames(x@dNs),cols=c("est","lower limit","upper limit"),dim=paste("p",nodes(x@tree))))

#if i decide to add var to the list, be sure to add the name in cols= ... above
    for(i in 1:nrow(x@dNs)){ ##loop through times

#        time <- as.numeric(rownames(x@dNs))[i]  #possibly needed for report

      for (j in as.numeric(nodes(x@tree))){ #loop through states

        res.ci <- strsplit(colnames(x@ps), " ") #string splits names
        col.idx <- which(sapply(res.ci, function(x) x[2]== j)) # looks at state transitioned to
        res.ci2 <- strsplit(colnames(x@cov.p), " ") #string splits names
        col.idx2 <- which(sapply(res.ci2, function(x) x[3]== j)) # looks at state transitioned to

        CI.p[i,1,j]<- PE.p <- x@all.ajs[1,col.idx,i]
        var.p <- x@cov.p[i,col.idx2]

        switch(ci.trans[1],
		"linear" = {
	        CI.p[i,2,j] <- PE.p - z.alpha * sqrt(var.p)
      	        CI.p[i,3,j] <- PE.p + z.alpha * sqrt(var.p)},
		"log" = {
	        CI.p[i,2,j] <- exp(log(PE.p) - z.alpha * sqrt(var.p) / PE.p)
      	  	CI.p[i,3,j] <- exp(log(PE.p) + z.alpha * sqrt(var.p) / PE.p)},
        	"cloglog" = {
	        CI.p[i,2,j] <- 1 - (1 - PE.p)^(exp(z.alpha * (sqrt(var.p) / ((1 - PE.p) * log(1 - PE.p)))))
      	  	CI.p[i,3,j] <- 1 - (1 - PE.p)^(exp(-z.alpha * (sqrt(var.p) / ((1 - PE.p) * log(1 - PE.p)))))},
		"log-log" = {
	        CI.p[i,2,j] <- PE.p^(exp(-z.alpha * (sqrt(var.p) / (PE.p * log(PE.p)))))
      	  	CI.p[i,3,j] <- PE.p^(exp(z.alpha * (sqrt(var.p) / (PE.p * log(PE.p)))))})

 	  CI.p[i,2,j] <- pmax(CI.p[i,2,j],0)
	  CI.p[i,3,j] <- pmin(CI.p[i,3,j],1)
     } #end states loop
    } #end times loop for CI.p

## CIs on transition probability matrices##

    #creating storage for CIs and labeling columns/rows
    #each array corresponds to a different possible transition (each pos.trans)
    #pos.trans are 11,12,13,22,23
    #each row in the array corresponds to a time point

    CI.trans <- array(0,dim=c(nrow(x@dNs),4,length(x@pos.trans)),dimnames=list(rows=rownames(x@dNs),cols=c("est","lower limit","upper limit","var.tp"),dim=paste(x@pos.trans,"transition")))

    for(i in 1:nrow(x@dNs)){ ##loop through times

     for(j in 1:length(x@pos.trans)){ #loop through possible transitions

        idx <- as.numeric(unlist(strsplit(x@pos.trans[j], " ")))
        CI.trans[i,1,j] <- PE <- x@all.ajs[idx[1], idx[2] ,i]
        CI.trans[i,4,j] <- var <- x@out[x@pos.trans[j], x@pos.trans[j], i]


        switch(ci.trans[1],
		"linear" = {
	        CI.trans[i,2,j] <- PE - z.alpha * sqrt(var)
      	  CI.trans[i,3,j] <- PE + z.alpha * sqrt(var)},
		"log" = {
	        CI.trans[i,2,j] <- exp(log(PE) - z.alpha * sqrt(var) / PE)
      	  CI.trans[i,3,j] <- exp(log(PE) + z.alpha * sqrt(var) / PE)},
       	"cloglog" = {
	        CI.trans[i,2,j] <- 1 - (1 - PE)^(exp(z.alpha * (sqrt(var) / ((1 - PE) * log(1 - PE)))))
      	  CI.trans[i,3,j] <- 1 - (1 - PE)^(exp(-z.alpha * (sqrt(var) / ((1 - PE) * log(1 - PE)))))},
		"log-log" = {
	        CI.trans[i,2,j] <- PE^(exp(-z.alpha * (sqrt(var) / (PE * log(PE)))))
      	  CI.trans[i,3,j] <- PE^(exp(z.alpha * (sqrt(var) / (PE * log(PE)))))})

 	  CI.trans[i,2,j] <- pmax(CI.trans[i,2,j],0)
	  CI.trans[i,3,j] <- pmin(CI.trans[i,3,j],1)

     } #end j loop
    } #end times loop

    list(CI.p=CI.p,CI.trans=CI.trans)
} #end of function


#################################################################################
##################### CIS for distributions #####################################
#################################################################################

Dist.CIs <- function(x,ci.level=0.95,ci.trans="linear"){

##ests are the point estimates
##varest are the variance estimates

	z.alpha <- qnorm(ci.level + (1 - ci.level) / 2)

	ci.trans <- match.arg(ci.trans,c("linear","log","cloglog","log-log"))


	CI.Fs <- array(0,dim=c(nrow(x@Fs),3,length(nodes(x@tree))),dimnames=list(rows=rownames(x@Fs),cols=c("est","lower limit","upper limit"),dim=paste("F",nodes(x@tree))))
	CI.Gs <- array(0,dim=c(nrow(x@Gs),3,length(nodes(x@tree))),dimnames=list(rows=rownames(x@Gs),cols=c("est","lower limit","upper limit"),dim=paste("G",nodes(x@tree))))

    for(i in 1:nrow(x@Fs)){ ##loop through times

      for (j in as.numeric(nodes(x@tree))){ #loop through states

        res.ci.F <- strsplit(colnames(x@Fs), " ") #string splits names
        col.idx.F <- which(sapply(res.ci.F, function(x) x[2]== j)) # looks at state transitioned to
        res.ci2.F <- strsplit(colnames(x@Fs.var), " ") #string splits names
        col.idx2.F <- which(sapply(res.ci2.F, function(x) x[2]== j)) # looks at state transitioned to


        res.ci.G <- strsplit(colnames(x@Gs), " ") #string splits names
        col.idx.G <- which(sapply(res.ci.G, function(x) x[2]== j)) # looks at state transitioned to
        res.ci2.G <- strsplit(colnames(x@Gs.var), " ") #string splits names
        col.idx2.G <- which(sapply(res.ci2.G, function(x) x[2]== j)) # looks at state transitioned to

        CI.Fs[i,1,j]<- PE.F <- x@Fs[i,col.idx.F]
        varestF <- x@Fs.var[i,col.idx2.F]

        CI.Gs[i,1,j]<- PE.G <- x@Gs[i,col.idx.G]
        varestG <- x@Gs.var[i,col.idx2.G]


        switch(ci.trans[1],
		"linear" = {
	        CI.Fs[i,2,j] <- PE.F - z.alpha * sqrt(varestF)
      	        CI.Fs[i,3,j] <- PE.F + z.alpha * sqrt(varestF)
	        CI.Gs[i,2,j] <- PE.G - z.alpha * sqrt(varestG)
      	        CI.Gs[i,3,j] <- PE.G + z.alpha * sqrt(varestG)},
		"log" = {
	        CI.Fs[i,2,j] <- exp(log(PE.F) - z.alpha * sqrt(varestF) / PE.F)
      	  	CI.Fs[i,3,j] <- exp(log(PE.F) + z.alpha * sqrt(varestF) / PE.F)
		CI.Gs[i,2,j] <- exp(log(PE.G) - z.alpha * sqrt(varestG) / PE.G)
      	  	CI.Gs[i,3,j] <- exp(log(PE.G) + z.alpha * sqrt(varestG) / PE.G)},
        	"cloglog" = {
	        CI.Fs[i,2,j] <- 1 - (1 - PE.F)^(exp(z.alpha * (sqrt(varestF) / ((1 - PE.F) * log(1 - PE.F)))))
      	  	CI.Fs[i,3,j] <- 1 - (1 - PE.F)^(exp(-z.alpha * (sqrt(varestF) / ((1 - PE.F) * log(1 - PE.F)))))
		CI.Gs[i,2,j] <- 1 - (1 - PE.G)^(exp(z.alpha * (sqrt(varestG) / ((1 - PE.G) * log(1 - PE.G)))))
      	  	CI.Gs[i,3,j] <- 1 - (1 - PE.G)^(exp(-z.alpha * (sqrt(varestG) / ((1 - PE.G) * log(1 - PE.G)))))},
		"log-log" = {
	        CI.Fs[i,2,j] <- PE.F^(exp(-z.alpha * (sqrt(varestF) / (PE.F * log(PE.F)))))
      	  	CI.Fs[i,3,j] <- PE.F^(exp(z.alpha * (sqrt(varestF) / (PE.F * log(PE.F)))))
	        CI.Gs[i,2,j] <- PE.G^(exp(-z.alpha * (sqrt(varestG) / (PE.G * log(PE.G)))))
      	  	CI.Gs[i,3,j] <- PE.G^(exp(z.alpha * (sqrt(varestG) / (PE.G * log(PE.G)))))})

 	  CI.Fs[i,2,j] <- pmax(CI.Fs[i,2,j],0)
	  CI.Fs[i,3,j] <- pmin(CI.Fs[i,3,j],1)

 	  CI.Gs[i,2,j] <- pmax(CI.Gs[i,2,j],0)
	  CI.Gs[i,3,j] <- pmin(CI.Gs[i,3,j],1)

     } #end states loop


    } #end times loop for CI

    list(CI.Fs=CI.Fs,CI.Gs=CI.Gs)

} #end of function
