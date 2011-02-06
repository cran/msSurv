setGeneric("tree",function(object,...) standardGeneric("tree"))
setMethod("tree",signature(object="msSurv"),
          function(object) return(object@tree))

####################################################################
####################################################################
setMethod("print", signature(x="msSurv"),
   function(x, covar = FALSE, ee.distn=TRUE, ...) {
     ## ee.distn=TRUE  ENTRY/EXIT Distributions
     ## if (!inherits(x, "msSurv "))
     ## stop("'x' must be of class 'msSurv' ")

     ## nonterminal states
     transient <- as.character(which(sapply(edges(tree(x)), function(x) length(x) > 0)))
     ## absorbing states
     absorb <- as.character(which(sapply(edges(tree(x)), function(x) length(x) == 0)))
     trans <- strsplit(colnames(x@dNs),"")
     idx1 <- sapply(trans, function(x) x[4]); idx2 <- sapply(trans, function(x) x[6])
     trans <- paste(idx1,idx2)


       cat(paste("The specified multistate model has", length(transient), "transient state(s) and ", length(absorb), "absorbing state(s)\n\n", sep = " "))

       cat("Possible States in this Model:\n")
       print(nodes(x@tree))
       cat("\n")

       cat("Possible Transitions for this Model:\n")
       print(trans)
       cat("\n")

## start of state occupation prob info
 
       cat("State Occupation Information at time ", max(as.numeric(rownames(x@dNs))),": \n",sep="")
       cat("\n")

       cat(paste("Estimates of State Occupation Probabilities","\n", sep = ""))
       print(round(x@ps[nrow(x@ps),],digits=4))
       cat("\n")

       if (!is.null(x@cov.p) & covar == TRUE) {
           cat("Estimates of Covariance of State Occupation Probabilities","\n")
           print(round(x@cov.p[nrow(x@cov.p),],digits=4))
           cat("\n")
       }

       if (ee.distn) {
           cat("Estimates of State Entry Time Distribution","\n")
           print(round(x@Fs[nrow(x@Fs),],digits=4))
           cat("\n")

           cat("Estimates of State Exit Time Distribution","\n")
           print(round(x@Gs[nrow(x@Gs),],digits=4))
           cat("\n")
       }

       cat("\n")

## transition probability info
       cat("Transition Probability Information:", "\n", "\n")
       cat(paste("Estimate of P(",0,",",max(as.numeric(rownames(x@dNs))),")\n", sep = ""))

       ## set up currently to do P(0,max(x@et))
       print(round(x@all.ajs[,,dim(x@all.ajs)[3]],digits=4))
       cat("\n")

       if (!is.null(x@out) & covar == TRUE) {
           cat(paste("Estimate of cov(P(",0,",",max(as.numeric(rownames(x@dNs))),"))\n", sep = ""))
           print(round(x@out[, , dim(x@out)[3]],digits=4))
       }

       cat("\n")
       invisible()
})

####################################################################
####################################################################

setMethod("show", "msSurv",
   function(object) {

     ## nonterminal states
     transient <- as.character(which(sapply(edges(tree(object)), function(object) length(object) > 0)))
     ## absorbing states
     absorb <- as.character(which(sapply(edges(tree(object)), function(object) length(object) == 0)))
     trans <- strsplit(colnames(object@dNs),"")
     idx1 <- sapply(trans, function(x) x[4]); idx2 <- sapply(trans, function(x) x[6])
     trans <- paste(idx1,idx2)


       cat(paste("The specified multistate model has", length(transient), "transient state(s) and ", length(absorb), "absorbing state(s)\n\n", sep = " "))

       cat("Possible States in this Model:\n")
       print(nodes(object@tree))
       cat("\n")

       cat("Possible Transitions for this Model:\n")
       print(trans)
       cat("\n")

## start of state occupation prob info
       cat("State Occupation Information:", "\n", "\n")

       cat(paste("Estimates of State Occupation Probabilities","\n", sep = ""))
       print(object@ps)
       cat("\n")

           cat("Estimates of State Entry Time Distribution","\n")
           print(object@Fs)
           cat("\n")

           cat("Estimates of State Exit Time Distribution","\n")
           print(object@Gs)
           cat("\n")
       

       cat("\n")

## transition probability info
       cat("Transition Probability Information:", "\n", "\n")
       cat(paste("Estimate of P(",0,",",max(as.numeric(rownames(object@dNs))),")\n", sep = ""))

       ## set up currently to do P(0,max(object@et))
       print(object@all.ajs[,,dim(object@all.ajs)[3]])
       cat("\n")

       cat("\n")
       invisible()
})
####################################################################
####################################################################

setMethod("summary","msSurv",
   function(object, digits=3,all = FALSE, ci.fun = "linear", ci.level = 0.95, stateocc=TRUE, trans.pr=TRUE) {

       if (ci.level <= 0 | ci.level > 1) {
           stop ("confidence level must be between 0 and 1")
       }

       tmp <- MSM.CIs(object,ci.level=0.95)
       times <-object@et

       if(!all){
           dt <- quantile(times, probs = c(0,0.25,0.5,0.75,1))
           ind <- findInterval(dt,times) #indicator of which information to print
       }

## State Occupation Probability Section
       if(stateocc){
           cat("State Occupation Information:", "\n", "\n")

           if(all){

               for(i in seq(object@ns)){
                   cat(paste("State ", i, "\n"))
                   sop.sum <- data.frame(time=times,estimate=object@ps[,i],variance=object@cov.p[,i],lower.ci=tmp$CI.p[,2,i],upper.ci=tmp$CI.p[,3,i],Fs=object@Fs[,i],Gs=object@Gs[,i])
                   print(sop.sum,row.names=FALSE,digits=digits)
                   cat("\n")
               } #end of for statement

           } #end of if(all

           else{
               for(i in seq(object@ns)){

                   cat(paste("State ", i, "\n"))
                   sop.sum <- data.frame(time=times[ind],estimate=object@ps[ind,i],variance=object@cov.p[ind,i],lower.ci=tmp$CI.p[ind,2,i], upper.ci=tmp$CI.p[ind,3,i],entry.d=object@Fs[ind,i],exit.d=object@Gs[ind,i])

                   print(sop.sum,row.names=FALSE,digits=digits)
                   cat("\n")
               } #end of for statement

           } #end of else
       } #end of if(stateocc)


## Transition Probability Matrix Section
       if(trans.pr){

           cat("Transition Probability Information:", "\n", "\n")

           lt <- length(object@pos.trans)
           tts <- strsplit(object@pos.trans, split = " ")

           if(all){

               i=seq_along(object@pos.trans)[2]
               for (i in seq_along(object@pos.trans)) {

                   cat(paste("Transition", tts[[i]][1], "->", tts[[i]][2], "\n", sep = " "))

                   ## code to add number at risk and number transitions
                   ## print # events for transitions out of stage, else print # left
                   dns.name <- ifelse(tts[[i]][1] == tts[[i]][2], paste("dN", tts[[i]][1], ".", sep=" "), paste("dN", object@pos.trans[i], sep=" "))

                   ifelse(dns.name %in% colnames(object@dNs),n.event <- object@dNs[, dns.name],n.event <- object@sum.dNs[, dns.name])
                   ifelse(dns.name %in% colnames(object@dNs.K),n.event.K <- object@dNs.K[, dns.name],n.event.K <- object@sum.dNs.K[, dns.name])

                   ys.name <- paste("y", tts[[i]][1], sep=" ")
                   n.risk <- object@Ys[, ys.name]
                   n.risk.K <- object@Ys.K[,ys.name]

                   if (dns.name %in% colnames(object@dNs)) {
                       tp.sum <- data.frame(time=times,estimate=tmp$CI.trans[,1,i],variance=tmp$CI.trans[,4,i],lower.ci=tmp$CI.trans[,2,i],upper.ci=tmp$CI.trans[,3,i],n.risk = n.risk, n.event=n.event,n.risk.K=n.risk.K,n.event.K=n.event.K)
                   } else {
                       tp.sum <- data.frame(time=times,estimate=tmp$CI.trans[,1,i],variance=tmp$CI.trans[,4,i],lower.ci=tmp$CI.trans[,2,i],upper.ci=tmp$CI.trans[,3,i],n.risk = n.risk, n.remain=n.risk-n.event,n.risk.K=n.risk.K,n.remain.K=n.risk.K-n.event.K)
                   }


                   print(tp.sum, row.names = FALSE,digits=digits)
                   cat("\n")
               } #end of for loop
           } #end of if(all)

           else{

               for (i in seq_along(object@pos.trans)) {

                   cat(paste("Transition", tts[[i]][1], "->", tts[[i]][2], "\n", sep = " "))

                   ## code to add number at risk and number transitions
                   ## print # events for transitions out of stage, else print # left
                   dns.name <- ifelse(tts[[i]][1] == tts[[i]][2], paste("dN", tts[[i]][1], ".", sep=" "),paste("dN", object@pos.trans[i], sep=" "))

                   ifelse(dns.name %in% colnames(object@dNs),n.event <- object@dNs[, dns.name],n.event <- object@sum.dNs[, dns.name])

                   ifelse(dns.name %in% colnames(object@dNs.K),n.event.K <- object@dNs.K[, dns.name],n.event.K <- object@sum.dNs.K[, dns.name])

                   ys.name <- paste("y", tts[[i]][1], sep=" ")
                   n.risk <- object@Ys[, ys.name]
                   n.risk.K <- object@Ys.K[,ys.name]

                   if (dns.name %in% colnames(object@dNs)) {
                       tp.sum <- data.frame(time=times[ind], estimate=tmp$CI.trans[ind,1,i],variance=tmp$CI.trans[ind,4,i],lower.ci=tmp$CI.trans[ind,2,i],upper.ci=tmp$CI.trans[ind,3,i],n.risk = n.risk[ind], n.event=n.event[ind],n.risk.K = n.risk.K[ind], n.event.K=n.event.K[ind])
                   } else {
                       tp.sum <- data.frame(time=times[ind],estimate=tmp$CI.trans[ind,1,i],variance=tmp$CI.trans[ind,4,i],lower.ci=tmp$CI.trans[ind,2,i],upper.ci=tmp$CI.trans[ind,3,i],n.risk = n.risk[ind],n.remain=n.risk[ind]-n.event[ind],n.risk.K=n.risk.K[ind],n.remain.K=n.risk.K[ind]-n.event.K[ind])
                   }

                   print(tp.sum, row.names = FALSE,digits=digits)
                   cat("\n")
               } #end of for statement
           } #end of else
       } #end of if trans.pr
} #end of summary function
) #ends setMethod

####################################################################
####################################################################

setMethod("plot", signature(x="msSurv", y="missing"),
           function (x, states="ALL", trans="ALL", plot.type="stateocc",
			CI=TRUE,ci.level=0.95,ci.trans="linear",...) {

##Please enter transitions WITHOUT spaces in between, ie:  "11", "12",...

	plot.type=match.arg(plot.type, c("stateocc", "transprob","entry.d","exit.d"))

#	tmp <- MSM.CIs(x,ci.level,ci.trans) #Calling CIs

	#factoring the states,trees
	if(plot.type=="stateocc"){

		tmp <- MSM.CIs(x,ci.level=0.95) #Calling CIs
		if(states[1]=="ALL") states<-nodes(x@tree)

		f.st <- factor(states)
		ls <- length(states)
		sl <- which(nodes(x@tree)%in%as.numeric(states)) #location of states in the matrix

		if(CI==TRUE & !is.null(x@cov.p)){
			rd <- tmp$CI.p
 			dimnames(rd)$dim=gsub("p", "State", dimnames(rd)$dim)
			y <- as.vector(rd[,1,sl])
			y2 <- as.vector(rd[,2,sl]) #lower limit
			y3 <- as.vector(rd[,3,sl]) #upper limit
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)$dim[sl], each=dim(rd)[1]))
			## NOTE: add '...' argument below
			st.plot <- xyplot(y + y2 + y3 ~ x | f.st, allow.multiple=TRUE,type="s",lty=c(1,2,2),col=c(1,2,2),...)
			st.plot <- update(st.plot,main="Plot of State Occupation Probabilites",
				 xlab="Event Times",ylab="State Occupation Probabilities",
				 key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)), 
					 text=list(c("Est", "Lower CI", "Upper CI")),
					 columns=3))
			print(st.plot)
				 
		} #end of CIs TRUE


		if(CI==FALSE){
			rd <- tmp$CI.p
			dimnames(rd)$dim=gsub("p", "State", dimnames(rd)$dim)
			y <- as.vector(rd[,1,sl])
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)$dim[sl], each=dim(rd)[1]))
			st.plot <- xyplot(y~x|f.st, type="s",col=1)
			st.plot <- update(st.plot,main="Plot of State Occupation Probabilites",
				 xlab="Event Times",ylab="State Occupation Probabilities",
				 key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
			print(st.plot)
		} #end of no CIs

	} #end of state occ plot

	if(plot.type=="transprob"){

		tmp <- MSM.CIs(x,ci.level,ci.trans) #Calling CIs
		all.trans <- x@pos.trans
		cc <- strsplit(all.trans," ")
		cc2 <- sapply(cc,function(x) x[1])
		cc3 <- sapply(cc,function(x) x[2])
		all.trans <- paste(cc2,cc3,sep="")	

		if(trans[1] =="ALL") trans <- all.trans

		rd <- tmp$CI.trans
		names(rd) <- paste(trans,"transition")
		tr <- which(all.trans%in%trans) #location of states in the matrix

		if(CI==TRUE & !is.null(x@out)){
				
			y <- as.vector(rd[,1,tr])
			y2 <- as.vector(rd[,2,tr]) #lower limit
			y3 <- as.vector(rd[,3,tr]) #upper limit
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(trans))
			f.tp <- as.factor(rep(dimnames(rd)$dim[tr], each=dim(rd)[1]))
			tr.plot <- xyplot(y + y2 + y3 ~ x | f.tp, allow.multiple=TRUE,type="s",lty=c(1,2,2),col=c(1,2,2),...)
			tr.plot <- update(tr.plot,main="Plot of Transition Probabilites",
				 xlab="Event Times",ylab="Transition Probabilites",
				 key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)), 
					 text=list(c("Est", "Lower CI", "Upper CI")),
					 columns=3))
			print(tr.plot)


		} #end of CIs TRUE


		if(CI==FALSE){
			rd <- tmp$CI.trans
			y <- as.vector(rd[,1,tr])
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(trans))
			f.tp <- as.factor(rep(dimnames(rd)$dim[tr], each=dim(rd)[1]))
			tr.plot <- xyplot(y~x|f.tp, allow.multiple=FALSE,type="s",col=1,...)
			tr.plot <- update(tr.plot,main="Plot of Transition Probabilites",xlab="Event Times",ylab="Transition Probabilities",
				 key = list(lines=list(col=1, lty=1), text=list("Est"),columns=1))
			print(tr.plot)

		} #end of no CIs

	} #end of trans prob plot


##1/10/11 -- adding plots for entry/exit distn
##1/19/11 -- adding CIs for entry/exit distn


	if(plot.type=="entry.d"){

	        enter <- as.character(which(!(sapply(inEdges(x@tree), function(x) length(x) == 0))))
		if(states[1]=="ALL") states<-enter

		f.st <- factor(states)
		ls <- length(states)
		sl <- which(nodes(x@tree)%in%as.numeric(states)) #location of states in the matrix


		if(CI==TRUE){

		   if(!is.null(x@Fs.var)){
			tmp <- Dist.CIs(x,ci.level,ci.trans) #Calling CIs
			rd <- tmp$CI.Fs
 			dimnames(rd)$dim=gsub("F", "State", dimnames(rd)$dim)
			y <- as.vector(rd[,1,sl])
			y2 <- as.vector(rd[,2,sl]) #lower limit
			y3 <- as.vector(rd[,3,sl]) #upper limit
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)$dim[sl], each=dim(rd)[1]))
			ent.plot <- xyplot(y + y2 + y3 ~ x | f.st, allow.multiple=TRUE,type="s",lty=c(1,2,2),col=c(1,2,2))
			ent.plot <- update(ent.plot,main="Plot of State Entry Time Distributions",
				 xlab="Event Times",ylab="State Entry Time Distributions",
				 key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)), 
					 text=list(c("Est", "Lower CI", "Upper CI")),
					 columns=3))

			print(ent.plot)
		   }

		   else {
			#cat("Warning:  Covariance for state entry distributions is NULL and therefore CIs not plotted. \n")
			rd <- x@Fs
			dimnames(rd)[[2]]=gsub("F", "State", dimnames(rd)[[2]])
			y <- as.vector(rd[,sl])
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
			ent.plot <- xyplot(y~x|f.st, type="s",col=1)
			ent.plot <- update(ent.plot,main="Plot of State Entry Time Distributions",
				 xlab="Event Times",ylab="State Entry Time Distributions",
				 key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
			print(ent.plot)

		   }
		} #end of CI False

		
#

		if(CI==FALSE){
			rd <- x@Fs
			dimnames(rd)[[2]]=gsub("F", "State", dimnames(rd)[[2]])
			y <- as.vector(rd[,sl])
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
			ent.plot <- xyplot(y~x|f.st, type="s",col=1)
			ent.plot <- update(ent.plot,main="Plot of State Entry Time Distributions",
				 xlab="Event Times",ylab="State Entry Time Distributions",
				 key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
			print(ent.plot)
		} #end of CI False

	} #end of entry distribution plot


	if(plot.type=="exit.d"){

		transient <- as.character(which(sapply(edges(x@tree), function(x) length(x) > 0)))
		if(states[1]=="ALL") states<-transient

		f.st <- factor(states)
		ls <- length(states)
		sl <- which(nodes(x@tree)%in%as.numeric(states)) #location of states in the matrix


		if(CI==TRUE){

		   if(!is.null(x@Gs.var)){
			tmp <- Dist.CIs(x,ci.level,ci.trans) #Calling CIs
			rd <- tmp$CI.Gs
			dimnames(rd)$dim=gsub("G", "State", dimnames(rd)$dim)
			y <- as.vector(rd[,1,sl])
			y2 <- as.vector(rd[,2,sl]) #lower limit
			y3 <- as.vector(rd[,3,sl]) #upper limit
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)$dim[sl], each=dim(rd)[1]))
			ent.plot <- xyplot(y + y2 + y3 ~ x | f.st, allow.multiple=TRUE,type="s",lty=c(1,2,2),col=c(1,2,2),...)
			ent.plot <- update(ent.plot,main="Plot of State Exit Time Distributions",
				 xlab="Event Times",ylab="State Exit Time Distributions",
				 key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)), 
					 text=list(c("Est", "Lower CI", "Upper CI")),
					 columns=3))

			print(ent.plot)
		   }


		   else{
#			cat("Warning:  Covariance for state exit distributions is NULL and therefore CIs not plotted. \n")
			rd <- x@Gs
			dimnames(rd)[[2]]=gsub("G", "State", dimnames(rd)[[2]])
			y <- as.vector(rd[,sl])
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
			exit.plot <- xyplot(y~x|f.st, type="s",col=1)
			exit.plot <- update(exit.plot,main="Plot of State Exit Time Distributions",
				 xlab="Event Times",ylab="State Exit Time Distributions",
				 key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
			print(exit.plot)
		   } #end of null variance loop
		} #end of CI FALSE loop


		if(CI==FALSE){
			rd <- x@Gs
			dimnames(rd)[[2]]=gsub("G", "State", dimnames(rd)[[2]])
			y <- as.vector(rd[,sl])
			x <- rep(as.numeric(dimnames(rd)[[1]]), length(states))
			f.st <- as.factor(rep(dimnames(rd)[[2]][sl], each=dim(rd)[1]))
			exit.plot <- xyplot(y~x|f.st, type="s",col=1)
			exit.plot <- update(exit.plot,main="Plot of State Exit Time Distributions",
				 xlab="Event Times",ylab="State Exit Time Distributions",
				 key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
			print(exit.plot)
		} #end of CI FALSE loop


	} #end of entry distribution plot


}#end of function

)
