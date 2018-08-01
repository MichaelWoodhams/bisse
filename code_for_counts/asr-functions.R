#Functions to hide ASR details

# This code is from http://blog.revolutionanalytics.com/2014/10/r-in-production-controlling-runtime.html
eval_fork <- function(..., timeout = 60) {
  myfork <- parallel::mcparallel({
    eval(...)
  }, silent = FALSE)
  
  # wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)
  # kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)
  
  # clean up:
  parallel::mccollect(myfork, wait = FALSE)
  # timeout?
  if (is.null(myresult))
    stop("reached elapsed time limit")
  
  # move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]]
  
  # send the buffered response
  return(myresult)
}


#modification to slot in place of asr.mk2
asr.mk2.randomStarts <- function(tree, pars, numStarts = 2) {
  messages <- c()
	#try and recover the ancestral states using maximum likelihood
  lik.m <- make.mk2(tree, tree$tip.state)
  #print("asr.mk2.randomStarts find.mle base call")
	best.fit.m <- find.mle(lik.m, pars, method="subplex") ##supply required initial rate estimates
	best.lnL <- best.fit.m$lnLik
	for (i in 1:(numStarts)) {
	  #print(paste("asr.mk2.randomStarts find.mle jitter call",i))
		fit.m <- find.mle(lik.m, jitter.pars(pars), method="subplex") #idea is to perturb the pars each time
		if (fit.m$lnLik > best.lnL){
		  messages <- append(messages,paste("mk2 accepting jitter call",i," for delta lnL=",fit.m$lnLik-best.lnL,"\n"))
			best.fit.m <- fit.m
			best.lnL <- fit.m$lnLik
		}
	}
	st.m <- asr.marginal(lik.m, coef(best.fit.m))

	return(list(data=st.m,messages=messages))
}

#modification to slot in place of asr.bisse
asr.bisse.randomStarts <- function(tree, pars, numStarts = 2) {
  messages <- c()
	## BiSSE ancestral state reconstructions under the ML model
	lik.s <- make.bisse(tree, tree$tip.state)
	#print("asr.bisse.randomStarts find.mle base call")
	fit.s <- find.mle(lik.s, pars, method="subplex")
	best.fit.s <- find.mle(lik.s, pars, method="subplex") ##supply required initial rate estimates
	best.lnL <- best.fit.s$lnLik
	for (i in 1:(numStarts)) {
	  #print(paste("asr.bisse.randomStarts find.mle jitter call",i))
		fit.s <- find.mle(lik.s, jitter.pars(pars), method="subplex") #idea is to perturb the pars each time
		if (fit.s$lnLik > best.lnL){
		  messages <- append(messages,paste("bisse accepting jitter call",i," for delta lnL=",fit.s$lnLik-best.lnL,"\n"))
			best.fit.s <- fit.s
			best.lnL <- fit.s$lnLik
		}
	}

	st.s <- asr.marginal(lik.s, coef(best.fit.s))

	return(list(data=st.s,messages=messages))
}

jitter.pars <- function(pars) {
	#not quite sure what a senisble thing here is
	#at the moment I am multiplying each entry of pars by a factor between 0.5 and 1.5
	multfac <- runif(length(pars))+0.5

	out <- pars
	for(i in 1:length(pars))
		out[[i]] = pars[[i]] * multfac[i]

	return(out)
}

asr.bisseno_opt <- function(tree, pars) {
	## bisse_opt ancestral state reconstructions under the ML model
	lik.ss <- make.bisse(tree, tree$tip.state)
	fit.ss <- find.mle(lik.ss, pars, method="subplex")
	st.ss <- asr.marginal(lik.ss, pars)

	return(st.ss)
}
asr.mpr <- function(tree, numTaxa) {
	#MPR expects the column of data, an unrooted tree, an arbitrary outgroup
	anc.pars.MPR <- MPR(tree$tip.state, unroot(tree), tree$tip.label[1])

	#get the output of anc.pars in the same format as st.m
	temp.MPR <- matrix(rep(0,2*(numTaxa-1)), nrow=(numTaxa-1))
	rownames(temp.MPR) <- tree$node.label
	root.node <- tree$node.label[1]
	temp.MPR[root.node,] = c(0,0) #tree was unrooted so this node disappeared
	for (node in rownames(anc.pars.MPR)){
		if (anc.pars.MPR[node,"lower"]==0 & anc.pars.MPR[node,"upper"]==0){
			temp.MPR[node,] = c(1,0)
		}
		if (anc.pars.MPR[node,"lower"]==1 & anc.pars.MPR[node,"upper"]==1){
			temp.MPR[node,] = c(0,1)
		}
		if (anc.pars.MPR[node,"lower"]==0 & anc.pars.MPR[node,"upper"]==1) {
			temp.MPR[node,] = c(0.5,0.5)
		}
	}

	mpr <- t(temp.MPR)

	return(mpr)
}

calculateNodeWeights <- function(states, probState) {
	#work out how correct each node is, i.e. the probability weight on right answer
	zeroMask <- 1 - states
	rightness <- probState[1,]*zeroMask + probState[2,]*states

	return(rightness)
}
