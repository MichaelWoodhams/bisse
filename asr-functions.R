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
  #try and recover the ancestral states using maximum likelihood
  lik.m <- make.mk2(tree, tree$tip.state) # lik.m is a function
  best.fit.m <- find.mle(lik.m, pars, method="subplex") ##supply required initial rate estimates
  best.lnL <- best.fit.m$lnLik
  fitParams <- array(c(best.lnL,best.fit.m$par))
  for (i in 1:(numStarts)) {
    #jittered <- jitter.pars(pars)
    #messages <- append(messages,paste("Jittered pars: ",jittered[[1]],jittered[[2]],"\n"))
    #fit.m <- find.mle(lik.m, jittered, method="subplex") #idea is to perturb the pars each time
    # Above 3 lines debugging. Delete and replace with:
    fit.m <- find.mle(lik.m, jitter.pars(pars), method="subplex") #idea is to perturb the pars each time
    fitParams <- rbind(fitParams,c(lnL=fit.m$lnLik,fit.m$par))
    if (fit.m$lnLik > best.lnL){
      best.fit.m <- fit.m
      best.lnL <- fit.m$lnLik
    }
  }
  st.m <- asr.marginal(lik.m, coef(best.fit.m))

  return(list(data=st.m,parameters=fitParams))
}

#modification to slot in place of asr.bisse
asr.bisse.randomStarts <- function(tree, pars, numStarts = 2) {
  ## BiSSE ancestral state reconstructions under the ML model
  lik.s <- make.bisse(tree, tree$tip.state)
  fit.s <- find.mle(lik.s, pars, method="subplex")
  best.fit.s <- find.mle(lik.s, pars, method="subplex") ##supply required initial rate estimates
  best.lnL <- best.fit.s$lnLik
  fitParams <- array(c(best.lnL,best.fit.s$par))
  for (i in 1:(numStarts)) {
    #jittered <- jitter.pars(pars)
    #messages <- append(messages,paste("Jittered pars: ",jittered[[1]],jittered[[2]],jittered[[3]],jittered[[4]],jittered[[5]],jittered[[6]],"\n"))
    #fit.s <- find.mle(lik.s, jittered, method="subplex") #idea is to perturb the pars each time
    # Above 3 lines debugging. Delete and replace with:
    fit.s <- find.mle(lik.s, jitter.pars(pars), method="subplex") #idea is to perturb the pars each time
    fitParams <- rbind(fitParams,c(lnL=fit.s$lnLik,fit.s$par))
    if (fit.s$lnLik > best.lnL){
      best.fit.s <- fit.s
      best.lnL <- fit.s$lnLik
    }
  }
  st.s <- asr.marginal(lik.s, coef(best.fit.s))

  return(list(data=st.s,parameters=fitParams))
}

jitter.pars <- function(pars) {
	#not quite sure what a senisble thing here is
	#at the moment I am multiplying each entry of pars by a factor between 0.5 and 1.5
	multfac <- runif(length(pars))+0.5

	out <- pars
	for(i in 1:length(pars))
		out[[i]] = pars[[i]] * multfac[i]
  # Debugging, expect to delete this line later
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
	# The unrooting of the tree makes root node disappear, so we need to
	# reconstruct it.
	# It might be that root is always "nd1", first two children always "nd2", "nd3"
	# If I were confident this was true, we could simplify greatly:
	#temp.MPR["nd1",]<-(temp.MPR["nd2",]+temp.MPR["nd3",])*0.5
	
	# This code is less compact than possible, for clarity.
	root.node <- getRoot(tree) # A number, always numTaxa+1?
	root.label <- tree$node.label[root.node-numTaxa]
	first.children <- Children(tree,root.node) # two numbers, either numTaxa+(node number) or, rarely, tip.number
	# node states can sum to 0 (from {0,0}), 0.5 ({0,0.5} or {0.5,0}), 1 ({0,1},{0.5,0.5},{1,0}), 1.5 ({0.5,1},{1,0.5}), or 2 ({1,1}).
	# We want root state to be set to 0, 0, 0.5, 1, 1 respectively in these cases. This 'round' call does that.
	# c(tree$tip.state,tree$node.state) allows for the possibility that one of the children is a tip
	root.state <- round(sum(c(tree$tip.state,tree$node.state)[first.children]))/2
	temp.MPR[root.label,] = c(1-root.state,root.state) 
	
	mpr <- t(temp.MPR)
	return(mpr)
}

printWhichBest <- function(parameters,method,rep) {
  likelihoods <- parameters[,"lnL"]
  minL <- min(likelihoods)
  maxL <- max(likelihoods)
  maxIndex <- match(maxL,likelihoods)
  print(sprintf("Repl %d method %s max lnL from attempt %d, lnL range %f",rep,method,maxIndex,maxL-minL))
}

calculateNodeWeights <- function(states, probState1) {
	#work out how correct each node is, i.e. the probability weight on right answer
	zeroMask <- 1 - states
	rightness <- probState1*states + (1-probState1)*(1-states)
	return(rightness)
}
