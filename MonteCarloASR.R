#setwd("/home/mdw2/Documents/botany/bisse/R/code_for_counts")

# Command line arguments:
# 1) tree summary output file name
# 2) Node summary output file name
# 3) lambda0
# 4) lambda1
# 5) mu0
# 6) mu1
# 7) q01
# 8) q10
# 9) id (just an identifying string naming this scenario)
# 10) (optional) number of trees to simulate.
# 11) (optional) number of tree to start at. 
# 12) (optional) number of tips on tree.

library(diversitree)
library(phangorn)
library(phytools)
#library(R.utils)

Rcpp::sourceCpp('optimised.cpp')
source("asr-functions.R")

# Standard values:
numTaxa = 400 #set how big we want the trees to be
numTrees = 200 #set how many trees we want per parameter setting. (Also overridable on cmd line)
numRandomStarts = 1 # This is in addition to one start from true parameters. 
max.asr.runtime = 1800 # Give up if ASR hasn't completed for a single tree in this time (seconds).
testForceAbortReps = c() # empty, unless you're testing what happens when a calculation times out.

# Override standard values here:
#numTrees = 10
#numRandomStarts = 3
#testForceAbortReps = c(2,4)

args <- commandArgs(trailingOnly = TRUE) # if called from cmd line, will have enough args here
if (length(args) < 9) {
  print("Too few args on cmd line, using hard coded args")
  args <- c("myTreeOutput/1-1-tree.txt", "myNodeOutput/1-1-node.txt", 0.5, 1.5, 0.01, 0.01, 0.01, 0.01, "1-1",10)
  #args <- c("myTreeOutput/1-4-tree.txt", "myNodeOutput/1-4-node.txt", 0.5, 1.5, 0.8, 0.01, 0.01, 0.01, "1-4", 4)
  #args <- c("myTreeOutput/1-4-tree.txt", "myNodeOutput/1-4-node.txt", 0.5, 1.5, 0.8, 0.01, 0.01, 0.01, "1-4", 1, 2)
  #args <- c("myTreeOutput/3-631-tree.rerun.txt", "myNodeOutput/3-631-node.rerun.txt", 1.8, 0.2, 0.5, 0.25, 0.01, 0.05, "3-631")
} 

output.file = args[[1]]
node.file = args[[2]]
outputFileHeadersWritten = FALSE

pars <- list(
	lambda0=as.numeric(args[[3]]),
	lambda1=as.numeric(args[[4]]),
	mu0=as.numeric(args[[5]]),
	mu1=as.numeric(args[[6]]),
	q01=as.numeric(args[[7]]),
	q10=as.numeric(args[[8]]),
	id=args[[9]])

if (length(args)>9) {
  numTrees=as.numeric(args[[10]])
}
if (length(args)>10) {
  firstTree=as.numeric(args[[11]])
} else {
  firstTree=1
}
if (length(args)>11) {
  numTaxa=as.numeric(args[[12]])
}


#Set the random seed using permutatuib identifier.
#all non-numeric characters in the ID are replaced with zeros. e.g., "2-27"
#becomes seed 2027
baseSeed <- as.numeric(gsub("[^0-9]", "0", pars$id))
#set.seed(baseSeed)

#evolve a set of trees under these parameter settings
#only subset of parameters used for bisse
pars.bisse <- c(pars$lambda0, pars$lambda1, pars$mu0, pars$mu1, pars$q01, pars$q10)
#someTrees <- trees(pars.bisse, "bisse", max.taxa=numTaxa, max.t=50, x0=0, n=numTrees, include.extinct=FALSE)

# I don't have any particular use for these counters, but might be interesting.
countExtinct <- 0 # How many trees did we generate where the lineage went extinct?
countInvar0 <- 0 # How many trees did we generate where all tips were in state 0?
countInvar1 <- 0 # How many trees did we generate where all tips were in state 1?
countTreeTooLong <- 0; # How many times did we reach max.t before numTaxa? Expect never.
# This one is different - it reduces the number of tree analyses we have.
countTimedOut <- 0 # How many trees had calculation aborted due to infinite loop?

allData <- c() # vector of records (lists), one record per tree.

decile <- ceiling((1:(numTaxa-1))/numTaxa*9.999) # Deciles of node depth: (nearly) equal numbers of 1's, 2's etc to 10's, in order.

# Because sometimes a calculation is finished by a time out, we need to separately count
# replicates (RNG seeds) and successful calculations.
rep <- firstTree
countCompleted <- 0
while (countCompleted < numTrees) {
  record <- list(id=pars$id,rep=rep)
  # Note: max integer is .Machine$integer.max == .Machine$integer.max.
  # Max basSeed is 30720. Mult by 10,000 avoids overflow of integer.
  set.seed(baseSeed*10000+rep)
  
  maxRepeat <- 10000
  for (i in 1:maxRepeat) {
    #oneTree <- trees(pars.bisse, "bisse", max.taxa=numTaxa, max.t=50, x0=0, n=1, include.extinct=FALSE)[[1]]
    oneTree <- tree.bisse(pars.bisse, max.taxa=numTaxa, max.t=50,x0=0)
    if (is.null(oneTree)) {
      countExtinct <- countExtinct+1
    } else if (oneTree$Nnode != numTaxa-1) {
      countTreeTooLong <- countTreeTooLong+1
    } else if (sum(oneTree$tip.state)==0) {
      countInvar0 <- countInvar0+1
    } else if (sum(oneTree$tip.state)==numTaxa) {
      countInvar1 <- countInvar1+1
    } else {
      break
    }
  }
  # Use 'for' rather than 'repeat' as sanity check
  if (i==maxRepeat) {
    # Should never happen, as maxRepeat is set fairly high.
    stop("Could not get valid tree after ",maxRepeat," tries")
  }
  
  print(paste("ID", pars$id,": Doing tree", rep))
	  
	#Do ASR.
	# I've had just one case where calculation went into endless loop, hence the timeout protection.
	# Changes to RNG seeding means I can't give an example where this happens for current code.

	analysis <- try(eval_fork({
	  # DEBUG: force occasional errors to test success counting
	  if (rep %in% testForceAbortReps) {stop("forced abort for testing")}
    # I don't know why, but parameter jittering is inconsistent (depends on how many
		# previous tree analyses have happened in this run) unless we reseed here.
		set.seed(baseSeed*10100+rep)
		mk2 <- asr.mk2.randomStarts(oneTree, c(pars$q01, pars$q10),numStarts=numRandomStarts);
		bisse <- asr.bisse.randomStarts(oneTree, pars.bisse,numStarts=numRandomStarts);
		bisseno_opt <- asr.bisseno_opt(oneTree, pars.bisse);
		mpr <- asr.mpr(oneTree, numTaxa);
		list(mk2=mk2,bisse=bisse,bisseno_opt=bisseno_opt,mpr=mpr)
	},timeout=max.asr.runtime),silent=TRUE)
	if (class(analysis) == "list") {
	  record$rep              <- rep
	  record$mk2.params       <- analysis$mk2$parameters
	  record$mk2.asr          <- analysis$mk2$data[2,] # Vector of ASR probability of being in state 1, by internal node
    record$bisse.params     <- analysis$bisse$parameters
    record$bisse.asr        <- analysis$bisse$data[2,]
    record$bisseno_opt.asr <- analysis$bisseno_opt[2,]
    record$mp.asr           <- analysis$mpr[2,]
    printWhichBest(record$mk2.params,"MK2",rep)
    printWhichBest(record$bisse.params,"BiSSE",rep)
    countCompleted <- countCompleted+1
	 } else {
	  print(paste("ID", pars$id,"Timed out tree", rep))
	   countTimedOut <- countTimedOut+1
	   rep <- rep+1
	  next
	}		
	allData <- append(allData,record)
	#Make a dataframe to store the tree results
	treeResult = data.frame(paramSetID=pars$id, treeID=rep)
	treeResult = cbind(treeResult,data.frame(pars[1:6])) # pars[7] is 'id' again.

	#Get true history of tree
	history.from.sim.discrete(oneTree, 0:1)
	treeResult$trueTransitions <- length(oneTree$hist$to)
	
	treeResult$numTips <- length(oneTree$tip.label)
	treeResult$treeLength<- sum(oneTree$edge.length) #new
	
	
	states <- oneTree$node.state
	numNodes = length(states)

	rightness.mk2 <- calculateNodeWeights(states, record$mk2.asr)
	rightness.bisse <- calculateNodeWeights(states, record$bisse.asr)
	rightness.bisseno_opt <- calculateNodeWeights(states, record$bisseno_opt.asr)
	rightness.mpr <- calculateNodeWeights(states, record$mp.asr)

	#summary measure of how accurate the ASRs are:
	#sum the probabilities assigned to the correct state across all internal nodes

	treeResult$Mk <- sum(rightness.mk2)/numNodes
	treeResult$BiSSE <- sum(rightness.bisse)/numNodes
	treeResult$BiSSEno_opt <- sum(rightness.bisse)/numNodes
	treeResult$MP <- sum(rightness.mpr)/numNodes

	#minCharChange needs the tip states to be separate from the tree for some
	#stupid reason, so remove them and calculate this last
	tips <- oneTree$tip.state
	oneTree$tip.state <- NULL
	minCharChanges <- minCharChangeFast(oneTree, trait = tips, printMinResult = FALSE)

	treeResult$infChanges <- length(minCharChanges$sumTransitions)

	#Append current tree results to output file (make file if headers not writted)
	write.table(treeResult, output.file, append= outputFileHeadersWritten, row.names=FALSE, col.names=!outputFileHeadersWritten, sep="\t")

	#for node output we need the node ID, node depth, and true state (calculated above in states)
	nodeResults <- data.frame(paramSetID = pars$id,
													  treeID = rep,
													  nodeID = oneTree$node.label,
													  nodeAge = branching.times(oneTree),
													  trueState = states,
														rightness.mk2 = rightness.mk2,
														rightness.bisse = rightness.bisse,
														rightness.bisseno_opt = rightness.bisseno_opt,
														rightness.mpr = rightness.mpr)
	nodeResults <- nodeResults[order(nodeResults$nodeAge),]
	nodeResults$decile <- decile

	#Append current tree node results to output file (make file if first tree)
	write.table(nodeResults, node.file, append=outputFileHeadersWritten, row.names=FALSE, col.names=!outputFileHeadersWritten, sep="\t")

	outputFileHeadersWritten=TRUE
	rep <- rep+1
}
print(paste(countExtinct,      " trees regenerated due to extinction."))
print(paste(countInvar0,       " trees regenerated due to all tips state 0."))
print(paste(countInvar1,       " trees regenerated due to all tips state 1."))
print(paste(countTreeTooLong,  " trees regenerated due to tree depth limit."))
print(paste(countTimedOut,     " trees had calculations aborted for taking too long."))
# Some extra complication to allow for possibility output.file does not end in '.txt'
save(allData,file=paste(gsub("txt.RData","RData",paste(output.file,".RData",sep=""))))
