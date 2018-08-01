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

library(diversitree)
library(phangorn)
library(phytools)
#library(R.utils)

Rcpp::sourceCpp('optimised.cpp')
source("asr-functions.R")

# Standard values, as used in randomStartsAmazon.R:
numTaxa = 400 #set how big we want the trees to be
numTrees = 200 #set how many trees we want per parameter setting. 
numRandomStarts = 1 # This is in addition to one start from true parameters. 
  # randomStartsAmazon uses value 2 but then OBOEs it down to 1.
  # I've fixed the OBOE so value 1 here replicates the randomStartsAmazon behaviour.
# Added by me:
max.asr.runtime = 600 # Give up if ASR hasn't completed for a single tree in this time.

# Override standard values here:
numTrees = 10
numRandomStarts = 1

args <- commandArgs(trailingOnly = TRUE) # if called from cmd line, will have enough args here
if (length(args) < 9) {
  print("Too few args on cmd line, using hard coded args")
  #args <- c("myTreeOutput/1-1-tree.txt", "myNodeOutput/1-1-node.txt", 0.5, 1.5, 0.01, 0.01, 0.01, 0.01, "1-1")
  args <- c("myTreeOutput/3-631-tree.rerun.txt", "myNodeOutput/3-631-node.rerun.txt", 1.8, 0.2, 0.5, 0.25, 0.01, 0.05, "3-631")
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

#Set the random seed using permutatuib identifier.
#all non-numeric characters in the ID are replaced with zeros. e.g., "2-27"
#becomes seed 2027
baseSeed <- as.numeric(gsub("[^0-9]", "0", pars$id))
set.seed(baseSeed)

#evolve a set of trees under these parameter settings
#only subset of parameters used for bisse
pars.bisse <- c(pars$lambda0, pars$lambda1, pars$mu0, pars$mu1, pars$q01, pars$q10)
someTrees <- trees(pars.bisse, "bisse", max.taxa=numTaxa, max.t=50, x0=0, n=numTrees, include.extinct=FALSE)

for (rep in 1:numTrees){
#for (rep in 1:7){
	oneTree <- someTrees[[rep]]
	#add a check to ensure that we don't use constant characters

	#If boring result is NA, else calculate
	if (table(oneTree$tip.state)[1]==0 | table(oneTree$tip.state)[1]==numTaxa
			| length(oneTree$tip.label) < numTaxa){

		print(paste("ID", pars$id,"Skipping tree", rep))

	} else {

		print(paste("ID", pars$id,": Doing tree", rep))
	  # Note: max integer is .Machine$integer.max == .Machine$integer.max.
	  # Max basSeed is 30720. Mult by 10,000 avoids overflow of integer.
	  set.seed(baseSeed*10000+rep)
	  

		#Make a dataframe to store the tree results
		treeResult = data.frame(pars)
		treeResult$paramSetID = pars$id
		treeResult$id = NULL
		treeResult$treeID = rep

		#Get true history of tree
		history.from.sim.discrete(oneTree, 0:1)
		treeResult$trueTransitions <- length(oneTree$hist$to)

		treeResult$numTips <- length(oneTree$tip.label)
		treeResult$treeLength<- sum(oneTree$edge.length) #new

		#Do ASR.
		# I've had just one case where calculation went into endless loop, hence the timeout protection. (scenario 3-631, tree #8 out of 10.)
#		timedout <- withTimeout({
#		  mk2 <- asr.mk2.randomStarts(oneTree, c(pars$q01, pars$q10),numStarts=numRandomStarts);
#		  bisse <- asr.bisse.randomStarts(oneTree, pars.bisse,numStarts=numRandomStarts);
#		  bisseno_opt <- asr.bisseno_opt(oneTree, pars.bisse);
#		  mpr <- asr.mpr(oneTree, numTaxa);
#		},timeout=max.asr.runtime,onTimeout="silent")
#		if (is.null(timedout)) {
#		  print(paste("ID", pars$id,"Timed out tree", rep))
#		  break
#		}		
#		mk2 <- 0; bisse <- 0; bisseno_opt <- 0; mpr <- 0;
		timedout <- try(eval_fork({
		  mk2 <- asr.mk2.randomStarts(oneTree, c(pars$q01, pars$q10),numStarts=numRandomStarts);
		  bisse <- asr.bisse.randomStarts(oneTree, pars.bisse,numStarts=numRandomStarts);
		  bisseno_opt <- asr.bisseno_opt(oneTree, pars.bisse);
		  mpr <- asr.mpr(oneTree, numTaxa);
		  list(mk2=mk2,bisse=bisse,bisseno_opt=bisseno_opt,mpr=mpr,seed=.Random.seed)
		},timeout=max.asr.runtime),silent=TRUE)
		if (class(timedout) == "list") {
      mk2         <- timedout$mk2$data
      bisse       <- timedout$bisse$data
      bisseno_opt <- timedout$bisseno_opt
      mpr         <- timedout$mpr
      #.Random.seed <- timedout$seed
      cat(paste(  timedout$mk2$messages,collapse="\n"))
      cat(paste(timedout$bisse$messages,collapse="\n"))
	  } else {
		  print(paste("ID", pars$id,"Timed out tree", rep))
		  next
		}		
		
		states <- oneTree$node.state
		numNodes = length(states)

		rightness.mk2 <- calculateNodeWeights(states, mk2)
		rightness.bisse <- calculateNodeWeights(states, bisse)
		rightness.bisseno_opt <- calculateNodeWeights(states, bisseno_opt)
		rightness.mpr <- calculateNodeWeights(states, mpr)

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

		#Append current tree node results to output file (make file if first tree)
		write.table(nodeResults, node.file, append=outputFileHeadersWritten, row.names=FALSE, col.names=!outputFileHeadersWritten, sep="\t")

		outputFileHeadersWritten=TRUE
	}
}

