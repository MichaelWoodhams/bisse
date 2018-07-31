library(diversitree)
library(phangorn)
library(phytools)

Rcpp::sourceCpp('optimised.cpp')
source("asr-functions.R")

numTaxa = 400 #set how big we want the trees to be
numTrees = 200 #set how many trees we want per parameter setting

#Read in output paths and parameters from the command line
args <- commandArgs(trailingOnly = TRUE)

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

#Set the random seed using permutatuib identifier.
#all non-numeric characters in the ID are replaced with zeros. e.g., "2-27"
#becomes seed 2027
set.seed(as.numeric(gsub("[^0-9]", "0", pars$id)))

#evolve a set of trees under these parameter settings
#only subset of parameters used for bisse
pars.bisse <- c(pars$lambda0, pars$lambda1, pars$mu0, pars$mu1, pars$q01, pars$q10)
someTrees <- trees(pars.bisse, "bisse", max.taxa=numTaxa, max.t=50, x0=0, n=numTrees, include.extinct=FALSE)

for (rep in 1:numTrees){
	oneTree <- someTrees[[rep]]

	#add a check to ensure that we don't use constant characters

	#If boring result is NA, else calculate
	if (table(oneTree$tip.state)[1]==0 | table(oneTree$tip.state)[1]==numTaxa
			| length(oneTree$tip.label) < numTaxa){

		print(paste("ID", pars$id,"Skipping tree", rep))

	} else {

		print(paste("ID", pars$id,": Doing tree", rep))

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

		#Do ASR
		mk2 <- asr.mk2.randomStarts(oneTree, c(pars$q01, pars$q10))
		bisse <- asr.bisse.randomStarts(oneTree, pars.bisse)
		bisseno_opt <- asr.bisseno_opt(oneTree, pars.bisse)
		mpr <- asr.mpr(oneTree, numTaxa)

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

