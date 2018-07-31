library(diversitree)
library(phangorn)
library(phytools)

countTransitions <- function(tree, tipStates = NULL, nodeStates = NULL, includeTips = TRUE)
{
	if(is.null(tipStates)) tipStates =  tree$tip.state
	if(is.null(nodeStates)) nodeStates =  tree$node.state
	if(!includeTips) tipStates = NULL

	t = tree$orig
	t = merge(data.frame(
		name2 = c(names(tipStates), names(nodeStates)),
		inferredState = c(tipStates, nodeStates)), t)
	t = t[order(t$idx2),]

	u = tree$orig[order(tree$orig$idx2),]

	t$transitionFromParent = mapply(function(parent, inferredState) {
		parentIndex = which(t$idx == parent)

		#There are hidden nodes in the tree for some reason which are included in
		#the internal structure, but not in nodes or tips, so if we encounter one
		#as a parent, find it's first non-hidden ancestor
		while(length(parentIndex)==0)
		{
			#root node has no ancestor, so return NA
			if(!(parent %in% u$idx)) return(NA)

			parent = u$parent[u$idx == parent]
			parentIndex = which(t$idx == parent)
		}

		return(t$inferredState[parentIndex] != inferredState)

	}, t$parent, t$inferredState)

	sum(t$transitionFromParent, na.rm = TRUE)
}

countMprTransitions <- function(tree, mpr, includeTips = TRUE)
{
	states = data.frame(
		name2 = rownames(mpr),
		lower = mpr[,1],
		upper = mpr[,2])

	if(includeTips)
	{
		tipStates = data.frame(
			name2 = names(tree$tip.state),
			lower = tree$tip.state,
			upper = tree$tip.state)

		states = rbind(tipStates, states)
	}

	t = tree$orig
	t = merge(states, t)
	t = t[order(t$idx2),c(1,4,6,8,2,3)]

	u = tree$orig[order(tree$orig$idx2),]

	t$transitionFromParent = mapply(function(parent, Xl, Xu) {
		parentIndex = which(t$idx == parent)

		#There are hidden nodes in the tree for some reason which are included in
		#the internal structure, but not in nodes or tips, so if we encounter one
		#as a parent, find it's first non-hidden ancestor
		while(length(parentIndex)==0)
		{
			#root node has no ancestor, so return NA
			if(!(parent %in% u$idx)) return(NA)

			parent = u$parent[u$idx == parent]
			parentIndex = which(t$idx == parent)
		}

		Yl = t$lower[parentIndex]
		Yu = t$upper[parentIndex]

		if(Xu==Xl && Yu==Yl && Xu==Yl) return(0)
		if(Xu!=Xl && Yu!=Yl) return(0)
		if(Xu==Xl && Yu==Yl && Xu!=Yl) return(1)
		if(Xu==Xl && Yu!=Yl) return(0.5)
		if(Xu!=Xl && Yu==Yl) return(0.5)

		return(NA)

	}, t$parent, t$lower, t$upper)

	sum(t$transitionFromParent, na.rm = TRUE)
}

numTaxa = 400 #set how big we want the trees to be
numTrees = 200 #set how many trees we want per parameter setting

#Read in output paths and parameters from the command line
args <- commandArgs(trailingOnly = TRUE)

output.file = args[[1]]
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
	tree <- someTrees[[rep]]

	#add a check to ensure that we don't use constant characters

	#If boring result is NA, else calculate
	if (table(tree$tip.state)[1]==0 | table(tree$tip.state)[1]==numTaxa
			| length(tree$tip.label) < numTaxa){

		message(paste("ID", pars$id,"Skipping tree", rep))

	} else {

		message(paste("ID", pars$id,": Doing tree", rep))

		#Make a dataframe to store the tree results
		treeResult = data.frame(pars)
		treeResult$paramSetID = pars$id
		treeResult$id = NULL
		treeResult$treeID = rep

		#Get true history of tree
		treeResult$trueTransitions <- length(tree$hist$to)

		treeResult$numTips <- length(tree$tip.label)
		treeResult$treeLength<- sum(tree$edge.length) #new

		treeResult$countTransitions = countTransitions(tree)
		treeResult$countTransitions_no_tips = countTransitions(tree, nodeStates = tree$node.state, includeTips = FALSE)

		mpr = MPR(tree$tip.state, unroot(tree), tree$tip.label[1])

		treeResult$mpTransitions = countMprTransitions(tree, mpr)
		treeResult$mpTransitions_no_tips = countMprTransitions(tree, mpr, includeTips = FALSE)

		treeResult$numTipsState0 = sum(tree$tip.state == 0)
		treeResult$numTipsState1 = sum(tree$tip.state == 1)

		#Append current tree results to output file (make file if headers not writted)
		write.table(treeResult, output.file, append= outputFileHeadersWritten, row.names=FALSE, col.names=!outputFileHeadersWritten, sep="\t")

		outputFileHeadersWritten=TRUE
	}
}


# pars <- list(
# 	lambda0=0.5,
# 	lambda1=1.5,
# 	mu0=0.01,
# 	mu1=0.25,
# 	q01=0.1,
# 	q10=0.1,
# 	id = "test-1")
#
# output.file = "test-1-trees.txt"

#
# which(!(t$trueTransitions >= t$countTransitions & t$countTransitions >= t$mpTransitions))
# which(!(t$countTransitions >= t$countTransitions_no_tips))
# which(!(t$mpTransitions >= t$mpTransitions_no_tips))
# which(!(t$countTransitions_no_tips >= t$mpTransitions_no_tips))
#
# for(i in which(!(t$countTransitions_no_tips >= t$mpTransitions_no_tips)))
# {
#
# 	tree = someTrees[[i]]
# 	mpr = MPR(tree$tip.state, unroot(tree), tree$tip.label[1])
#
# 	cat(i, ":", which(mpr[,2] != mpr[, 1]), "\n", sep = " ")
# }

