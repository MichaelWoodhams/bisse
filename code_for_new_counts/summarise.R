library(data.table)
library(reshape2)

treeFiles = list.files("deep-dark-forest/treeOutput", pattern="-tree.txt", full.names = TRUE)
treeSummary = do.call("rbind", lapply(treeFiles, function(file) {
	message(paste("processing", file))
	read.delim(file)
}))

transitionFiles = list.files("plantation/treeOutput", pattern="-tree.txt", full.names = TRUE)
transitionSummary = do.call("rbind", lapply(transitionFiles, function(file) {
	message(paste("processing", file))
	read.delim(file)
}))

merged = merge(treeSummary, transitionSummary)

summariseNodes = function(nodeFile) {
	message(paste("processing", nodeFile))

  nodes = read.delim(nodeFile)

  nodes$incorrect.MP=ifelse(nodes$rightness.mpr<0.4,1,0)
  nodes$incorrect.Mk=ifelse(nodes$rightness.mk2<0.4,1,0)
  nodes$incorrect.BiSSE=ifelse(nodes$rightness.bisse<0.4,1,0)
  nodes$incorrect.BiSSE_no_opt=ifelse(nodes$rightness.bisseno_opt<0.4,1,0)

  #Rank nodes by age
  nodes$nodeRank = ave(-nodes$nodeAge, nodes$treeID, FUN=rank)

	#Summarise counts of incorrect notes and incorrect nodes deeper than 41
    summary = data.table(nodes)[,list(
  	incorrect.MP = sum(incorrect.MP) - 1, #root node for MP is not reconstructed and always incorrect, so ignore it
  	incorrect.Mk = sum(incorrect.Mk),
  	incorrect.BiSSE = sum(incorrect.BiSSE),
  	incorrect.BiSSE_no_opt = sum(incorrect.BiSSE_no_opt),

  	num.nodes = length(nodeID),

  	incorrect.MP.deep = sum(incorrect.MP[nodeRank < 41]) - 1, #root node for MP is not reconstructed and always incorrect, so ignore it
  	incorrect.Mk.deep = sum(incorrect.Mk[nodeRank < 41]),
  	incorrect.BiSSE.deep = sum(incorrect.BiSSE[nodeRank < 41]),
  	incorrect.BiSSE_no_opt.deep = sum(incorrect.BiSSE_no_opt[nodeRank < 41]),
  	num.nodes.deep = length(nodeID[nodeRank < 41])

  ), by=c("treeID","paramSetID")]

	return(summary)
}

nodeFiles = list.files("deep-dark-forest/nodeOutput", pattern="-node.txt", full.names = TRUE)

nodeSummary = do.call("rbind", lapply(nodeFiles, summariseNodes))

summary = merge(merged, nodeSummary, by=c("treeID", "paramSetID"))

transitionBins = c(0, 5, 10, 15, 20, 30, 50, Inf)
transitionBinLabels = c("1_5", "6_10", "11_15", "16_20", "21_30", "31_50", "50")
summary$transitionGroup = cut(summary$trueTransitions, transitionBins, labels = transitionBinLabels)

deepIncorrectBins = c(-1, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, Inf)
deepIncorrectBinLabels = (0:10) * 10

summary$MP.deep.cat = cut(summary$incorrect.MP.deep, deepIncorrectBins, labels = deepIncorrectBinLabels)
summary$Mk.deep.cat = cut(summary$incorrect.Mk.deep, deepIncorrectBins, labels = deepIncorrectBinLabels)
summary$BiSSE.deep.cat = cut(summary$incorrect.BiSSE.deep, deepIncorrectBins, labels = deepIncorrectBinLabels)
summary$BiSSE_no_opt.deep.cat = cut(summary$incorrect.BiSSE_no_opt.deep, deepIncorrectBins, labels = deepIncorrectBinLabels)

#we actually calculated some extra, unfeasible parameter sets, so throw them away
mu0s = 0.01
mu1s = c(0.01, 0.25, 0.5, 0.8)
qs = c(0.01, 0.05, 0.1)
lambdas = c("0.2 1.8")

out = expand.grid(mu0 = mu0s, mu1 = mu1s, q01 = qs, q10 = qs, lambdas = lambdas)

mus = c(0.01, 0.25, 0.5, 0.8)
qs = c(0.01, 0.05, 0.1)
lambdas = c("0.5 1.5", "1 1", "1.5 0.5")

out = rbind(out, expand.grid(mu0 = mus, mu1 = mus, q01 = qs, q10 = qs, lambdas = lambdas))

mu0s = c(0.01, 0.25, 0.5, 0.8)
mu1s = 0.01
qs = c(0.01, 0.05, 0.1)
lambdas = c("1.8 0.2")

out = rbind(out, expand.grid(mu0 = mu0s, mu1 = mu1s, q01 = qs, q10 = qs, lambdas = lambdas))
out$lambda0 = sapply(out$lambdas, function(lambdas) {strsplit(as.character(lambdas), ' ')[[1]][[1]]})
out$lambda1 = sapply(out$lambdas, function(lambdas) {strsplit(as.character(lambdas), ' ')[[1]][[2]]})
out$lambdas = NULL
out = merge(out, summary, all=FALSE)

#write the summary
write.csv(out, file="tree.summary.csv",row.names = FALSE)

summary = out