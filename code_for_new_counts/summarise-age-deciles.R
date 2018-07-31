library(data.table)

summariseNodes = function(nodeFile) {
	message(paste("processing", nodeFile))

	nodes = read.delim(nodeFile)

	nodes$incorrect.MP=ifelse(nodes$rightness.mpr<0.4,1,0)
	nodes$incorrect.Mk=ifelse(nodes$rightness.mk2<0.4,1,0)
	nodes$incorrect.BiSSE=ifelse(nodes$rightness.bisse<0.4,1,0)
	nodes$incorrect.BiSSE_no_opt=ifelse(nodes$rightness.bisseno_opt<0.4,1,0)

	nodes$ageDecile = cut(-nodes$nodeAge, breaks = quantile(-nodes$nodeAge, probs = seq(0, 1, 0.1)),
												 include.lowest = TRUE, labels = 1:10)

	summary = data.table(nodes)[,list(
		incorrect.MP = sum(incorrect.MP),
		incorrect.Mk = sum(incorrect.Mk),
		incorrect.BiSSE = sum(incorrect.BiSSE),
		incorrect.BiSSE_no_opt = sum(incorrect.BiSSE_no_opt),

		num.nodes = length(nodeID)
	), by=c("treeID","paramSetID", "ageDecile")]

	#Remove root node from MP, since it's not reconstructed
	summary$incorrect.MP[summary$ageDecile == 1] = summary$incorrect.MP[summary$ageDecile == 1] - 1

	return(summary)
}

nodeFiles = list.files("deep-dark-forest/nodeOutput", pattern="-node.txt", full.names = TRUE)

nodeSummary = do.call("rbind", lapply(nodeFiles, summariseNodes))

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

treeSummary = merge(treeSummary, transitionSummary)
summary = merge(treeSummary, nodeSummary, by=c("treeID", "paramSetID"))

write.csv(summary, file="tree.summary.age.deciles.csv",row.names = FALSE)
