setwd("/home/mdw2/Documents/botany/bisse/results/taxa100")

library(diversitree) # needed

# Standard values:
numTaxa = 400 #set how big we want the trees to be

# Code to take paramSetID or filename starting with a paramSetID and reverse-engineer the parameters:
lookup.lambda0 = c(0.5,1,1.5,0.2,1.8) # indexed by (num-1)/144
lookup.q10     = c(0.01,0.05,0.1) # indexed by (num-1)/48%3
lookup.q01     = c(0.01,0.05,0.1) # indexed by (num-1)/16%3
lookup.mu1     = c(0.01,0.25,0.5,0.8) # indexed by (num-1)/4%4
lookup.mu0     = c(0.01,0.25,0.5,0.8) # indexed by (num-1)%4
numToParams <- function(num) {
  num <- num-1
  params = list()
  params$lambda0 <- lookup.lambda0[num%/%144+1]
  params$lambda1 <- 2-params$lambda0
  params$mu0 <- lookup.mu0[num%%4+1]
  params$mu1 <- lookup.mu1[num%/%4%%4+1]
  params$q01 <- lookup.q01[num%/%16%%3+1]
  params$q10 <- lookup.q10[num%/%48%%3+1]
  return(params)
}
# Given paramSetID or filename starting with paramSetID, extract
# the parameter set numbers and turn them back to parameters.
nameToParams <- function(name) {
  #print(name)
  name <- tail(strsplit(name,"/")[[1]],1) # strips of directory path, if any
  name <- substring(name,3) # remove 1-, 2- or 3- from start of name.
  position <- regexpr(pattern='-',name)[1] # Find second '-' (if it exists)
  if (position>0) {
    name <- substring(name,1,position-1)
  }
  return(numToParams(as.numeric(name)))
}

# I neglected to store tip states. Happily it is fairly straight forward to reconstruct the true
# tree for a given run.

trueTree <- function(id,rep) {
  pars <- nameToParams(id) # "id" is e.g. "1-1"
  pars.bisse <- as.numeric(pars)
  ## Commented code was for when id was passed as an integer.
  #idPrefix <- 1
  #if (pars$lambda0<pars$mu0) {idPrefix <- 2}
  #if (pars$lambda1<pars$mu1) {idPrefix <- 3}
  #id <- sprintf("%d-%d",idPrefix,paramSetID)
  baseSeed <- as.numeric(gsub("[^0-9]", "0", id))
  set.seed(baseSeed*10000+rep)
  maxRepeat <- 1000
  for (i in 1:maxRepeat) {
    #oneTree <- trees(pars.bisse, "bisse", max.taxa=numTaxa, max.t=50, x0=0, n=1, include.extinct=FALSE)[[1]]
    oneTree <- tree.bisse(pars.bisse, max.taxa=numTaxa, max.t=50,x0=0)
    if (is.null(oneTree)) {
      #countExtinct <- countExtinct+1
    } else if (oneTree$Nnode != numTaxa-1) {
      #countTreeTooLong <- countTreeTooLong+1
    } else if (sum(oneTree$tip.state)==0) {
      #countInvar0 <- countInvar0+1
    } else if (sum(oneTree$tip.state)==numTaxa) {
      #countInvar1 <- countInvar1+1
    } else {
      break
    }
  }
  # Use 'for' rather than 'repeat' as sanity check
  if (i==maxRepeat) {
    # Should never happen, as maxRepeat is set fairly high.
    stop("Could not get valid tree after ",maxRepeat," tries")
  }
  return(oneTree)
}

countState1Tips <- function(paramSetID,rep) {
  tree <- trueTree(paramSetID,rep)
  return(sum(tree$tip.state))
}

# Effectively we are rounding the predicted probability:
# prediction in [0,0.3) -> 0
#               [0.3,0.7] -> 0.5
#               (0.7,1] -> 1
# Then compare this to trueState, which is 0 or 1. This gives
# all rightnesses equal to 0, 0.5 or 1.
# So continuous rightnesses are mapped:
# [0,0.3) -> 0
# [0.3,0.7] -> 0.5
# (0.7,1] -> 1
quantizedRightness <- function(rightness,threshold=0.3) {
  quantRight <- rep(0.5,length(rightness))
  quantRight[rightness <   threshold] <- 0
  quantRight[rightness > 1-threshold] <- 1
  quantRight[is.na(rightness)] <- NA
  return(quantRight)
}

# Returns a vector: (num scored 0 and true 0; num scored 1 and true 1; num scored 0.5; num scored 0 and true 1; num scored 1 and true 0)
quantRightCounts <- function(rightness,trueState,threshold=0.3) {
  quantRight <- quantizedRightness(rightness,threshold)
  counts <- c(0,0,0,0,0)
  counts[1] <- sum(quantRight==1 & trueState==0)
  counts[2] <- sum(quantRight==1 & trueState==1)
  counts[3] <- sum(quantRight==0.5)
  counts[4] <- sum(quantRight==0 & trueState==0)
  counts[5] <- sum(quantRight==0 & trueState==1)
  return(counts)
}

# Add a column 'deepest' which is true for the deepest 4 nodes for each tree.
findDeepest <- function(nodes) {
  nodes$deepest <- FALSE
  for (i in unique(nodes$treeID)) {
    nodes$deepest[nodes$treeID==i & order(nodes[nodes$treeID==1,"nodeAge"])+4>=numTaxa] <- TRUE
  }
  return(nodes)
}

# There are no HiSSE results for 100 or 1600 taxon trees, but leave code in for completeness.

# Reads nodes file, adds 'deepest' column and rightness.hisse column,
# which will either come from a Hisse nodes file, or will be full of NAs.
# Also extends 'decile' factor to have a value 100 (used for 'deepest' nodes)
addHisseNodes <- function(file) {
  nodes <- findDeepest(read.table(paste("nodeOutput/",file,sep=""),header=TRUE))
  nodes$decile <- factor(nodes$decile,levels=c(1:10,100)) # add 100 as a factor level
  hisseFile <- paste("hisse/nodeOutput/",file,sep="") # nonexistant dir, will find nothing
  if (file.exists(hisseFile)) {
    hisseNodes <- read.table(hisseFile,header=TRUE)
    nodes <- merge(nodes,hisseNodes[,c("paramSetID","treeID","nodeID","rightness.hisse")],by=c("paramSetID","treeID","nodeID"),all.x=TRUE)
  } else {
    nodes$rightness.hisse = as.numeric(NA)
  }
  return(nodes)
}

# add hisse column to trees dataframes, and add column numTipState1 counting tips in state 1
# Add wrongness counts for quantized predictions
addHisseTrees <- function(file,count.wrong) {
  dataset <- read.table(paste("treeOutput/",file,sep=""),header=TRUE)
  dataset$numTipState1 <- -1
  for (i in 1:length(dataset$treeID)) {
    # 'treeID' may have missing values due to timeouts, so can't assume treeID[i]==i.
    dataset$numTipState1[i] <- countState1Tips(as.character(dataset$paramSetID[i]),dataset$treeID[i])
  }
  hisseFile <- paste("hisse/treeOutput/",file,sep="") # nonexistant dir, will find nothing
  if (file.exists(hisseFile)) {
    hisseTree <- read.table(hisseFile,header=TRUE)
    dataset <- merge(dataset,hisseTree[,c("paramSetID","treeID","HiSSE")],by=c("paramSetID","treeID"),all.x=TRUE)
  } else {
    dataset$HiSSE = as.numeric(NA)
  }
  dataset <- merge(dataset,count.wrong,by=c("paramSetID","treeID"))
}

# Really nasty hack: summarizeFileDeciles puts info into count.wrong which is later used
# in tree summary table.
count.wrong = data.frame()

summarizeFileDeciles <- function(nodes,file) {
  nodes$decile <- factor(nodes$decile,levels=c(1:10,100)) # add 100 as a factor level
  dec <- nodes$decile # for brevity
  params <- nameToParams(file)
  sum <- data.frame(paramSetID           =nodes$paramSetID[[1]],
                    lambda0              =params$lambda0,
                    lambda1              =params$lambda1,
                    mu0                  =params$mu0,
                    mu1                  =params$mu1,
                    q01                  =params$q01,
                    q10                  =params$q10,
                    decile               =levels(dec),
                    nodeAge              =tapply(nodes$nodeAge,              dec,mean),
                    rightness.mk2        =tapply(nodes$rightness.mk2,        dec,mean),
                    rightness.bisse      =tapply(nodes$rightness.bisse,      dec,mean),
                    rightness.bisseno_opt=tapply(nodes$rightness.bisseno_opt,dec,mean),
                    rightness.mpr        =tapply(nodes$rightness.mpr,        dec,mean),
                    rightness.hisse      =tapply(nodes$rightness.hisse,      dec,function(x) {mean(x,na.rm=TRUE)}),
                    q.rightness.mk2        =tapply(quantizedRightness(nodes$rightness.mk2),        dec,mean),
                    q.rightness.bisse      =tapply(quantizedRightness(nodes$rightness.bisse),      dec,mean),
                    q.rightness.bisseno_opt=tapply(quantizedRightness(nodes$rightness.bisseno_opt),dec,mean),
                    q.rightness.hisse      =tapply(quantizedRightness(nodes$rightness.hisse),      dec,function(x) {mean(x,na.rm=TRUE)})
  )
  # And turn currently empty 'decile==100' row into row for the deepest 1% (i.e. four) nodes
  deep <- nodes[nodes$deepest,]
  sum["100","nodeAge"]                <- mean(deep$nodeAge)
  sum["100","rightness.mk2"]          <- mean(deep$rightness.mk2)
  sum["100","rightness.bisse"]        <- mean(deep$rightness.bisse)
  sum["100","rightness.bisseno_opt"]  <- mean(deep$rightness.bisseno_opt)
  sum["100","rightness.mpr"]          <- mean(deep$rightness.mpr)
  sum["100","rightness.hisse"]        <- mean(deep$rightness.hisse,na.rm=TRUE)
  sum["100","q.rightness.mk2"]        <- mean(quantizedRightness(deep$rightness.mk2))
  sum["100","q.rightness.bisse"]      <- mean(quantizedRightness(deep$rightness.bisse))
  sum["100","q.rightness.bisseno_opt"]<- mean(quantizedRightness(deep$rightness.bisseno_opt))
  sum["100","q.rightness.hisse"]      <- mean(quantizedRightness(deep$rightness.hisse),na.rm=TRUE)
  row.names(sum) <- NULL
  return(sum)
}

summarize.wrong.counts <- function(nodes) {
  # Get counts of wrong predictions. Needs to be done here because
  # this is where we have true states, but is saved for output in tree summary file.
  wrong <- data.frame()
  for (rep in unique(nodes$treeID)) {
    thisTree <- nodes[nodes$treeID==rep,]
    q.r.count.mk2         <- quantRightCounts(thisTree$rightness.mk2,        thisTree$trueState)
    q.r.count.bisse       <- quantRightCounts(thisTree$rightness.bisse,      thisTree$trueState)
    q.r.count.bisseno_opt <- quantRightCounts(thisTree$rightness.bisseno_opt,thisTree$trueState)
    q.r.count.mpr         <- quantRightCounts(thisTree$rightness.mpr,        thisTree$trueState)
    q.r.count.hisse       <- quantRightCounts(thisTree$rightness.hisse,      thisTree$trueState)
    q.count.wrong.mk2        =q.r.count.mk2[4]        +q.r.count.mk2[5]
    q.count.wrong.bisse      =q.r.count.bisse[4]      +q.r.count.bisse[5]
    q.count.wrong.bisseno_opt=q.r.count.bisseno_opt[4]+q.r.count.bisseno_opt[5]
    q.count.wrong.mpr        =q.r.count.mpr[4]        +q.r.count.mpr[5]
    q.count.wrong.hisse      =q.r.count.hisse[4]      +q.r.count.hisse[5]
    q.count.ambig.mk2        =q.r.count.mk2[3]        
    q.count.ambig.bisse      =q.r.count.bisse[3]      
    q.count.ambig.bisseno_opt=q.r.count.bisseno_opt[3]
    q.count.ambig.mpr        =q.r.count.mpr[3]        
    q.count.ambig.hisse      =q.r.count.hisse[3]      
    
    
    # Not working: Want totals by tree, this is giving totals by param set.
    wrong.row = data.frame(paramSetID=nodes$paramSetID[[1]],
                           treeID=rep,
                           q.count.wrong.mk2        =q.count.wrong.mk2,
                           q.count.wrong.bisse      =q.count.wrong.bisse,
                           q.count.wrong.bisseno_opt=q.count.wrong.bisseno_opt,
                           q.count.wrong.mpr        =q.count.wrong.mpr,
                           q.count.wrong.hisse      =q.count.wrong.hisse,
                           q.count.ambig.mk2        =q.count.ambig.mk2,
                           q.count.ambig.bisse      =q.count.ambig.bisse,
                           q.count.ambig.bisseno_opt=q.count.ambig.bisseno_opt,
                           q.count.ambig.mpr        =q.count.ambig.mpr,
                           q.count.ambig.hisse      =q.count.ambig.hisse
    )
    wrong = rbind(wrong,wrong.row)
  }
  #names(count.wrong)=c("paramSetID","treeID","q.count.wrong.mk2","q.count.wrong.bisse","q.count.wrong.bisseno_opt","q.count.wrong.mpr")
  return(wrong)
}


# And yet another summary:
# Similar to decileSummary, but with a bit more info:
# rightness for node depths in ranges.
# Counts of nodes in those same ranges.

rangeNames <- c("[0.0,0.1)","[0.1,0.2)","[0.2,0.3)","[0.3,0.5)","[0.5,0.7)","[0.7,1.0)","[1.0,2.0)","[2.0,inf)")
rangeBounds <- c(0,0.1,0.2,0.3,0.5,0.7,1.0,2.0,9999)


summarizeFileRanges <- function(nodes,file) {
  params <- nameToParams(file)
  ageBin <- 0*nodes$treeID
  for (i in 1:length(rangeNames)) {
    ageBin[nodes$nodeAge>=rangeBounds[i] & nodes$nodeAge<rangeBounds[i+1]] = rangeNames[i]
  } 

  ageBin <- as.factor(ageBin)
  sum <- data.frame(paramSetID           =nodes$paramSetID[[1]],
                    lambda0              =params$lambda0,
                    lambda1              =params$lambda1,
                    mu0                  =params$mu0,
                    mu1                  =params$mu1,
                    q01                  =params$q01,
                    q10                  =params$q10,
                    ageBin               =levels(ageBin),
                    rightness.mk2        =tapply(nodes$rightness.mk2,        ageBin,mean),
                    rightness.bisse      =tapply(nodes$rightness.bisse,      ageBin,mean),
                    rightness.bisseno_opt=tapply(nodes$rightness.bisseno_opt,ageBin,mean),
                    rightness.mpr        =tapply(nodes$rightness.mpr,        ageBin,mean),
                    rightness.hisse      =tapply(nodes$rightness.hisse,      ageBin,function(x) {mean(x,na.rm=TRUE)}),
                    q.rightness.mk2        =tapply(quantizedRightness(nodes$rightness.mk2),        ageBin,mean),
                    q.rightness.bisse      =tapply(quantizedRightness(nodes$rightness.bisse),      ageBin,mean),
                    q.rightness.bisseno_opt=tapply(quantizedRightness(nodes$rightness.bisseno_opt),ageBin,mean),
                    q.rightness.hisse      =tapply(quantizedRightness(nodes$rightness.hisse),      ageBin,function(x) {mean(x,na.rm=TRUE)}),
                    binCounts            =as.numeric(table(ageBin))
  )
  row.names(sum) <- NULL
  return(sum)
}


# count.wrong is calculated here, but used in making the trees summary.
deciles <- data.frame()
count.wrong = data.frame()
ranges <- data.frame()
for (file in list.files("nodeOutput")) {
  print(file)
  nodes <- addHisseNodes(file)
  deciles <- rbind(deciles,summarizeFileDeciles(nodes,file))
  ranges <- rbind(ranges,summarizeFileRanges(nodes,file))
  count.wrong <- rbind(count.wrong,summarize.wrong.counts(nodes))
}
save(deciles,file="decileSummary.RData")
write.csv(deciles,file="decileSummary.csv",row.names=FALSE)
save(count.wrong,file="count.wrong.RData") # not used by any external program, just saved for manual inspection and reloading for debugging.
save(ranges,file="nodeDepthSummary.RData")
write.csv(ranges,file="nodeDepthSummary.csv",row.names=FALSE)
print("**** Finished deciles and ranges summaries")

trees <- data.frame()
for (file in list.files("treeOutput",pattern=".*txt")) {
  print(file)
  dataset <- addHisseTrees(file,count.wrong)
  trees <- rbind(trees,dataset)
}
save(trees,file="treeSummary.RData")
write.csv(trees,file="treeSummary.csv",row.names=FALSE)
print("**** Finished trees summary")
