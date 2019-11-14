#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool allColsUnique(CharacterMatrix solMat)
{
	/*
	 * all(!sapply(1:ncol(solMat), function(x) any(sapply((1:ncol(solMat))[-x], function(y) identical(solMat[, x], solMat[, y])))))
	 */

	const int nCol = solMat.ncol();
	const int nRow = solMat.nrow();

	//Check for duplicate columns
	for(int colA=0;colA<nCol-1;colA++)
	{
		for(int colB=colA+1;colB<nCol;colB++)
		{
			int match = true;
			for(int row=0;row<nRow;row++)
			{
				if(solMat(row,colA) != solMat(row,colB))
				{
					match = false;
					break;
				}
			}

			if(match)
				return false;
		}
	}

	//All columns are unique
	return true;
}

// [[Rcpp::export]]
SEXP getDuplicateColumns(CharacterMatrix solMat)
{
	/*
	* c(FALSE, sapply(2:ncol(solMat), function(x) any(sapply((1:ncol(solMat))[1:(x - 1)], function(y) identical(solMat[, x], solMat[,
	 y])))))
	*/

	const int nCol = solMat.ncol();
	const int nRow = solMat.nrow();

	std::vector<bool> isDuplicate (nCol, false);

	for(int x=1;x<nCol;x++)
		for(int y=0;y<x-1;y++)
		{
			int match = true;
			for(int row=0;row<nRow;row++)
			{
				if(solMat(row,x) != solMat(row,y))
				{
					match = false;
					break;
				}
			}

			if(match)
				isDuplicate[x] = true;
		}

	return wrap(isDuplicate);
}

/*** R
library(paleotree)

minCharChangeFast <- function (trait, tree, randomMax = 10000, maxParsimony = TRUE,
					orderedChar = FALSE, type = "MPR", cost = NULL, printMinResult = TRUE,
					ambiguity = c(NA, "?"), dropAmbiguity = FALSE, polySymbol = "&",
					contrast = NULL)
{
	ancMat <- ancPropStateMat(trait, tree, orderedChar = orderedChar, type = type, cost = cost)
	taxSol <- apply(ancMat, 1, function(x) sum(x > 0))
	nSol <- prod(taxSol)
	charN <- colnames(ancMat)
	if (nSol > randomMax) {
		solMat <- t(apply(ancMat, 1, function(x) sample(charN[x > 0], randomMax, replace = T)))
	}	else {
		noChange <- which(taxSol == 1)
		solMat <- matrix(sapply(noChange, function(x) charN[ancMat[x, ] > 0]), , 1)
		rownames(solMat) <- noChange
		if (nSol > 1) {
			for (i in 2:max(taxSol)) {
				changers <- which(taxSol == i)
				for (j in changers) {
					solMat2 <- lapply(charN[ancMat[j, ] > 0], function(x) rbind(solMat, x))
					solMat1 <- solMat2[[1]]
					for (k in 2:length(solMat2)) {
						solMat1 <- cbind(solMat1, solMat2[[k]])
					}
					colnames(solMat1) <- NULL
					rownames(solMat1) <- c(rownames(solMat), j)
						solMat <- solMat1
				}
			}
		}
		solMat <- solMat[order(as.numeric(rownames(solMat))),
                   , drop = FALSE]
	}
	solUnq <- allColsUnique(solMat)
		if (!solUnq) {
			if (nSol > randomMax) {
				solDup <- getDuplicateColumns(solMat)
				solMat <- solMat[, !solDup]
			}
			else {
				stop("Not all solutions are unique, as calculated, despite random permutations not used. Please investigate or contact Dave Bapst.")
			}
		}
		edgeSol <- array(, dim = c(Nedge(tree), 2, ncol(solMat)))
			for (i in 1:ncol(solMat)) {
				xSol <- solMat[, i]
				edgeSol[, , i] <- cbind(xSol[sapply(tree$edge[, 1], function(x) which(x ==
					names(xSol)))], xSol[sapply(tree$edge[, 2], function(x) which(x ==
					names(xSol)))])
			}
			tranMat <- array(, dim = c(length(charN), length(charN),
                              ncol(solMat)))
				rownames(tranMat) <- paste("anc.", colnames(ancMat), sep = "")
				colnames(tranMat) <- paste("desc.", colnames(ancMat), sep = "")
				sumTran <- numeric()
				for (i in 1:ncol(solMat)) {
					edgeTran <- edgeSol[, , i]
					tranMat1 <- t(sapply(charN, function(x) sapply(charN,
                                          function(y) sum(edgeTran[, 1] == x & edgeTran[, 2] ==
                                                   y))))
					tranMat[, , i] <- tranMat1
					diag(tranMat1) <- 0
					sumTran[i] <- sum(tranMat1)
				}
				if (nSol == 1) {
					maxPars <- 1
					tranMat <- tranMat[, , 1, drop = FALSE]
					edgeSol <- edgeSol[, , 1, drop = FALSE]
				}
				else {
					maxPars <- which(sumTran == min(sumTran))
					if (maxParsimony) {
						solMat <- solMat[, maxPars, drop = FALSE]
						tranMat <- tranMat[, , maxPars, drop = FALSE]
						sumTran <- sumTran[maxPars]
					}
				}
				tranSumChange <- t(sapply(lapply(1:dim(tranMat)[3], function(y) tranMat[,
                                     , y]), function(x) c(sum(x[upper.tri(x)]), sum(diag(x)),
                                     sum(x[lower.tri(x)]))))
					colnames(tranSumChange) <- c("Gains", "NoChanges", "Losses")
					minTran <- apply(tranMat, c(1, 2), min)
					funcMess <- c(paste0(nSol, " potential solutions under ",
                          type, ", ", length(maxPars), " most parsimonious solutions found"),
                          ifelse(nSol > randomMax, "Solutions sampled stochastically",
                                 "Solutions exhaustively checked"))
					if (printMinResult) {
						if (length(maxPars) < 6) {
							print(list(message = funcMess, sumTransitions = sumTran,
                  transitionArray = tranMat, minTransitions = minTran))
						}
						else {
							print(list(message = funcMess, sumTransitions = sumTran,
                  minTransitions = minTran))
						}
					}
					return(invisible(list(message = funcMess, sumTransitions = sumTran,
                           minTransitions = minTran, solutionArray = edgeSol, transitionArray = tranMat,
                           transitionSumChanges = tranSumChange)))
}

*/
