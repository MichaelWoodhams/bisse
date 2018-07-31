mus = c(0.01, 0.25, 0.5, 0.8)
qs = c(0.01, 0.05, 0.1)
lambdas = c("0.5 1.5", "1 1", "1.5 0.5")

out = expand.grid(mu0 = mus, mu1 = mus, q01 = qs, q10 = qs, lambdas = lambdas, idPrefix = 1)

mu0s = 0.01
mu1s = c(0.01, 0.25, 0.5, 0.8)
qs = c(0.01, 0.05, 0.1)
lambdas = c("0.2 1.8")

out = rbind(out, expand.grid(mu0 = mus, mu1 = mus, q01 = qs, q10 = qs, lambdas = lambdas, idPrefix = 2))

mu0s = c(0.01, 0.25, 0.5, 0.8)
mu1s = 0.01
qs = c(0.01, 0.05, 0.1)
lambdas = c("1.8 0.2")

out = rbind(out, expand.grid(mu0 = mus, mu1 = mus, q01 = qs, q10 = qs, lambdas = lambdas, idPrefix = 3))

scriptPath = "calculate.R"
r = "/usr/local/bin/R"

cat(file = "jobs.jobs")

for(i in 1:nrow(out))
{
	id = paste(out$idPrefix[i], i, sep="-")
	treeOutput = paste("treeOutput/",id,"-tree.txt", sep = "")

	cat("[ -e", treeOutput," ] ||",r,"<", scriptPath, "--slave --args", treeOutput, "none", as.character(out$lambdas[i]), out$mu0[i], out$mu1[i], out$q01[i], out$q10[i], id, "\n",
			file = "jobs.jobs",
			append = TRUE)
}