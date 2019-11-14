To replicate our results (under Unix):
Install the needed R libraries (diversitree, phangorn, phytools

$ mkdir treeOutput nodeOutput textOutput
$ makeJobs.pl > run.jobs
$ parallel :::: run.jobs
(This will take a long while - days or weeks)
From R, run
nodeSummary.R. (This will take hours)

makeJobs.pl:
Outputs Unix shell commands suitable to run the simulations over the grid
of model parameters. These commands will write results into subdirectories
treeOutput, nodeOutput, textOutput.

MonteCarloASR.R:
Top level simulation code. Read the code for command line arguments.

asr-functions.R:
Low level simulation code.

optimised.cpp:
C++ code to speed up some operations

nodeSummary.R
Reads simulation results from treeOutput, nodeOUtput and textOutput
and generates R dataframe files treeSummary.{csv|RData},
decileSummary.{csv|RData}, nodeDepthSummary.{csv|RData}
