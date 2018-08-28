install.packages("devtools")
devtools::install_github('spedygiorgio/markovchain')
require(graphics)
library(igraph)
library(expm)
require(markovchain)
byRow <- TRUE
numberStates <- 6
m <- c(0.34,0.00,0.00,0.33,0.33,0.00,0.00,0.50,0.00,0.25,0.25,0.00,1.00,0.00,0.00,0.00,0.00,
       0.00,0.34,0.00,0.00,0.33,0.33,0.00,0.00,0.20,0.00,0.00,0.20,0.60,0.00,0.25,0.25,0.00,
       0.25,0.25)
matrixT <- matrix(data = m, byrow = byRow, nrow = numberStates)
rowSums(matrixT)
markovRun <- new("markovchain", byrow = byRow,transitionMatrix = matrixT, name = "Markov Chain")
myMc<-as(matrixT, "markovchain")
plot(myMc, color="black", edge.arrow.size=0.3,state.name.size=0.5)
show(markovRun)
summary(markovRun)
absorbingStates(markovRun) # it has or not?
transientStates(markovRun) # finite number
steadyStates(markovRun) # stationaries states

