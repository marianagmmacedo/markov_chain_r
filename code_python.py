%matplotlib inline
# config InlineBackend.figure_format = 'retina'
%load_ext autoreload
%autoreload 2

import pandas as pd
import sys
sys.path.insert(0,'../../')
import numpy as np
import matplotlib.pyplot as plt
from swarm_analyzer import SwarmAnalyzer
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import euclidean_distances
import networkx as nx
import scipy
from swarm_plotter import SwarmPlotter
import seaborn as sns 
from random import randint

number_of_runs = 30
topologies = ["regular30"]
function = 27
input_subdir = "path"
filenames = []
particles = 50
dimensions = 1000

dist_matrix_one = []
dist_matrix = []
hist_non = []
simulations = []
import scipy
for topology in topologies:
    for run in range(number_of_runs):
        print(run)
        filename_hdf = "%s/%s_F%02d_%02d_fluctuations_50000_e_100_d_50_p.hdf" % (input_subdir, topology, function, run)
        df_folder = "run_%02d" % run
        df = pd.read_hdf(filename_hdf)
        df = df.transpose()
        matrix_df = df.as_matrix()
        simulations.append(matrix_df)

import math


			#
			#	Markov Chain based on
			# 	good coordinations
			#
			#
			sequence_states = [[] for run in range(number_of_runs)]
			number_particles_connected = []
			for run in range(number_of_runs):
			    for each in range(len(simulations[run])): 
			        sum_histogram = np.sum(simulations[run][each])
			        greater_0_8_ = len(simulations[run][each][ np.where( abs(simulations[run][each]) > 0.8 ) ])
			        greater_0_5_ = len(simulations[run][each][ np.where( abs(simulations[run][each]) > 0.5 ) ]) - greater_0_8_
			        smaller_0_1_ = len(simulations[run][each][ np.where( abs(simulations[run][each]) < 0.1 ) ])
			        smaller_0_3_ = len(simulations[run][each][ np.where( abs(simulations[run][each]) < 0.3 ) ]) - smaller_0_1[run]            
			        if greater_0_8_ > 150:
			            sequence_states[run].append(0) 
			        elif greater_0_8_ > 100: 
			            sequence_states[run].append(1) 
			        elif greater_0_5_ > 100: 
			            sequence_states[run].append(2) 
			        elif smaller_0_1_ > 1000: 
			            sequence_states[run].append(3) 



						#
						#	Markov Chain based on
						# 	sum of coordinations
						#
						#
						sequence_states = [[] for run in range(number_of_runs)]
						sum_histogram_all = []
						for run in range(number_of_runs):
						    for each in range(len(simulations[run])): 
						        sum_histogram = np.sum(simulations[run][each])
						        sum_histogram_all.append(sum_histogram)
						        if sum_histogram > 1000:
						            sequence_states[run].append(0) 
						        elif sum_histogram > 500: 
						            sequence_states[run].append(1) 
						        elif sum_histogram > 126: 
						            sequence_states[run].append(2) 
						        elif sum_histogram > 50: 
						            sequence_states[run].append(3) 
						        else:
						            sequence_states[run].append(4) 



def transition_matrix(transitions):
    n = 1+ max(transitions) #number of states
    M = [[0]*n for _ in range(n)]

    for (i,j) in zip(transitions,transitions[1:]):
        M[i][j] += 1

    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M



matrix_transition = []
number_states = 4
number_runs = 30
m = np.zeros((number_states,number_states))
for run in range(number_runs):
    matrix_transition.append(transition_matrix(sequence_states[run]))
    m = np.add(m, matrix_transition[run])


each_colunm = [[[] for s in range(number_states)] for s in range(number_states)] 
each_state = [[] for s in range(number_runs)] 
for s in range(number_states):
    for run in range(number_runs):
        each_state[s].append(matrix_transition[run][s])
        for l in range(len(each_state[s][run])):
            each_colunm[s][l].append(each_state[s][run][l])

i = 0
for ss in range(number_states):
    for s in range(number_states):
        print('each_state[', end='')
        print(i+1, end='')
        print(',] <- c(', end='')
        print(','.join('{0:.2f}'.format(x) for x in each_colunm[ss][s]), end=',')
        print("\b",end="")
        print(')', end='')
        print("\n")
        i += 1
    

m = np.divide(m,number_runs)
m = m.tolist()
print('numberStates <- %d'%(len(m)))
print(m)


 
output = "%s/regular30_small_apriori_4_states.csv" % (input_subdir)
f1 = open(output, "w")
f1.write("library(igraph)")
f1.write("\nlibrary(expm)")
f1.write("\nrequire(graphics)")
f1.write("\nrequire(markovchain)")
f1.write("\nbyRow <- TRUE")
f1.write('\nnumberStates <- %d'%(len(m)))
f1.write("\nm <- c(")
f1.write(" ,".join(str(x).replace("[", "").replace("]", "") for x in m))
f1.write(")")
f1.write("\nmatrixT <- matrix(data = m, byrow = byRow, nrow = numberStates)")
f1.write("\nrowSums(matrixT)")
f1.write("\nmarkovRun <- new('markovchain', byrow = byRow,transitionMatrix = matrixT, name = 'Markov Chain')")
f1.write("\nmyMc<-as(matrixT, 'markovchain')")
f1.write("\nplot(myMc, color='black', edge.arrow.size=0.3,state.name.size=0.5)")
f1.write("\nshow(markovRun)")
f1.write("\nsummary(markovRun)")
f1.write("\nabsorbingStates(markovRun)")
f1.write("\ntransientStates(markovRun)")
f1.write("\nsteadyStates(markovRun)")
f1.close()
