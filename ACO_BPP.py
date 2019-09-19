import numpy as np
from random import shuffle
import pdb
#This script discribe a application of ant colony algorithm on bin packing problem.
MAX_GEN = 50
DECAY = 0.85
Q = 2
DIM=9
m=np.array([6,6,3,3,2,2,2,2,2]) 

def fitness(candidate):
        k=0
        bin=1
        m=np.array([6,6,3,3,2,2,2,2,2])
        #packages to pack
        for i in range (0,DIM):
           k+=m[candidate[i]]
           if k>10:
                k=m[candidate[i]]
                bin+=1
                #caculate the number of bins needed.
        return bin

def ant_tours(pheromone):
	ants = []
	for a in range(DIM):
		ant = [a]
		taboo = set([a])

		for _ in range(DIM-1):
			i = ant[-1]
			probs = np.array(pheromone[i], copy=True)
			for j in taboo:
				probs[j] = 0
			for j in range(len(probs)):
				g = 1/m[j]
				probs[j] *= 0 if not g else (1/g) ** 2
			probs /= np.sum(probs)
			
			j = np.random.choice(DIM, p=probs)
			ant.append(j)
			taboo.add(j)

		ants.append(ant)

	return ants

def aco():
	# pheromone = np.full((DIM, DIM), 1.0/DIM)
	pheromone = np.random.random((DIM, DIM))
	pheromone /= np.sum(pheromone, axis=0)

	# Seed population
	# pop = [randomSolution() for _ in range(POP_SIZE)]
	# pop_fits = None
	best_score = float('inf')
	# best_score = None
	generation = 0

	while generation < MAX_GEN:
		# compute ant tours
		pop = ant_tours(pheromone)
		scores = [fitness(ant) for ant in pop]
		if not all(s for s in scores):
			pdb.set_trace()

		best_score = min(scores)
		best_ant = pop[np.argmin(scores)]

		if generation % 5 == 0:
			print(best_score)

		# update pheromone matrix 
		for i in range(DIM):
			for j in range(DIM):
				pheromone[i][j] *= DECAY
				# Update with best ant's trail
				if (i,j) in zip(best_ant, best_ant[1:]):
					pheromone[i][j] += Q / best_score

				# OR update with all ants' trails
				# for ant in zip(pop, scores):
				# 	if (i,j) in zip(ant[0], ant[0][1:]):
				# 		pheromone[i][j] += Q / ant[1]

		generation += 1

	return best_score

def main():
	#import QAP
	#d, f = QAP.from_file('qapdata/nug21.dat')
	#optimum, sol = QAP.sol_from_file('qapsoln/nug21.sln')

	#global FLOW, DIST, DIM
	#FLOW = f
	#DIST = d
	#DIM = len(f)

	#print('Goal = ' + str(optimum))
	result = aco()
	#print('\nFINISHED: ' + str(result))


if __name__ == '__main__':
	main()
