# This is an application of Chemical Reaction Optimization on bin packing problem.

import random
import numpy as np
import pdb

POPSIZE = 9
KE_LOSSRATE = 0.8
MOLEC_COLLRATE = 0.2
KE_INITIAL = 1000000

ALPHA = 1300
BETA = 10000

ITERATIONS = 50

NUMBER_OF_DIMENSION = 9

FLOW = None
DIST = None
m=np.array([6,6,3,3,2,2,2,2,2])      
m=np.array(m)
class M(object):
	def __init__(self, o, k, min_hit=0, num_hit=0):
		self.omega = o
		self.PE = fitness(o)
		self.KE = k
		self.min_hit = min_hit
		self.num_hit = num_hit
		


def fitness(omega):
        k=0
        bin=1
        m=np.array([6,6,3,3,2,2,2,2,2]) 
        #package to pack
        for i in range (0,NUMBER_OF_DIMENSION):
           k+=m[omega[i]]
           if k>10:
                k=m[omega[i]]
                bin+=1
                #add another bin once it the item beyond bin size
        return bin

def randomSolution():
        return np.random.permutation(NUMBER_OF_DIMENSION)

def neighbor(mol):
	a1,a2 = np.random.choice(NUMBER_OF_DIMENSION,2,replace=False)
	omega_p = list(mol.omega)
	omega_p[a1] = mol.omega[a2]
	omega_p[a2] = mol.omega[a1]

	M1=M(omega_p, mol.KE, num_hit=mol.num_hit+1)

	if (M1.PE < mol.PE):
		M1.min_hit = M1.num_hit
	else:
		M1.min_hit = mol.min_hit
	
	return M1

def decompose(mol):
	a1,a2 = np.random.randint(-NUMBER_OF_DIMENSION, NUMBER_OF_DIMENSION,2)
	M1=M(np.roll(mol.omega,a1).tolist(), mol.KE, num_hit=mol.num_hit+1)
	M2=M(np.roll(mol.omega,a2).tolist(), mol.KE, num_hit=mol.num_hit+1)

	if (M2.PE < mol.PE):
		M2.min_hit = M2.num_hit
	else:
		M2.min_hit = mol.min_hit

	if (M2.PE < mol.PE):
		M2.min_hit = M2.num_hit
	else:
		M2.min_hit = mol.min_hit

	return M1, M2

def collide(M1, M2):
	return neighbor(M1), neighbor(M2)

def combine(mol1, mol2):
	dif_list = []
	o_prime = np.full(NUMBER_OF_DIMENSION, -1)

	for i in range(NUMBER_OF_DIMENSION):
		if mol1.omega[i] == mol2.omega[i]:
			o_prime[i] = mol1.omega[i]
			dif_list.append(-1)
		else:
			dif_list.append(mol1.omega[i])

	while any(a == b or a == c for (a,b,c) in zip(dif_list, mol1.omega, mol2.omega)):
		random.shuffle(dif_list)

	for i in range(NUMBER_OF_DIMENSION):
		if o_prime[i] == -1:
			# while dif_list[0] == mol1.omega[i] or dif_list[0] == mol2.omega[i]:
			# 	item = dif_list.pop(0)
			# 	dif_list.append(item)
			while dif_list[0] == -1:
				dif_list.pop(0)
			o_prime[i] = dif_list.pop(0)

	M_prime = M(o_prime, mol1.KE, num_hit=mol1.num_hit+mol2.num_hit)

	if (M_prime.PE < mol1.PE and M_prime.PE < mol2.PE):
		M_prime.min_hit = M_prime.num_hit
	else:
		M_prime.min_hit = (mol1.min_hit + mol2.min_hit)

	return M_prime

def ineff_coll_on_wall(M, CEB):
	# obtain neighbor solution omega_p of M
	# compute objective of omega_p PE
	# if new solution PE <= old PE + KE:
	# 	compute new energies
	#	add energy to the CEB
	#	update profile of M

	M_prime=neighbor(M)
	if M.PE+M.KE >= M_prime.PE:
		p=random.uniform(KE_LOSSRATE,1)
		new_energy = M.PE+M.KE-M_prime.PE
		M_prime.KE = new_energy*p
		CEB = CEB+new_energy*(1-p)
		M = M_prime
	return M, CEB

def decomposition(M, CEB):
	if (M.num_hit - M.min_hit) < ALPHA:
		return False, None, None, CEB

	# obtain neighbor solutions omega_p1 and omega_p2
	# calculate PE for candidates
	# calculate leftover energy from candidates
	# if able to decompose without CEB:
	#	compute new energies for candidates
	#	create new molecules M1 and M2
	# else if CEB needed to decompose and is capable:
	#	compute new energies for candidates using CEB energy
	#	update CEB
	#	create new molecules M1 and M2
	# else:
	#	fail

	M_prime1, M_prime2 = decompose(M)
	new_energy = M.PE+M.KE-M_prime1.PE-M_prime2.PE
	success = False

	if new_energy>0:
		success=True
		k=random.random()
		M_prime1.KE=new_energy*k
		M_prime2.KE=new_energy*(1-k)
	elif new_energy+CEB >= 0:
		success=True
		rs=[random.random() for _ in range(4)]
		M_prime1.KE=(new_energy+CEB)*rs[0]*rs[1]
		M_prime2.KE=(new_energy+CEB-M_prime1.KE)*rs[2]*rs[3]
		CEB=(new_energy+CEB)-M_prime1.KE-M_prime2.KE		
	return success, M_prime1, M_prime2, CEB

def inter_ineff_col(M1, M2):
	# obtain neighbor solutions omega_p1 and omega_p2
	# calculate PE for candidates
	# calculate leftover energy from candidates
	# if able to collide:
	#	compute new energies for candidates
	#	update profiles of molecules M1 and M2

	M_prime1, M_prime2 = collide(M1, M2)
	new_energy = (M1.PE+M1.KE+M2.PE+M2.KE)-(M_prime1.PE+M_prime2.PE)
	if new_energy >= 0:
		p=random.uniform(0,1)
		M_prime1.KE=new_energy*p
		M_prime2.KE=new_energy*(1-p)
		M1 = M_prime1
		M2 = M_prime2
	return M1, M2

def synthesis(M1, M2):
	if M1.KE > BETA or M2.KE > BETA:
		return False, None

	# obtain new solution omega_p from omega_1 and omega_2
	# calculate PE for candidate
	# calculate leftover energy from candidate
	# if able to synthesize:
	#	compute new energy for candidate
	#	create new molecule M
	# else:
	#	fail

	M_prime = combine(M1, M2)
	success = False
	new_energy = (M1.PE+M1.KE+M2.PE+M2.KE)
	if new_energy >= M_prime.PE:
		success = True
		M_prime.KE = new_energy - M_prime.PE
		M = M_prime
	return success, M

def CRO():
	# initialize variables
	pop = [M(randomSolution(), KE_INITIAL) for _ in range(POPSIZE)]
	CEB = 0
	min_M = None
	min_F = float('inf')

	# while reacting...
	reactions = 0
	reaction_score = {'uni_inef':0, 'uni_dec':0, 'bi_syn':0, 'bi_inef':0}
	reaction_count = {'uni_inef':0, 'uni_dec':0, 'bi_syn':0, 'bi_inef':0}
	reason = ''
	while reactions <= ITERATIONS:
		if reactions % 5 == 0:
                        print(min_F)
			#print('ITER: ' + str(reactions) + '\tBEST: ' + str(min_F) + '\tPOP: ' +str(len(pop)))
			# print('  oldest molecule: ' + str(sorted(pop, key=lambda x: x.num_hit, reverse=True)[0].num_hit))
			# print('improvements: ' + str(reaction_score))
			# print('reactions: ' + str(reaction_count))
	#	if intermolecular collision:
		if random.random() < MOLEC_COLLRATE and len(pop) > 1:
			m1, m2 = random.sample(pop, 2)
	#		try synthesizing
			syn_success, m_prime = synthesis(m1, m2)
			if syn_success:
				reaction_count['bi_syn'] += 1
				pop.remove(m1)
				pop.remove(m2)
				pop.append(m_prime)
				reason = 'bi_syn'
	#		otherwise inter_ineff_col
			else:
				reaction_count['bi_inef'] += 1
				m_prime1, m_prime2 = inter_ineff_col(m1, m2)
				pop.remove(m1)
				pop.remove(m2)
				pop.append(m_prime1)
				pop.append(m_prime2)
				reason = 'bi_inef'

		else:
			m1 = random.sample(pop, 1)[0]
	#		try decomposition
			dec_success, m_prime1, m_prime2, CEB = decomposition(m1, CEB)
			if dec_success:
				reaction_count['uni_dec'] += 1
				pop.remove(m1)
				pop.append(m_prime1)
				pop.append(m_prime2)
				reason = 'uni_dec'
	#		otherwise ineff_col_on_wall
			else:
				reaction_count['uni_inef'] += 1
				m_prime, CEB = ineff_coll_on_wall(m1, CEB)
				pop.remove(m1)
				pop.append(m_prime)
				reason = 'uni_inef'

	#	check for new min point
		current_min = min(pop, key=lambda m: m.PE)
		if current_min.PE < min_F:
			min_F = current_min.PE
			min_M = current_min
			reaction_score[reason] += 1

		reactions += 1

	#	if stopping criteria met:
	#		stop
	# print('improvements: ' + str(reaction_score))
	# print('reactions: ' + str(reaction_count))
	return min_M.omega, min_F


def main():
	# accept arguments
	# begin CRO
	#import QAP
	#d, f = QAP.from_file('qapdata/nug22.dat')
	#optimum, sol = QAP.sol_from_file('qapsoln/nug22.sln')

	#global FLOW, DIST, NUMBER_OF_DIMENSION
	#FLOW = f
	#DIST = d
	#NUMBER_OF_DIMENSION = len(f)
        #optimum=0
	#print('Goal = ' + str(optimum))
        result = CRO()
	# print(result)




if __name__ == '__main__':
	main()
