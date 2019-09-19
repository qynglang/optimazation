import numpy as np
from random import shuffle
#This file show an application of GA algorithm on Rastrigin Function.
POPSIZE = 10
MAX_GEN = 1000
P_CROSS = 0.5
#cross overate
P_MUT = 0.001
#mutation rate
m=np.array([128, 64, 32, 16, 8, 4, 2, 1])
#the matrix used to transfer binary strings into real number.
NUMBER_OF_DIMENSION=8
#the length of string. This is a one demetionalproblem
min=1000
class GENE(object):
	def __init__(self,gene):
		self.gene = gene
		self.f = fitness(self.gene)

def fitness(gene):
        n=np.array(gene)
        if m.shape!=n.shape:
                print('error gene',gene)
        
        gene=randomSolution()
        s=n.dot(m)
        s=s.sum()
        x=-5.12+10.24*s/255
        f=10+pow(x,2)-10*np.cos(2*np.pi*x)
        #caculate fitness fuction, as in CRO.
        return f

def mut(gene):
        #mutation
        for i in range(len(gene)):
                if np.random.random() < P_MUT:
                    gene[i]=1-gene[i]
                return GENE(gene)
def cros(gene1,gene2):
        #crossover
        if np.random.random() < P_CROSS:
            a,b=np.random.randint(0,NUMBER_OF_DIMENSION,2)
            if a<b:
                gene_m=np.hstack((gene2[0:a],gene1[a:b],gene2[b:NUMBER_OF_DIMENSION]))
                gene_n=np.hstack((gene1[0:a],gene2[a:b],gene1[b:NUMBER_OF_DIMENSION]))
            else:
                gene_m=np.hstack((gene2[0:b],gene1[b:a],gene2[a:NUMBER_OF_DIMENSION]))
                gene_n=np.hstack((gene1[0:b],gene2[b:a],gene1[a:NUMBER_OF_DIMENSION]))
        else:
            gene_n=gene1
            gene_m=gene2
        if np.array(gene_n).shape!=m.shape:
                print('error cros',gene_n)
        if np.array(gene_m).shape!=m.shape:
                print('error cros',gene_m)
        M=GENE(gene_n)
        N=GENE(gene_m)

        return M,N
def randomSolution():
        return np.random.randint(0,2,NUMBER_OF_DIMENSION)
def main():
        NUMBER_OF_DIMENSION=8
        min=1000
        m=np.array([128, 64, 32, 16, 8, 4, 2, 1])
        pop = [GENE(randomSolution()) for _ in range(POPSIZE)]
        for q in range (0,POPSIZE):
                if pop[q].f<min:
                        min=np.min(pop[q].f)
        for w in range (0,MAX_GEN):
                if w%50==0:
                        print(min)
                k1,k2=np.random.randint(0,POPSIZE,2)
                pop1,pop2=cros(pop[k1].gene,pop[k2].gene)
                
                pop.append(pop1)
                pop.append(pop2)
                if pop[k1].f<min:
                        min=pop[k1].f
                if pop[k2].f<min:
                        min=pop[k2].f    
                k3=np.random.randint(POPSIZE)
                pop3=mut(pop[k3].gene)
                pop.append(pop3)
                if pop[k3].f<min:
                        min=pop[k3].f
                o=[]
                for q in range (0,len(pop)):
                    r=pop[q].f
                    o.append(r)
                o=np.array(o)
                index=o.argsort()
                pop_a=[]
                for a in range (0,POPSIZE):
                    pop_a.append(pop[index[a]])
                pop=pop_a
                
if __name__ == '__main__':
	main()





            
            
    
        
