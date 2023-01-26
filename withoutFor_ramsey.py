import itertools
import sys
import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support

import networkx as nx
from networkx.classes.function import degree
import numpy as np
#from memory_profiler import profile
from networkx.algorithms.cuts import edge_expansion
from numpy.random import *
from scipy.linalg import circulant
from scipy.special import comb
from tqdm import tqdm
import random
import math

n = 13
m = 1
d = int(n / m)

N = int(n*(n-1)/2)

k1 = 3
k2 = 5

# Degree of the vertices
degree = [0] * (n+1)

# Stores the vertices
store = [0]* n

thr = 24

def integer_to_binary(p,n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    return bistr

def calcGraph(i, m):
    G1 = nx.Graph()
    G2 = nx.Graph()
    graphMat = np.zeros((n, n))

    #make binary Vector
    if m == 1:
        bistr = integer_to_binary(i, math.floor(n//2))
    else:
        length = m * math.floor(d//2) + int(comb(m,2))*(math.floor(d//2) + 1)
        bistr = integer_to_binary(i, length)

    #print(bistr)

    #make Circulant Matrix
    cirIdxList = [1]
    cirIdxNum = m + 1
    for j in range(1, m):
        cirIdxList.append(cirIdxNum)
        cirIdxNum = cirIdxNum + (m-j)
    
    getIdxNum = 0
    putIdxNum = 0
    cirMatList = []
    if m != 1:
        choose_Number_of_Color = int((len(bistr) - len(cirIdxList)*math.floor(d/2)) / (m*(m-1)/2))

    for j in range(1, int(comb(m+1, 2))+1):
        cirVec = [0]*d
        if j in cirIdxList:
            b = bistr[getIdxNum:getIdxNum + math.floor(d/2)]
            getIdxNum += len(b)
            for k in range(1, len(b)+1):
                c = b.pop(0)
                cirVec[k] = c
                cirVec[-k] = c
        else:
            b = bistr[getIdxNum:getIdxNum + choose_Number_of_Color]
            getIdxNum += len(b)
            for k in range(len(b)):
                c = b.pop(0)
                cirVec[k] = c
                cirVec[-k] = c

        cirMat = circulant(cirVec)
        cirMatList.append(cirMat)
    
    idx = 0
    for k in range(m):
        for l in range(k, m):
            if k == l:
                graphMat[k*d:(k+1)*d, k*d:(k+1)*d] = cirMatList[idx]
                idx += 1
            else:
                graphMat[k*d:(k+1)*d, l*d:(l+1)*d] = cirMatList[idx]
                graphMat[l*d:(l+1)*d, k*d:(k+1)*d] = cirMatList[idx].T
                idx +=1

    print(graphMat)
    secondColorgraph = np.ones((n,n)) - graphMat

    edgesList1 = matrixToedges(graphMat)
    G1.add_edges_from(edgesList1)

    edgesList2 = matrixToedges(secondColorgraph)
    G2.add_edges_from(edgesList2)

    #print(graphMat)
    #print(secondColorgraph)

    amount1 = find_cliques_size_k(G1, k1)
    #print(amount1)

    amount2 = find_cliques_size_k(G2, k2)
    #print(amount2)
    print("k1: ", amount1)
    print("k2: ", amount2)
    print("amount: " + str(amount1 + amount2))

    if amount1 + amount2 == 0:
        print("0 clique: ",i)
    #v = np.array(A)
    #print(v)
    #return A

def find_cliques_size_k(G, k):
    i = 0
    #for j in nx.find_cliques(G):
        #print("clique" ,j)
    for clique in nx.find_cliques(G):
        if len(clique) == k:
            i += 1
        elif len(clique) > k:
            #i += len(list(itertools.combinations(clique, k)))
            i += int(comb(len(clique), k))
    return i

def matrixToedges(adj_mat):
    #print(adj_mat)
    same_edge = np.where(adj_mat == 1)
    i_component = same_edge[0].tolist()
    j_component = same_edge[1].tolist()
    i2_component = []
    j2_component = []
    edges = []
    
    for k in range(len(i_component)):
        if i_component[k] <= j_component[k]:
            i2_component.append(i_component[k])
            j2_component.append(j_component[k])
    
    for k in range(len(i2_component)):
        edges.append([i2_component[k], j2_component[k]])

    size = len(edges)

    for i in range(size):
        degree[edges[i][0]] += 1
        degree[edges[i][1]] += 1
    return edges

def main():
    print(N)

    start = time.time()
    #for i in tqdm(range(8000)):
    #coloringVector = random.randint(0, 2**n-1) + 1
    #coloringVector = 10
    coloringVector = int('011110', 2)
    print(coloringVector)
            #pro = p.map(calcGraph, range(2**N))
        #for i in tqdm(range(2**N)):
    calcGraph(coloringVector, m)

        #with Pool(thr) as p:
            #pro = p.map(calcGraph, tqdm(range(2**N)))
        
    elapsed_time = time.time() - start
    print(elapsed_time)

if __name__ == "__main__":
    main()
