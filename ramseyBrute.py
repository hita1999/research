import itertools
import sys
import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support

import networkx as nx
import numpy as np
from memory_profiler import profile
from networkx.algorithms.cuts import edge_expansion
from numpy.random import *
from scipy.linalg import circulant
from scipy.special import comb
from tqdm import tqdm

n = 6
N = int(n*(n-1)/2)
#MAX = 5

#integerVector
coloringVector = 7

k1 = 3
k2 = 3

# Degree of the vertices
d = [0] * (n+1)

# Stores the vertices
store = [0]* n

thr = 12


def integer_to_binary(p,n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    return bistr

def calcGraph(i):
    G1 = nx.Graph()
    G2 = nx.Graph()
    graphMat = np.zeros((n, n))
    count = 0
   
    bistr = integer_to_binary(i,N)
    for j in range(0,n):
        for k in range(j+1,n):
            b = bistr.pop(0)
            graphMat[j][k] = b
            graphMat[k][j] = b
    
    secondColorgraph = np.ones((n,n)) - graphMat

    edgesList1 = matrixToedges(graphMat)
    G1.add_edges_from(edgesList1)

    edgesList2 = matrixToedges(secondColorgraph)
    G2.add_edges_from(edgesList2)

    amount1 = find_cliques_size_k(G1, k1)

    amount2 = find_cliques_size_k(G2, k2)
    #print("amount: " + str(i+1) + "番目 " + str(amount1 + amount2))

    if amount1+amount2 == 0:
        count += 1
    return amount1+amount2

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
        d[edges[i][0]] += 1
        d[edges[i][1]] += 1
    return edges

def main():
    print(N)
    print("the number of all graphs:",2**N)
    count = 0
    
    start = time.time()

    with Pool(thr) as p:
        res = p.map(calcGraph, tqdm(range(2**N)))
    
    print(res)

    elapsed_time = time.time() - start
    print("the number of 0 clique: ", sum(res))
    print(elapsed_time)

if __name__ == "__main__":
    main()