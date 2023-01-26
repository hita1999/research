import itertools
import math
import sys
import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support

import networkx as nx
import numpy as np
from memory_profiler import profile
from networkx.algorithms.cuts import edge_expansion
from networkx.classes.function import degree
from numpy.random import *
from scipy.linalg import circulant
from scipy.special import comb
from tqdm import tqdm
import random

#the number of colors
r = 2
#the size of tabu list
L = 250
K = 50
#the order of graph
n = 14
#the number of edges
N = int(n*(n-1)/2)
#the number of division number
m = 7
#block circulant matrix size
d = int(n / m)
#subgraph of first color
k1 = 3
#subgraph of second color
k2 = 5
# Degree of the vertices
degree = [0] * (n+1)
# Stores the vertices
store = [0]* n
#the number of threads for multiprocessing
thr = 24

def integer_to_binary(p,n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    return bistr

def generateColoringVector(i, m):
    #make binary Vector
    if m == 1:
        bistr = integer_to_binary(i, math.floor(n//2))
    else:
        length = m * math.floor(d//2) + int(comb(m,2))*(math.floor(d//2) + 1)
        bistr = integer_to_binary(i, length)
    return bistr

def setGraph(bistr):
    graphMat = np.zeros((n, n))
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
    return graphMat

def score(settedGraph, w1, w2):
    G1 = nx.Graph()
    G2 = nx.Graph()

    secondColorgraph = np.ones((n,n)) - settedGraph

    edgesList1 = matrixToedges(settedGraph)
    G1.add_edges_from(edgesList1)

    edgesList2 = matrixToedges(secondColorgraph)
    G2.add_edges_from(edgesList2)

    #print(settedGraph)
    #print(secondColorgraph)

    #the number of monochromatic complete graph of first color
    f1 = find_cliques_size_k(G1, k1)
    #print(f1)

    #the number of monochromatic complete graph of second color
    f2 = find_cliques_size_k(G2, k2)

    scoreV = w1 * f1 + w2 * f2

    #print("k1: ", f1)
    #print("k2: ", f2)
    return scoreV, w1, w2, f1, f2
    #print(f2)

    #print("amount: " + str(f1 + f2))

    #if f1 + f2 == 0:
    #print("0 clique: ",i)
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

def my_index_multi(l, x):
    return [i for i, _x in enumerate(l) if _x == x]

def adjustWeight(w1,w2, f1, f2):
    w1 = (K*w1 + f1) / (K+1)*(f1+f2)
    w2 = (K*w2 + f2) / (K+1)*(f1+f2)
    #print(w1, w2)
    return w1, w2

def main():
    print(N)
    #coloring vector
    """
    if m == 1:
        vector_number = random.randint(0, 2**math.floor(n//2)-1)
    else:
        length = m * math.floor(d//2) + int(comb(m,2))*(math.floor(d//2) + 1)
        vector_number = random.randint(0, 2**length)
    """
    #入力!
    vector_number = 5
    print(vector_number)
    start = time.time()
    #with Pool(thr) as p:
            #pro = p.map(calcGraph, tqdm(range(2**N)))
    #initial weight
    w1 = 1
    w2 = 1
    coloringVector = vector_number
    #Set elements of V to random colors
    rmVector = generateColoringVector(coloringVector, m)
    print("length of vector:",len(rmVector))
    H = set()
    M = set()
    scoreV = 1.0
    
    while scoreV > 0.0 or len(M) > 0:
        W = []
        scoreList = []
        for i in range(len(rmVector)/2, len(rmVector)/4 + len(rmVector)/2):
            for c in range(r):
                Vic = []
                Vic = rmVector
                tmp = rmVector[i]
                Vic[i] = c
                if tuple(Vic) not in H:
                    settedGraph = setGraph(Vic)
                    W.append(Vic.copy())
                    scoreV, weight1, weight2, f1, f2 = score(settedGraph, w1, w2)
                    scoreList.append(scoreV)

                    #print("score",scoreV)
                Vic[i] = tmp
            if f1 + f2 == 0:
                break
        if len(scoreList) > 0:  
            minScoreList = my_index_multi(scoreList, min(scoreList))
        if len(W) > 0:
            for j in minScoreList:
                M.add(tuple(W[j]))
        else:
            break
            
        if len(M) != 0:
            rmVector = M.pop()
            H.add(rmVector)
            rmVector = list(rmVector)

            if len(H) > L:
                H.pop()
            
            w1, w2 = adjustWeight(weight1, weight2, f1, f2)
    
    print(rmVector)
    elapsed_time = time.time() - start
    print(elapsed_time)

if __name__ == "__main__":
    main()
