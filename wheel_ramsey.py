import itertools
import sys
import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support

import networkx as nx
from networkx.classes.function import degree
import numpy as np
from networkx.algorithms.cuts import edge_expansion
from numpy.random import *
from scipy.linalg import circulant
from scipy.special import comb
from tqdm import tqdm
import random
import math

def integer_to_binary(p,n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    return bistr

def calcGraph(i, m, n, degree):
    d = int(n/m)

    G1 = nx.Graph()
    G2 = nx.Graph()
    graphMat = np.zeros((n, n), dtype=int)

    #make binary Vector
    if m == 1:
        bistr = integer_to_binary(i, math.floor(n//2))
    else:
        length = m * math.floor(d//2) + int(comb(m,2))*(math.floor(d//2) + 1)
        bistr = integer_to_binary(i, length)


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

    #print(graphMat)
    secondColorgraph = np.ones((n,n), dtype=int) - graphMat

    edgesList1 = matrixToedges(graphMat, degree)
    G1.add_edges_from(edgesList1)

    #print(secondColorgraph)
    edgesList2 = matrixToedges(secondColorgraph, degree)
    G2.add_edges_from(edgesList2)

    return edgesList1, edgesList2, graphMat, secondColorgraph, G1, G2

def matrixToedges(adj_mat, degree):
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

def findCycles(path, graph, cycle):
    start_node = path[0]
    next_node= None
    sub = []

    #visit each edge and each node of each edge
    for edge in graph:
        node1, node2 = edge
        if start_node in edge:
                if node1 == start_node:
                    next_node = node2
                else:
                    next_node = node1
                if not visited(next_node, path):
                        # neighbor node not on path yet
                        sub = [next_node]
                        sub.extend(path)
                        # explore extended path
                        findCycles(sub, graph, cycle);
                elif len(path) > 2  and next_node == path[-1]:
                        # cycle found
                        p = rotate_to_smallest(path);
                        inv = invert(p)
                        if isNew(p, cycle) and isNew(inv, cycle):
                            cycle.append(p)

def invert(path):
    return rotate_to_smallest(path[::-1])

#  rotate cycle path such that it begins with the smallest node
def rotate_to_smallest(path):
    n = path.index(min(path))
    return path[n:]+path[:n]

def isNew(path, cycle):
    return not path in cycle

def visited(node, path):
    return node in path

def find_cliques_size_k(G, k):
    i = 0
    for clique in nx.find_cliques(G):
        if len(clique) == k:
            i += 1
        elif len(clique) > k:
            i += int(comb(len(clique), k))
    return i

def searchWheel(int_pathList, graph, n, w):
    wheel_count = 0
    vertex_all = [int(i) for i in range(n)]
    for i in int_pathList:
        vertex_all.remove(i)
    #print("vertex_candidate:",vertex_all)

    for i in vertex_all:
        count = 0
        #print(graph[i])
        for j in int_pathList:
            #print(graph[i][j])
            count += graph[i][j]
            if count == w-1:
                wheel_count += 1
    return wheel_count

def calcWheelRamsey(i, n, m, degree, w1, w2):
    #print("vectorNumber:",i)
    coloringVector = i
    graph1, graph2, firstgraph, secondgraph, G1, G2 = calcGraph(coloringVector, m, n, degree)
    cycles = []
    cycle2 = []
    wheel1_count = 0
    wheel2_count = 0
    amount1 = 0
    amount2 = 0

    #wheel1の探索
    if w1 > 4:
        for edge in graph1:
            for node in edge:
                findCycles([node], graph1, cycles)
                for cy in cycles:
                    path = [str(node) for node in cy]
                    int_path = [int(node) for node in cy]
                    s = ",".join(path)
                    if len(int_path) == w1-1:
                        #print("cycle1:",int_path)
                        wheel1_count += searchWheel(int_path, firstgraph, n, w1)

        #print("s1:",s,"length:",len(path))

    else:
        amount1 += find_cliques_size_k(G1, w1)
        #print("amount1:",amount1)

    #wheel2の探索
    if w2 > 4:
        for edge2 in graph2:
            for node_second in edge2:
                findCycles([node_second], graph2, cycle2)
                for cy2 in cycle2:
                    path2 = [str(node_second) for node_second in cy2]
                    int_path2 = [int(node_second) for node_second in cy2]
                    s2 = ",".join(path2)
                    if len(path2) == w2-1:
                        wheel2_count += searchWheel(int_path2, secondgraph, n, w2)
                        #print("cycle2:",int_path2)
        #print("s2:",s2,"length:",len(path2))
    else:
        amount2 += find_cliques_size_k(G2,w2)
        #print("amount2:",amount2)
    """
    if wheel1_count + amount1 == 0:
        print("wheel1_count:",wheel1_count + amount1)
    if wheel2_count + amount2 == 0:
        print("wheel2_count:",wheel2_count + amount2)
    """
    #if wheel1_count + amount1 + wheel2_count + amount2 == 0:
        #print("vectorNumber:",i)
        #print("sum:", wheel1_count + amount1 + wheel2_count + amount2, "\n")
    
    return wheel1_count + amount1 + wheel2_count + amount2

def wrap_arg(param):
    return calcWheelRamsey(*param)

"""
def multi_process(sampleList):
    p = Pool(24)
    result = p.imap(wrap_arg, sampleList)
    p.close()
    return result
"""

def main():
    global graph1, graph2
    global cycles, cycle2

    n,m,w1,w2 = map(int, input("n, m, wheel1,wheel2: ").split())
    start = time.time()
    #n = 8
    #m = 1
    d = int(n / m)
    N = int(n*(n-1)/2)

    #w1 = 3
    #w2 = 3

    # Degree of the vertices
    degree = [0] * (n+1)

    #coloringVector = 5
    if m == 1:
        coloringVector = random.randint(0, 2**math.floor(n//2)-1)
        colorMax = 2**math.floor(n//2)-1
        print("colorMax:",2**math.floor(n//2)-1)
    else:
        length = m * math.floor(d//2) + int(comb(m,2))*(math.floor(d//2) + 1)
        colorMax = 2**length
        coloringVector = random.randint(0, 2**length)
        print("colorMax:",2**length)

    #coloringVector = 7
    #print("coloringVector:",coloringVector)

    paramList = [(i, n, m, degree, w1, w2) for i in range(colorMax)]
            #print("wheel1_count:",wheel1_count + amount1)
            #print("wheel2_count:",wheel2_count + amount2)
            #print("sum:", wheel1_count + amount1 + wheel2_count + amount2, "\n")
    """
    with tqdm(total=colorMax) as t:
        res = multi_process(paramList)
        print(res)
    """
    p = Pool(24)
    result = []

    with tqdm(total=colorMax) as t:
        for i in p.imap_unordered(wrap_arg, paramList):
            t.update(1)
            result.append(i)

    print(result)
    print(0 in result)
    elapsed_time = time.time() - start
    print(elapsed_time)

if __name__ == "__main__":
    main()