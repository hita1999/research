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


def integer_to_binary(p, n):
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

    # make binary Vector
    if m == 1:
        bistr = integer_to_binary(i, math.floor(n//2))
    else:
        length = m * math.floor(d//2) + int(comb(m, 2))*(math.floor(d//2) + 1)
        bistr = integer_to_binary(i, length)

    # make Circulant Matrix
    cirIdxList = [1]
    cirIdxNum = m + 1
    for j in range(1, m):
        cirIdxList.append(cirIdxNum)
        cirIdxNum = cirIdxNum + (m-j)

    getIdxNum = 0
    putIdxNum = 0
    cirMatList = []
    if m != 1:
        choose_Number_of_Color = int(
            (len(bistr) - len(cirIdxList)*math.floor(d/2)) / (m*(m-1)/2))

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
                idx += 1

    #print(graphMat)
    secondColorgraph = np.ones((n, n), dtype=int) - graphMat

    edgesList1 = matrixToedges(graphMat, degree)
    G1.add_edges_from(edgesList1)

    #print(secondColorgraph)
    edgesList2 = matrixToedges(secondColorgraph, degree)
    G2.add_edges_from(edgesList2)

    return graphMat, secondColorgraph, G1, G2


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


def DFS(graph, marked, w, n, vert, start, cycle_list):
    # mark the vertex vert as visited
    marked[vert] = True
    # if the path of length (n-1) is found
    if w == 0:
        # mark vert as un-visited to make
        # it usable again.
        marked[vert] = False

        # Check if vertex vert can end with
        # vertex start
        if graph[vert][start] == 1:
            cycle_list.append(vert)
            cycle_list.append(start)
            cycle_list.append("start")
            cycle_list.append(start)
            return cycle_list
        else:
            return cycle_list

    # For searching every possible path of
    # length (n-1)
    for i in range(n):
        if marked[i] == False and graph[vert][i] == 1:

            # DFS for searching path by decreasing
            # length by 1
            cycle_list.append(vert)
            cycle_list = DFS(graph, marked, w-1, n, i, start, cycle_list)

    # marking vert as unvisited to make it
    # usable again.
    marked[vert] = False
    cycle_list.append("end")
    return cycle_list


def find_cliques_size_k(G, k):
    i = 0
    for clique in nx.find_cliques(G):
        if len(clique) == k:
            i += 1
        elif len(clique) > k:
            i += int(comb(len(clique), k))
    return i

def get_unique_list(vertex_cycleList):
    seen = []
    return [i for i in vertex_cycleList if i not in seen and not seen.append(i)]


def searchWheel(cycleList, graph, w, n):
    wheel_count = 0
    cycleCandidate = []
    vertex_cycleList = []
    sortedCandidate = []
    vertex_all = []
    #print("cycleList:",cycleList)
    #print(cycleList.count("start"))
   
    for j in range(len(cycleList)):
        if cycleList[j] == "start" and j-w > -1:
            is_duplicate = cycleList[j-w:j]
            if is_duplicate[0] == is_duplicate[-1] and len(is_duplicate)-1 == len(set(is_duplicate)) and not "end" in is_duplicate:
                #cycleCandidate += cycleList[j-w:j]
                #print(is_duplicate)
                vertex_cycleList.append(cycleList[j-w+1:j])
    """
    for j in range(0,len(cycleCandidate),w):
        vertex_cycle = cycleCandidate[j:j+w]
        print("vertex_cycle:",vertex_cycle)
    """
    del cycleList
    for j in range(0,len(cycleCandidate),w):
        #vertex_all = [int(i) for i in range(n)]
        sortedCandidate = cycleCandidate[j:j+w-1]
        vertex_cycleList.append(sortedCandidate)
        #print(sortedCandidate)

    for vertex in vertex_cycleList:
        vertex_all = [int(i) for i in range(n)]
        for i in vertex:    
            vertex_all.remove(i)
        #print("vertex_candidate:",vertex_all)

        for k in vertex_all:
            count = 0
            #print("vertex_all:",vertex_all)
            #print(graph[k])

            for l in vertex:
                count += graph[k][l]
                if count == w-1:
                    #print("one point:",k)
                    #print("wheel:",vertex)
                    wheel_count += 1
    return wheel_count

# Counts cycles of length
# N in an undirected
# and connected graph.


def countCycles(graph, w, n):

    # all vertex are marked un-visited initially.
    marked = [False] * n

    # Searching for cycle by using v-n+1 vertices
    cycle_list = []
    for j in range(n-(w-1)):
        cycle_list = DFS(graph, marked, w-2, n, j, j, cycle_list)

        # ith vertex is marked as visited and
        # will not be visited again.
        marked[j] = True

    #print("cycleList:",cycle_list)
    return cycle_list


def calcWheelRamsey(i, n, m, degree, w1, w2):
    # print("vectorNumber:",i)
    coloringVector = i
    firstgraph, secondgraph, G1, G2 = calcGraph(coloringVector, m, n, degree)
    cycle1 = []
    cycle2 = []
    wheel1_count = 0
    wheel2_count = 0
    amount1 = 0
    amount2 = 0

    # wheel1の探索
    if w1 > 4:
        cycle1 = countCycles(firstgraph, w1, n)
        wheel1_count += int(searchWheel(cycle1, firstgraph, w1, n)/2)
        #print("wheel1_count:",wheel1_count)

    else:
        amount1 += find_cliques_size_k(G1, w1)
        # print("amount1:",amount1)

    # wheel2の探索
    if w2 > 4:
        #cycle2 = countCycles(secondgraph, w2, n)
        cycle2 = countCycles(secondgraph, w2, n)
        wheel2_count += int(searchWheel(cycle2, secondgraph, w2, n)/2)

    else:
        amount2 += find_cliques_size_k(G2, w2)
        # print("amount2:",amount2)
    """
    if wheel1_count + amount1 == 0:
        print("wheel1_count:",wheel1_count + amount1)
    if wheel2_count + amount2 == 0:
        print("wheel2_count:",wheel2_count + amount2)
    """
    # if wheel1_count + amount1 + wheel2_count + amount2 == 0:
    # print("vectorNumber:",i)
    #print("sum:", wheel1_count + amount1 + wheel2_count + amount2, "\n")

    #return [i,wheel1_count + amount1 + wheel2_count + amount2, wheel1_count, amount1, wheel2_count, amount2]
    sum = wheel1_count + amount1 + wheel2_count + amount2
    return [i,sum]

def wrap_arg(param):
    return calcWheelRamsey(*param)


def main():
    p = Pool(1)
    #global graph1, graph2
    #global cycles, cycle2

    n, m, w1, w2 = map(int, input("n, m, wheel1,wheel2: ").split())
    start = time.time()
    #n = 8
    #m = 1
    d = int(n / m)
    N = int(n*(n-1)/2)

    #w1 = 3
    #w2 = 3

    # Degree of the vertices
    degree = [0] * (n+1)

    #coloringVector = 7
    if m == 1:
        coloringVector = random.randint(0, 2**math.floor(n//2)-1)
        colorMax = 2**math.floor(n//2)-1
        print("colorMax:", 2**math.floor(n//2)-1)
    else:
        length = m * math.floor(d//2) + int(comb(m, 2))*(math.floor(d//2) + 1)
        colorMax = 2**length
        coloringVector = random.randint(0, 2**length)
        print("colorMax:", 2**length)

    paramList = [(i, n, m, degree, w1, w2) for i in range(int((colorMax+1)/3))]
    paramList2 = [(i, n, m, degree, w1, w2) for i in range(int((colorMax+1)/3), int(2*(colorMax+1)/3))]
    paramList3 = [(i, n, m, degree, w1, w2) for i in range(int(2*(colorMax+1)/3), colorMax+1)]

    result = []

    with tqdm(total=int(colorMax/3)) as t:
        for i in p.imap_unordered(wrap_arg, paramList):
            t.update(1)
            result.append(i)

    with tqdm(total=int(colorMax/3)) as t:
        for i in p.imap_unordered(wrap_arg, paramList2):
            t.update(1)
            result.append(i)
    
    with tqdm(total=int(colorMax/3)) as t:
        for i in p.imap_unordered(wrap_arg, paramList3):
            t.update(1)
            result.append(i)
    #print(result)
    npresult = np.array(result)
    #if 0 in npresult[:,1]:
    #print(*result, sep='\n')
    print(0 in npresult[:,1])
    count0 = npresult[np.any(npresult == 0, axis=1)]
    print("0monochromatic copy:",len(count0))
    #ファイル名
    filename = str(n) + "_" + str(m) + "_W" + str(w1) + "_W" + str(w2) + ".txt"
    #出力

    res = '\n'.join('\t'.join(map(str, x)) for x in result)
    with open(filename, 'w')as f:
        f.write(res)

    elapsed_time = time.time() - start
    print(elapsed_time)


if __name__ == "__main__":
    main()
"""
# main :
graph = [[0, 1, 0, 1, 0],
		[1 ,0 ,1 ,0, 1],
		[0, 1, 0, 1, 0],
		[1, 0, 1, 0, 1],
		[0, 1, 0, 1, 0]]
		
n = 4
print("Total cycles of length ",n," are ",countCycles(graph, n))
"""
# this code is contributed by Shivani Ghughtyal
