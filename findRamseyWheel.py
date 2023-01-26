import itertools
import sys
import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support
from matplotlib.style import available

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
    for j in range(n):
        secondColorgraph[j][j] = 0

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


def find_cliques_size_k(G, k):
    i = 0
    for clique in nx.find_cliques(G):
        if len(clique) == k:
            i += 1
        elif len(clique) > k:
            i += int(comb(len(clique), k))
    
    if i == 0:
        return False
    else:
        return True

def get_unique_list(vertex_cycleList):
    seen = []
    return [i for i in vertex_cycleList if i not in seen and not seen.append(i)]

def nearList(vertexarray):
    return np.where(vertexarray == 1)

def findingWheel(G,size, graphmatrix):
    if size <= 4:
        return find_cliques_size_k(G,size)

    else:
        usedSet = np.zeros(1)
        vertexSet = np.array([i for i, x in enumerate(graphmatrix[0]) if x == 1])
        #print("before UsedSet:",usedSet)
        #print("vertexSet:",vertexSet)
        near0 = nearList(graphmatrix[0])
        #print("availableNear0:",near0)
        for v in vertexSet:
            #print("v:",v,"\n")
            usedSet = np.union1d(usedSet, v)
            #print("usedList:",usedSet)
            #print("vより大きいインデックスの抽出")
            vertexCandidate = np.array([vertexSet[i] for i, x in enumerate(vertexSet) if x > v])
            #print(vertexCandidate)
            availableSet = np.intersect1d(vertexSet, vertexCandidate)
            #print("availableSet:",availableSet)
            if cycleNear0(size-1, usedSet, availableSet, v, graphmatrix):
                return True
        return False

def cycleNear0(length, used, availableList, target, graphmatrix):
    #print("now length;",length)
    if length == 2:
        nearTarget = nearList(graphmatrix[target])
        #print("nearTarget:",nearTarget[0])
        #print("last Available:",availableList)
        #print(np.intersect1d(nearTarget[0], availableList))
        if len(np.intersect1d(nearTarget[0], availableList)) > 0:
            return True
        else:
            return False

    else:
        for v in availableList:
            #print("available v:",v)
            newused = np.union1d(used, v)
            #print("newused:",newused)

            nearV = nearList(graphmatrix[v])
            #print("nearV",nearV[0])
            near0 = nearList(graphmatrix[0])
            #print("near0:",near0[0])
            newAvailable = np.intersect1d(nearV[0], near0[0])
            newAvailable = np.setdiff1d(newAvailable, used)
            #print("newAvailable:",newAvailable)
            if cycleNear0(length-1, newused, newAvailable, target, graphmatrix):
                return True

def calcWheelRamsey(i, n, m, degree, w1, w2):
    # print("vectorNumber:",i)
    coloringVector = i
    firstgraph, secondgraph, G1, G2 = calcGraph(coloringVector, m, n, degree)

    #Wheelの探索
    firstColorResult = findingWheel(G1, w1, firstgraph)
    #print("result:",firstColorResult,"\n")
    secondColorResult = findingWheel(G2, w2, secondgraph)
    allResult = firstColorResult or secondColorResult
    #return [i,firstColorResult,secondColorResult, allResult]
    return [i, allResult]

def wrap_arg(param):
    return calcWheelRamsey(*param)


def main():
    p = Pool(24)
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

    """
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
    """
    paramList = [(i, n, m, degree, w1, w2) for i in range(colorMax+1)]
    result = []
    with tqdm(total=colorMax) as t:
        for i in p.imap_unordered(wrap_arg, paramList):
            t.update(1)
            result.append(i)

    #print(result)
    count = 0
    for j in range(len(result)):
        if result[j][1] == False:
            count +=1

    print("0 monochromatic copy:",count)

    #ファイル名
    filename = str(n) + "_" + str(m) + "_W" + str(w1) + "_W" + str(w2) + ".txt"
    #出力

    res = '\n'.join(' '.join(map(str, x)) for x in result)
    with open(filename, 'w')as f:
        f.write(res)

    elapsed_time = time.time() - start
    print(elapsed_time)


if __name__ == "__main__":
    main()
