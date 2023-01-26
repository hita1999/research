import math
import random
import time
from multiprocessing import Pool

import networkx as nx
import numpy as np
from numpy.random import *
from scipy.linalg import circulant
from scipy.special import comb
from tqdm import tqdm


#10進数から2進数彩色ベクトルへ変換
def integer_to_binary(p, n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    return bistr


#グラフの行列作成
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

    secondColorgraph = np.ones((n, n), dtype=int) - graphMat
    for j in range(n):
        secondColorgraph[j][j] = 0

    edgesList1 = matrixToedges(graphMat, degree)
    G1.add_edges_from(edgesList1)

    edgesList2 = matrixToedges(secondColorgraph, degree)
    G2.add_edges_from(edgesList2)

    return graphMat, secondColorgraph, G1, G2

#行列から2頂点間の辺情報を返す
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

#クリークの探索
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

#入力された頂点の近傍を返す
def nearList(vertexarray):
    return np.where(vertexarray == 1)

#車輪グラフの探索
def findingWheel(G,size, graphmatrix):
    if size <= 4:
        return find_cliques_size_k(G,size)

    else:
        usedSet = np.zeros(1)
        vertexSet = np.array([i for i, x in enumerate(graphmatrix[0]) if x == 1])
        for v in vertexSet:
            usedSet = np.union1d(usedSet, v)
            vertexCandidate = np.array([vertexSet[i] for i, x in enumerate(vertexSet) if x > v])
            availableSet = np.intersect1d(vertexSet, vertexCandidate)
            if cycleNear0(size-1, usedSet, availableSet, v, graphmatrix):
                return True
        return False

#閉路の探索
def cycleNear0(length, used, availableList, target, graphmatrix):
    if length == 2:
        nearTarget = nearList(graphmatrix[target])
        if len(np.intersect1d(nearTarget[0], availableList)) > 0:
            return True
        else:
            return False

    else:
        for v in availableList:
            newused = np.union1d(used, v)
            nearV = nearList(graphmatrix[v])
            near0 = nearList(graphmatrix[0])
            newAvailable = np.intersect1d(nearV[0], near0[0])
            newAvailable = np.setdiff1d(newAvailable, used)
            if cycleNear0(length-1, newused, newAvailable, target, graphmatrix):
                return True

#車輪グラフの探索結果を返す
def calcWheelRamsey(i, n, m, degree, w1, w2):
    coloringVector = i
    firstgraph, secondgraph, G1, G2 = calcGraph(coloringVector, m, n, degree)

    #Wheelの探索
    firstColorResult = findingWheel(G1, w1, firstgraph)
    secondColorResult = findingWheel(G2, w2, secondgraph)
    allResult = firstColorResult or secondColorResult
    return [i, allResult]

def wrap_arg(param):
    return calcWheelRamsey(*param)


def main():
    p = Pool(24)

    n, m, w1, w2 = map(int, input("n, m, wheel1,wheel2: ").split())
    start = time.time()
    d = int(n / m)
    degree = [0] * (n+1)

    if m == 1:
        colorMax = 2**math.floor(n//2)-1
        print("colorMax:", 2**math.floor(n//2)-1)
    else:
        length = m * math.floor(d//2) + int(comb(m, 2))*(math.floor(d//2) + 1)
        colorMax = 2**length
        print("colorMax:", 2**length)

    paramList = [(i, n, m, degree, w1, w2) for i in range(colorMax+1)]
    result = []
    with tqdm(total=colorMax) as t:
        for i in p.imap_unordered(wrap_arg, paramList):
            t.update(1)
            result.append(i)

    elapsed_time = time.time() - start
    count = 0
    for j in range(len(result)):
        if result[j][1] == False:
            count +=1

    print("0 monochromatic copy:",count)

    filename = str(n) + "_" + str(m) + "_W" + str(w1) + "_W" + str(w2) + ".txt"
    #出力
    res = '\n'.join(' '.join(map(str, x)) for x in result)
    with open(filename, 'w')as f:
        f.write(res)

    print(elapsed_time)


if __name__ == "__main__":
    main()
