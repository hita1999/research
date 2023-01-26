import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support

import numpy as np
from numpy.random import *
from scipy.linalg import circulant
from tqdm import tqdm
from memory_profiler import profile
import sys

n = 64
N = int(n*(n-1)/2)
thr = 1

"""
def integer_to_binary(p,n):
    bistr = bin(p)[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    print(bistr)
    return bistr
"""
def integer_to_binary(p,n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    bistr = np.array(bistr)
    #print(bistr)
    return bistr

def calcGraph(i):
    A = np.zeros((n,n), dtype=np.int64)
    bistr = integer_to_binary(i,N)
    print(bistr)
    """
    for j in range(0,n):
        for k in range(j+1,n):
            b = bistr[0]
            bistr = np.delete(bistr, 0)
            A[j,k] = b
            A[k,j] = b
    """
    v = np.tri(bistr, 1)
    print(v)
    #print(type(A))
    return A

def findClique(subset, k1, k2):
    for j in range(0,n):
        for k in range(j+1,n):
            k.append()

def main():
    print(N)
    start = time.time()
    with Pool(thr) as p:
        pro = p.map(calcGraph, tqdm(range(2**2)))
    #findClique(pro,3,3)
        #pro = p.map(calcGraph, range(2**N))
    #for i in tqdm(range(2**N)):
        #calcGraph(i)
    elapsed_time = time.time() - start
    print(elapsed_time)

if __name__ == "__main__":
    main()