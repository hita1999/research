import time
from multiprocessing import Pool
from multiprocessing.spawn import freeze_support

import numpy as np
from numpy.random import *
from scipy.linalg import circulant
from tqdm import tqdm
from memory_profiler import profile
import sys

n = 6
N = int(n*(n-1)/2)
thr = 24

def integer_to_binary(p,n):
    bistr = bin(p)
    bistr = bistr[2:]
    bistr = '0'*(n-len(bistr)) + bistr
    bistr = [int(s) for s in bistr]
    return bistr

def calcGraph(i):
    A = [[0 for j in range(n)] for i in range(n)]
    bistr = integer_to_binary(i,N)
    for j in range(0,n):
        for k in range(j+1,n):
            b = bistr.pop()
            A[j][k] = b
            A[k][j] = b
    v = np.array(A)

def main():
    print(N)
    count = 0
    start = time.time()
    with Pool(thr) as p:
        while count < 2**N:
            pro = p.map(calcGraph, tqdm(range(count, count+2**28)))
            count += 2**28
        #pro = p.map(calcGraph, range(2**N))
    #for i in tqdm(range(2**N)):
        #calcGraph(i)
    elapsed_time = time.time() - start
    print(elapsed_time)

if __name__ == "__main__":
    main()