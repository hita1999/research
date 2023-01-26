import numpy as np

g = np.loadtxt("resultMatrix.txt")

visit = [0] * len(g)  # 訪問済 or 未訪問

def DFS(u):
    visit[u] = 1
    adj = g[u]
    for i in range(0, len(adj)):
        if adj[i] == 1 and visit[i] == 0:
            print('edge=' + str([u, i]))
            DFS(i)

DFS(0)