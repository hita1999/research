import networkx as nx
from networkx.classes.function import degree
from networkx.algorithms.cuts import edge_expansion
import numpy as np
import matplotlib.pyplot as plt
import itertools


def matrixToedges(adj_mat, degree, color):
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
        edges.append([i2_component[k], j2_component[k], {"color": color}])

    size = len(edges)

    for i in range(size):
        degree[edges[i][0]] += 1
        degree[edges[i][1]] += 1
    return edges

g = np.loadtxt("resultMatrix.txt")
g2 = np.loadtxt("resultMatrix2.txt")

G = nx.Graph()
#G = nx.from_numpy_matrix(g)
G2 = nx.from_numpy_matrix(g2)

n = int(input("the number of vertices "))
"""
verList = [i for i in range(n)]
combList = []
for cycle in itertools.combinations(verList, 4):
    combList.append(list(cycle))
"""
# Degree of the vertices
degree = [0] * (n+1)

edgeList1 = matrixToedges(g, degree, "red")
edgeList2 = matrixToedges(g2, degree, "blue")
#print(*edgeList2, sep="\n")

G.add_edges_from(edgeList1)
G.add_edges_from(edgeList2)

edge_color = [edge["color"] for edge in G.edges.values()]
nx.draw(G, pos=nx.circular_layout(G), with_labels=True, edge_color=edge_color)
#nx.draw(G2, with_labels=True)

plt.show()
