# Python3 implementation of the approach
import numpy as np

MAX = 6

# Stores the vertices
store = [0]* MAX

# Graph
graph = np.zeros((MAX, MAX))

# Degree of the vertices
d = [0] * MAX

# Function to check if the given set of vertices
# in store array is a clique or not
def is_clique(b) :

	# Run a loop for all the set of edges
	# for the select vertex
	for i in range(1, b) :
		for j in range(i + 1, b) :

			# If any edge is missing
			if (graph[store[i]][store[j]] == 0) :
				return False
	
	return True

# Function to print the clique
def print_cli(n) :

	for i in range(1, n) :
		print(store[i], end = " ")
	print(",", end=" ")

# Function to find all the cliques of size s
def findCliques(i, l, s) :

	# Check if any vertices from i+1 can be inserted
	for j in range( i + 1, n -(s - l) + 1) :

		# If the degree of the graph is sufficient
		if (d[j] >= s - 1) :

			# Add the vertex to store
			store[l] = j

			# If the graph is not a clique of size k
			# then it cannot be a clique
			# by adding another edge
			if (is_clique(l + 1)) :

				# If the length of the clique is
				# still less than the desired size
				if (l < s) :

					# Recursion to add vertices
					findCliques(j, l + 1, s)

				# Size is met
				else :
					print_cli(l + 1)

# Driver code
if __name__ == "__main__" :

	edges = [ [ 1, 2 ],
			[ 2, 3 ],
			[ 3, 1 ],
			[ 4, 3 ],
			[ 4, 5 ],
			[ 5, 3 ] ]
	k = 3
	size = len(edges)
	n = 5

	for i in range(size) :
		graph[edges[i][0]][edges[i][1]] = 1
		graph[edges[i][1]][edges[i][0]] = 1
		d[edges[i][0]] += 1
		d[edges[i][1]] += 1

print(graph)
print(len(graph))
findCliques(0, 1, k)

# This code is contributed by AnkitRai01
