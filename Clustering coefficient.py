# AVERAGE CLUSTERING COEFFICIENT

def clustering_coefficient(node, node_neighbours):

    K = len(node_neighbours) # Number of neighbors of the node

    N = 0
    edges = []  

    for key, value in graph.items():

        if key != node and key in node_neighbours:
            for n in value:
                if n != node and [key, n] not in edges and [n, key] not in edges:
                    edges.append([key, n])
                    N = N + 1  # Number of connections among the neighbors of the node
    
    C = 2*N/(K*(K-1))

    return C

# Dictionary with the connections between different nodes in the network
graph = {'a': ['b', 'c', 'd', 'e'],
         'b': ['a', 'c'],
         'c': ['a', 'b', 'd'],
         'd': ['a', 'c', 'e'],
         'e': ['a', 'd'],
         'f': ['j', 'k']}

average_C = 0

for node, node_neighbours in graph.items():
    C = clustering_coefficient(node, node_neighbours)
    average_C = average_C + C
    print("\nNode:", node, "\nClustering coefficinent:", C)  

average_C = average_C/ len(graph)
print("\nAverage clustering coefficinent:", round(average_C,2))    
