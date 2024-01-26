# EXERCISE Given a protein network, propose an algorithm to check if it is connected.
# To check if the elements of a network are connected, we can use the DFS algorithm.

def create_graph(edges, n_nodes):
    
    graph = [[] for _ in range(n_nodes)]
    for i in range(0, len(edges)):
        graph[edges[i][0]].append(edges[i][1])
    
    return graph


# Depth First Search
def dfs(visited, node):

    visited[node] = True
    visited_nodes.append(node)

    for neighbour in graph[node]:
        if not visited[neighbour]:
            dfs(visited, neighbour)


edges = [(0, 4), (1, 0), (1, 2), (2, 1), (2, 4), (3, 1), (3, 2), (4, 3)]
n_nodes = 5

graph = create_graph(edges, n_nodes)
visited = [False] * n_nodes
visited_nodes = []

dfs(visited, 0)

print(visited)
print(visited_nodes)

connected = True
for n in visited:
    if n == False:
        connected = False

print("\nRede conexa:", connected)

"""
# Para ver se a rede é fortemente conectada, repete-se o mesmo processo com as direções invertidas
edges = [(j, i) for i in range(n_nodes) for j in graph[i]]
graph = create_graph(edges, n_nodes)
visited = [False] * n_nodes

dfs(visited, 0) 
"""

