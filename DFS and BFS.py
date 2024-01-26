######################## Depth First Search (DFS) ########################

# Depth-First Search Method: Explore the graph/tree (in depth) as far as possible,
# if not, backtrack and alternatively try other paths.
# Does not guarantee the shortest path or the path with the fewest steps first.

graph = {'A': ['B', 'C'],
         'B': ['D', 'E'],
         'C': ['F'],
         'D': [],
         'E': ['F'],
         'F': []
         }

visited = []


def dfs(visited, graph, node):

    if node not in visited:
        visited.append(node)
        for neighbour in graph[node]:
            dfs(visited, graph, neighbour)


dfs(visited, graph, 'A')
print("DFS:", visited)

######################## Breadth First Search (BFS) ########################

# Breadth-First Search Method: Starting from a node, all adjacent nodes are explored, and then
# nodes accessible through the adjacent nodes (next level) are explored, and so on. It is necessary
# to store all paths that can still be expanded. Guarantees that the path with the fewest steps is
# found first.


graph = {'A': ['B', 'C'],
         'B': ['D', 'E'],
         'C': ['F'],
         'D': [],
         'E': ['F'],
         'F': []
         }

visited = []
queue = []


def bfs(visited, graph, node):
    visited.append(node)
    queue.append(node)

    while queue:
        s = queue.pop(0)

        for neighbour in graph[s]:
            if neighbour not in visited:
                visited.append(neighbour)
                queue.append(neighbour)


bfs(visited, graph, 'A')
print("\nBFS:", visited)
