# DIJKSTRA ALGORITHM: Calculates the shortest path (minimum cost) between vertices in a graph.

def dijkstra(graph, start, goal):

    shortest_distance = {}
    predecessor = {}
    unseen_nodes = graph
    infinity = 99999999999
    path = []

    for node in unseen_nodes:
        shortest_distance[node] = infinity
    shortest_distance[start] = 0

    while unseen_nodes:
        
        min_node = None
        
        for node in unseen_nodes:
            if min_node is None:
                min_node = node
            elif shortest_distance[node] < shortest_distance[min_node]:
                min_node = node

        for child_node, weight in graph[min_node].items():
            if weight + shortest_distance[min_node] < shortest_distance[child_node]:
                shortest_distance[child_node] = weight + shortest_distance[min_node]
                predecessor[child_node] = min_node

        unseen_nodes.pop(min_node)

    current_node = goal
    while current_node != start:
        try:
            path.insert(0, current_node)
            current_node = predecessor[current_node]
        except KeyError:
            print('ERROR: There is no path between the input vertices!')
            break

    path.insert(0, start)
    if shortest_distance[goal] != infinity:
        print('Shortest Path: ' + str(path))
        print('Distance:', shortest_distance[goal])



graph = {'a': {'b': 10, 'c': 3},
         'b': {'c': 1, 'd': 2},
         'c': {'b': 4, 'd': 8, 'e': 2},
         'd': {'e': 7},
         'e': {'d': 9}}

dijkstra(graph, 'a', 'd')
