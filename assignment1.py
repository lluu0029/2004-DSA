import heapq
import math




# Q1
def fuse(fitmons):
    """
        Function Description:
        Returns an integer representing the fitmon with the highest cuteness score after fusing all fitmons based on the input array.

        Approach Description:
        Makes use of a memo table with a list of lists to represent the structure of a diagonal which stores the maximum cuteness score 
        obtained for each fuse combination.
        Uses a bottom-up approach to check all possible combinations by using values calculated in previous iterations of the diagonal. 
        The final element in memo[0] will be the maximum possible cuteness score.
        Eg, fitmons = [1, 2, 3]
        [1, max(1*2), max(1*max(2*3), max(1*2)*3)]
        [2,max(2*3]
        [3]

        Input:
        fitmons, a list of fitmons, where fitmons are represented as [affinity_left, cuteness_score, affinity_right].

        Output:
        integer, the cuteness score of the maximum cuteness score fitmon.

        Time Complexity: O(n^3), n = number of fitmons.
        Analysis:
        The worst possible case for the nested for loop would be when halfway through iterating fitmons.
        Iterations would be n/2, thus 'i' and 'n' would iterate over n/2 elements. 
        Since 'a' will always iterate over n-1, the worst case can be calculated as O((n-1) * (n/2) * (n/2))
        which will yield O(n^3) dominantly. 
        Adding the for loop will result in O(n^3 + n), but when only considering dominant elements: O(n^3).

        Space Complexity: O(n^2), n = number of fitmons.
        Analysis:
        The input fitmons takes n space, where n is the number of fitmons.
        The memo table creates a list for all fitmons, n space.
        For each list created, 'n' elements are stored. As such, the first fitmon will store n elements,
        second will store n-1 elements, ... the final fitmon will store 1 element. Thus the average length
        will be n/2 elements stored. 
        n number of lists of size n/2 = n * n/2 = n^2/2 = n^2.
        Thus O(n^2 + n) is the space required, resulting in O(n^2) space. 
    """
    memo = []

    # Creating memo table with the initial structure of a diagonal.
    for i in range(len(fitmons)): #O(n) time.
        row = [None] * (len(fitmons) - i)
        row[0] = fitmons[i][1]
        memo.append(row)

    iteration = 1 # Initially set to 1 where all the cuteness scores of fitmons are initially placed into the diagonal.
    for a in range(len(fitmons) - 1): # Iterates over O(n-1), n = num fitmons
        for i in range(len(fitmons) - iteration): # Worst case, iterates over O(n-1). This loop is used to iterate over the first fitmon to be fused.
            max_fuse = None
            for n in range(iteration): # Worst case, iteration = n, thus O(n). Loop used to iterate over the second possible fitmon to be fused.
                possible_fuse = int(fitmons[i+n][2] * (memo[i][n] + memo[i+n+1][iteration - n-1])) # Fusing calculation where previous values are retrieved from the memo table.
                if max_fuse is None: # Initial case where max_fuse is set to None.
                    max_fuse = possible_fuse
                else: # Checking whether the possible_fuse is greater than the currently stored max_fuse, reassigning max_fuse if required.
                    if possible_fuse > max_fuse:
                        max_fuse = possible_fuse
            memo[i][iteration] = max_fuse # Storing the maximum possible fuse into the memo table.
        iteration += 1

    return memo[0][-1] # Returns the final fitmon with the maximum cuteness score.



# Q2
class PriorityQueue:
    def __init__(self):
        """
            Function Description: 
            Initialising the heap as an empty list.

            Time and Space Complexity: O(1)
        """
        self.heap = []

    def push(self, item, key):
        """
            Function Description:
            Pushes an item into the heap.

            Input: item, key: integers

            Time Complexity: O(log(n)), n is the number of elements in the heap.
            Analysis: Pushing an item into a heap involves appending to the end of the list (self.heap)
            and may require swaps to be made to restore the correct order. 

            Space Complexity: O(1)
            Analysis: This function only appends a single element, a tuple of (item, key).
        """
        heapq.heappush(self.heap, (item, key))

    def pop_min(self):
        """
            Function Description:
            Removes and returns the smallest element based on 'key' from the heap.

            Output: Tuple, (item, key)

            Time Complexity: O(log(n)), n is the number of elements in the heap.
            Analysis: Removing the smallest element must require swaps to be made to the
            resulting heap to restore the correct 'heap order'.

            Space Complexity: O(1)
            Analysis: No additional space required.
        """
        return heapq.heappop(self.heap)

    def empty(self):
        """
            Function Description:
            Returns whether len(self.heap) is 0.

            Output: True/False

            Time Complexity: O(1)
            Analysis: Constant time to check length of a list.

            Space Complexity: O(1)
            Analysis: No additional space required.
        """
        return len(self.heap) == 0


class TreeMap():
    def __init__(self, roads, solulus):
        """
            Function Description:
            Creates a graph with roads representing edges and trees representing vertices.
            Uses an adjacency list to represent edges and their weights associated with vertex.

            Approach Description:
            The problem can be modelled as two connected 'mirror graphs', with the same roads and trees,
            however only the solulu trees which are found on the '1st' part of the graph are connected to the 
            '2nd' part of the graph which contains the exits. 
            This forms a single graph with 2R (roads) and 2T (trees), and ensures that a solulu
            must be visited to reach the '2nd' part of the graph where the exits are. 
            This allows the use of a single Dijkstra algorithm to find the shortest path to all other vertices,
            meeting the O(Rlog(T)) time complexity for escape.

            Input:
            roads, a list of tuples (u, v, w) containing non-negative integers.
            solulus, a list of tuples (x, y, z) containing non-negative integers.

            Time Complexity: 
            O(R + T), R is the number of elements in roads, T is the number of unique vertices (trees) found in roads.
            Analysis: 
            Summing the complexities in the function, we get O(R) + O(2T) + O(R) + O(T) = O(2R + 3T).
            As a result, when focusing on dominant terms we get O(R + T).

            Space Complexity: O(T + R)
            Analysis:
            The adjacency list which stores vertices (trees) and their edges (roads) creates a list for each tree O(T),
            and stores corresponding roads for each tree in each list O(R).
            When adding the initial input lists for roads and solulus, O(T + R + T + R) = O(2T + 2R) = O(T + R)
        """
        # Greatest value for a vertex used to determine the number of unique vertices (trees).
        vertices = None
        for road in roads: # O(R)
            if vertices is None:
                vertices = max(road[0], road[1])
            else:
                vertices = max(vertices, road[0], road[1])

        self.num_trees = 2*(vertices + 1)

        # Adjacency list to store roads and time (edges and weights) for each tree (vertex).
        self.adj_list = [[] for _ in range(self.num_trees)] # Time: O(2T) Space: O(T + R)

        # Adds road (edge) to the adjacency list.
        for road in roads: # O(R)
            # Adding edge to '1st' graph.
            self.add_edge(road[0], road[1], road[2])
            
            # Add the same edge in the '2nd' graph.
            self.add_edge(road[0] + (self.num_trees//2), road[1] + (self.num_trees//2), road[2])

        # Adds teleport of solulus as an edge from the '1st' graph to the '2nd' graph.
        for solulu in solulus: # O(T)
            self.add_edge(solulu[0], solulu[2] + (self.num_trees//2), solulu[1])

    def add_edge(self, u, v, weight):
        """
            Function Description:
            Adds an edge to vertex tree.

            Input:
            u, integer representing the tree the road is being added to.
            v, integer for the ending tree for the road.
            weight, integer representing the time taken travel a road or destroy a solulu.

            Time and Space Complexity: O(1)
            Analysis: append is a constant operation. 
        """
        self.adj_list[u].append((v, weight))

    def escape(self, start, exits):
        """
            Function Description:
            This function returns the shortest possible time and path to reach one of the exit trees found in exits.
            
            Approach Description:
            An implementation of Dijkstra's algorithm is used to return a list of shortest time from 'start'
            to all other trees, allowing for the shortest time to one of the exits to be derived as well as a 
            list of predecessor trees which is used to generate the path.

            Input:
            start, integer representing the starting tree for the path.
            exits, list of integers representing the possible trees which are possible exits.

            Output:
            None or Tuple: (time, path) where time is an integer representing the time taken to exit the forest. 
            path is a list of integers representing the order of trees to take to exit the forest.
            None is returned when there is no path to reach the exits.

            Time Complexity: O(Rlog(T)), R is the number of roads, T is the number of trees.
            Analysis:
            The implementation of function dijkstra yields a time complexity of O(E log(V)), where in this case E is the number
            of roads, and T is the number of trees. Based on the adjacency list used to represent TreeMap,
            E = 2R and V = 2T given the implementation of 2 'mirror' graphs. Thus this can be represented as O(R log(T)).

            Auxiliary Space Complexity: O(T)
            The use of dijkstra requires O(v) space (see function dijkstra), in this case v is 2T.
            For list exits, the worst case, all trees are exits and thus can have a maximum size of 2T.
            For list path, the worst case all trees are visited, and thus can have a maximum size of 2T.
            Thus O(2T + 2T + 2T) = O(6T) = O(T) auxiliary space required.

        """
        dist, pred = dijkstra(self.adj_list, start)

        min_time = None
        min_exit = None

        # Reindex the index of exit trees to represent them on the '2nd' part of the graph.
        exit_list = [x + (self.num_trees//2) for x in exits]

        # Finding shortest distance to exit.
        for exit_index in exit_list:
            if min_time is None:
                min_time = dist[exit_index]
                min_exit = exit_index
            else:
                if min_time > dist[exit_index]:
                    min_time = dist[exit_index]
                    min_exit = exit_index

        path = []
        current_tree = min_exit

        # Iterating through list pred to find shortest path.
        while current_tree is not None:
            if current_tree > ((self.num_trees // 2) - 1):
                path.append(current_tree - (self.num_trees // 2))
            else:
                if path[-1] != current_tree:
                    path.append(current_tree)
            current_tree = pred[current_tree]

        # If min_time is infinity, there is no path, thus return None.
        if math.isinf(min_time):
            return None
        else:
            return (min_time, path[::-1])


def dijkstra(graph, source):
    """
        Function Description:
        Implementation of Dijkstra's algorithm using a min-heap priority queue to 
        find the shortest path to all other vertices in a graph. 

        Input:
        graph, list of lists containing integers to represent an adjacency list of a graph.
        source, integer representing the starting vertex.

        Output:
        dist, list of integers representing the shortest distance to other vertices based on index.
        pred, list representing the predecessor vertex for the shortest path. Can be None or an integer.

        Time Complexity: O(Elog(V)), E is the number of edges, V is the number of vertices.
        Analysis: 
        The most dominant part of the function would be within the while loop. In the worst case, all vertices are connected,
        thus all vertices and their edges will be iterated through. As a result, the number of iterations can be E (total number of edges).
        Each loop, Q.pop_min() and Q.push() is performed, thus complexity can be represented as O(E(log(n) + log(n))), resulting in O(E log(n))
        where n is the number of vertices.

        Auxiliary Space Complexity: O(v), where v is the number of vertices.
        Analysis:
        v represents the number of vertices. dist, pred both have size O(v). In the worst case where all vertices need to be iterated,
        PriorityQueue will require O(v) space for the heap. As a result, O(3v) = O(v) auxiliary space is required.
    """
    n = len(graph)  # Number of vertices. O(1)
    dist = [float('inf')] * n # List storing distances. O(n)
    pred = [None] * n  # List storing predecessor vertices, initially None. O(n)

    dist[source] = 0
    Q = PriorityQueue() # Initialising priority queue. O(1)
    Q.push(source, key=0) # Pushing the source element into the priority queue. O(log(n))

    while not Q.empty(): # Q.empty() takes O(1) time. Worst case will iterate over all vertices.
        u, key = Q.pop_min() # O(log(n)) time
        if dist[u] == key:
            for v, weight in graph[u]: # Iterate through all edges for vertex u. Worst case, when iterating over all vertices, will iterate over all edges.
                if dist[v] > dist[u] + weight:
                    dist[v] = dist[u] + weight
                    pred[v] = u
                    Q.push(v, key=dist[v]) # push takes O(log(n)) time.

    return dist, pred