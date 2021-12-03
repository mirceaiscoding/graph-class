#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;
ifstream fin("royfloyd.in");
ofstream fout("royfloyd.out");

#define NO_PATH -1
#define NO_PARENT_NODE -1
#define TASK_NUMBER_UNITE_SETS 1
#define TASK_NUMBER_QUERY_SAME_SET 2
#define NO_EDGE 0
#define MAX_DISTANCE 1000000000

/**
 * @brief Graph class that implements algorithms
 * 
 */
class Graph
{
protected:
    // Number of nodes in the Graph
    int numberOfNodes;

    // If the Graph is oriented
    bool isOriented;

    // Edges in the Graph
    vector<vector<int> > edges;

    /**
     * @brief Does a depth first search and marks the visited nodes
     * 
     * @param node node from which to start the search
     * @param isVisited array that knows whether a node is visited
     */
    void DFS(int node, bool isVisited[]);

    /**
     * @brief Get a list of biconnected components (lists of nodes) in `biconnectedComponents`
     * 
     * @param node current node
     * @param currentDepth current depth
     * @param parentNode parent node of current node
     * @param biconnectedComponents the biconnected components will be stored here
     * @param depth depth of nodes (distance from root)
     * @param low minimum level a node can reach (without going back through parent nodes)
     * @param isVisited if a node is visited or not
     * @param visitedNodes order in which nodes are visited in the DFS
     */
    void findBiconnectedComponents(int node, int currentDepth, int parentNode, vector<vector<int> > &biconnectedComponents,
                                   int depth[], int low[], bool isVisited[], stack<int> &visitedNodes);

    /**
     * @brief Get a list of strongly connected components (lists of nodes) in `stronglyConnectedComponents`
     * 
     * @param node current node
     * @param index smallest unused index
     * @param parentNode parent node of current node
     * @param stronglyConnectedComponents the strongly connected components will be stored here
     * @param depth depth of nodes (distance from root)
     * @param low minimum level a node can reach (without going back through parent nodes)
     * @param isVisited if a node is visited or not
     * @param inStack if the node is in the visitedNodes stack
     * @param visitedNodes order in which nodes are visited in the DFS
     */
    void findStronglyConnectedComponents(int node, int &index, int parentNode, vector<vector<int> > &stronglyConnectedComponents,
                                         int depth[], int low[], bool isVisited[], bool inStack[], stack<int> &visitedNodes);

    /**
     * @brief Get a list of critical edges (pair of nodes) in `criticalEdges`
     * 
     * @param node current node
     * @param index smallest unused index
     * @param parentNode parent node of current node
     * @param criticalEdges the critical edges will be stored here
     * @param depth depth of nodes (distance from root)
     * @param low minimum level a node can reach (without going back through parent nodes)
     * @param isVisited if a node is visited or not
     */
    void findCriticalEdges(int node, int &index, int parentNode, vector<vector<int> > &criticalEdges,
                           int depth[], int low[], bool isVisited[]);

    /**
     * @brief Get the nodes in topological order in `topologicalOrder`
     * 
     * @param node current node in DFS
     * @param topologicalOrder the nodes in topological order will be stored here
     * @param isVisited if a node is visited or not
     */
    void findTopologicalOrder(int node, vector<int> &topologicalOrder, bool isVisited[]);

    /**
     * @brief Get the Root of a node and update the root of all nodes we go through
     * 
     * @param node 
     * @param root
     * @return int Root Node
     */
    int getRootUpdatePath(int node, int root[]);

    /**
     * @brief Unite the sets by setting the root of one to the other root
     * 
     * @param node1Root 
     * @param node2Root 
     * @param root 
     * @param height 
     */
    void uniteSets(int node1Root, int node2Root, int root[], int height[]);

public:
    /**
     * @brief Construct a new Graph object
     * 
     * @param numberOfNodes 
     * @param isOriented 
     */
    Graph(int numberOfNodes, bool isOriented)
    {
        this->numberOfNodes = numberOfNodes;
        this->isOriented = isOriented;

        // Create vectors for every node
        for (int i = 0; i < numberOfNodes; i++)
        {
            vector<int> targetNodes;
            edges.push_back(targetNodes);
        }
    }

    /**
     * @brief Set the edges
     * 
     * @param connections 
     */
    void setEdges(vector<vector<int> > connections);

    /**
     * @brief Add an edge
     * 
     * @param node 
     * @param targetNode 
     */
    void addEdge(int node, int targetNode);

    /**
     * @brief Read the edges from a stream
     * 
     * @param in 
     * @param numberOfEdges 
     * @param isZeroBased 
     */
    virtual void readEdges(istream &in, int numberOfEdges, bool isZeroBased);

    /**
     * @brief Print the edges to a output stream
     */
    void printEdges(ostream &out, bool isZeroBased);

    /**
     * @brief Get the minimum distances from startNode to all nodes
     * (BFS implementation)
     * 
     * @param startNode base node from which the distances are calculated
     */
    virtual vector<int> getMinimumDistances(int startNode);

    /**
     * @brief Get the number of conex components of the Graph
     * (DFS for each unvisited node)
     * 
     * @return int number of conex components
     */
    int getNumberOfConexComponents();

    /**
     * @brief Get the biconnected components of the Graph
     * 
     * @param startNode node from which to start looking
     * @return vector<vector<int> > vector of biconnected components (list of nodes)
     */
    vector<vector<int> > getBiconnectedComponents(int startNode);

    /**
     * @brief Get the strongly connected components of the Graph
     * 
     * @return vector<vector<int> > vector of strongly connected components (list of nodes)
     */
    vector<vector<int> > getStronglyConnectedComponents();

    /**
     * @brief Get the critical edges of the Graph
     * 
     * @return vector<vector<int> > vector of critical edges (pair of nodes)
     */
    vector<vector<int> > getCriticalEdges();

    /**
     * @brief Find if the node degrees can form a graph 
     * (Havel–Hakimi algorithm)
     * 
     * @param nodeDegrees vector of node degrees
     * @return whether the degrees can form a graph
     */
    static bool isGraph(vector<int> nodeDegrees);

    /**
     * @brief Get the nodes of the Graph in topological order
     * (DFS that adds nodes to vector when the recursive call is finished)
     * 
     * @return vector<int> nodes in topological order
     */
    vector<int> getNodesInTopologicalOrder();

    /**
     * @brief Solve the tasks and return answers to queries
     * Tasks can be:
     * - UNITE SETS (By linking the root of a node to the root of the other node)
     * - QUERY SAME SET (If nodes are in the same set or not - adds the answer to answers vector)
     * 
     * @param tasks 
     * @return vector<string> answers to queries ("DA" / "NU")
     */
    vector<string> solveDisjointSetsTasks(vector<pair<int, pair<int, int> > > tasks);

    /**
     * @brief Get the Diameter of this Tree 
     * (presumes that this graph is a tree)
     * 
     * @param rootNode root of the tree (zero based!!!)
     * @return int 
     */
    int getTreeDiameter(int rootNode);
};

#pragma region GraphClassImplementation
void Graph::printEdges(ostream &out, bool isZeroBased)
{
    for (int node = 0; node < numberOfNodes; node++)
    {
        for (int i = 0; i < edges[node].size(); i++)
        {
            int baseNode = node;
            int targetNode = edges[node][i];
            if (!isZeroBased)
            {
                baseNode++;
                targetNode++;
            }
            out << baseNode << " " << targetNode << "\n";
        }
    }
}

void Graph::DFS(int node, bool isVisited[])
{
    isVisited[node] = true;
    for (int i = 0; i < edges[node].size(); i++)
    {
        int targetNode = edges[node][i];
        if (!isVisited[targetNode])
        {
            DFS(targetNode, isVisited);
        }
    }
}

void Graph::findBiconnectedComponents(int node, int currentDepth, int parentNode, vector<vector<int> > &biconnectedComponents,
                                      int depth[], int low[], bool isVisited[], stack<int> &visitedNodes)
{
    isVisited[node] = true;
    depth[node] = currentDepth;
    low[node] = currentDepth;
    visitedNodes.push(node);

    for (int i = 0; i < edges[node].size(); i++)
    {
        int targetNode = edges[node][i];
        if (targetNode != parentNode)
        {
            if (isVisited[targetNode])
            {
                low[node] = min(low[node], depth[targetNode]);
            }
            else
            {
                findBiconnectedComponents(targetNode, currentDepth + 1, node, biconnectedComponents, depth, low, isVisited, visitedNodes);
                low[node] = min(low[node], low[targetNode]);

                if (low[targetNode] >= depth[node])
                {
                    vector<int> biconnectedComponent;
                    int currentNode;
                    do
                    {
                        currentNode = visitedNodes.top();
                        biconnectedComponent.push_back(currentNode);
                        visitedNodes.pop();
                    } while (currentNode != targetNode);
                    biconnectedComponent.push_back(node);
                    biconnectedComponents.push_back(biconnectedComponent);
                }
            }
        }
    }
}

void Graph::findStronglyConnectedComponents(int node, int &index, int parentNode, vector<vector<int> > &stronglyConnectedComponents,
                                            int depth[], int low[], bool isVisited[], bool inStack[], stack<int> &visitedNodes)
{
    isVisited[node] = true;
    depth[node] = index;
    low[node] = index;
    index++;
    visitedNodes.push(node);
    inStack[node] = true;
    for (int i = 0; i < edges[node].size(); i++)
    {
        int targetNode = edges[node][i];
        if (!isVisited[targetNode])
        {
            findStronglyConnectedComponents(targetNode, index, node, stronglyConnectedComponents, depth, low, isVisited, inStack, visitedNodes);
            low[node] = min(low[node], low[targetNode]);
        }
        else if (inStack[targetNode])
        {
            low[node] = min(low[node], low[targetNode]);
        }
    }

    if (low[node] == depth[node])
    {
        vector<int> stronglyConnectedComponent;
        int currentNode;
        do
        {
            currentNode = visitedNodes.top();
            stronglyConnectedComponent.push_back(currentNode);
            visitedNodes.pop();
            inStack[currentNode] = false;
        } while (currentNode != node);
        stronglyConnectedComponents.push_back(stronglyConnectedComponent);
    }
}

void Graph::findCriticalEdges(int node, int &index, int parentNode, vector<vector<int> > &criticalEdges,
                              int depth[], int low[], bool isVisited[])
{
    isVisited[node] = true;
    depth[node] = index;
    low[node] = index;
    index++;
    for (int i = 0; i < edges[node].size(); i++)
    {
        int targetNode = edges[node][i];
        if (!isVisited[targetNode])
        {
            findCriticalEdges(targetNode, index, node, criticalEdges, depth, low, isVisited);
            low[node] = min(low[node], low[targetNode]);

            // Only if there is no other way to reach targetNode from node
            // So {node, targetNode} is a critical edge
            if (low[targetNode] == depth[targetNode])
            {
                vector<int> edge;
                edge.push_back(node);
                edge.push_back(targetNode);
                criticalEdges.push_back(edge);
            }
        }
        else if (targetNode != parentNode)
        {
            low[node] = min(low[node], low[targetNode]);
        }
    }
}

void Graph::findTopologicalOrder(int node, vector<int> &topologicalOrder, bool isVisited[])
{
    isVisited[node] = true;
    for (int i = 0; i < edges[node].size(); i++)
    {
        int targetNode = edges[node][i];
        if (!isVisited[targetNode])
        {
            findTopologicalOrder(targetNode, topologicalOrder, isVisited);
        }
    }
    topologicalOrder.push_back(node);
}

void Graph::setEdges(vector<vector<int> > connections)
{
    for (int i = 0; i < connections.size(); i++)
    {
        int baseNode = connections[i][0];
        int targetNode = connections[i][1];

        // Add edges
        edges[baseNode].push_back(targetNode);
        if (!isOriented)
        {
            edges[targetNode].push_back(baseNode);
        }
    }
}

void Graph::addEdge(int node, int targetNode)
{
    edges[node].push_back(targetNode);
}

void Graph::readEdges(istream &in, int numberOfEdges, bool isZeroBased)
{
    for (int i = 0; i < numberOfEdges; i++)
    {
        int baseNode, targetNode;
        in >> baseNode >> targetNode;

        // Make nodes zero-based
        if (!isZeroBased)
        {
            baseNode--;
            targetNode--;
        }

        // Add edges
        edges[baseNode].push_back(targetNode);
        if (!isOriented)
        {
            edges[targetNode].push_back(baseNode);
        }
    }
}

vector<int> Graph::getMinimumDistances(int startNode)
{
    queue<int> bfsNodesQueue;
    bool isVisited[numberOfNodes];
    vector<int> distances(numberOfNodes);
    for (int i = 0; i < numberOfNodes; i++)
    {
        isVisited[i] = false;
        distances[i] = NO_PATH;
    }

    // Add the start node to the distances queue. The distance is 0
    bfsNodesQueue.push(startNode);
    distances[startNode] = 0;
    isVisited[startNode] = true;

    // BFS
    while (!bfsNodesQueue.empty())
    {
        int currentNode = bfsNodesQueue.front();
        int currentDistance = distances[currentNode];

        for (int i = 0; i < edges[currentNode].size(); i++)
        {
            int targetNode = edges[currentNode][i];
            if (!isVisited[targetNode])
            {
                // If node is not visited add it to the queue and update the distance
                bfsNodesQueue.push(targetNode);
                distances[targetNode] = currentDistance + 1;
                isVisited[targetNode] = true;
            }
        }
        bfsNodesQueue.pop();
    }

    return distances;
}

int Graph::getNumberOfConexComponents()
{
    bool isVisited[numberOfNodes];
    for (int i = 0; i < numberOfNodes; i++)
    {
        isVisited[i] = false;
    }
    int numberOfConexComponents = 0;

    for (int node = 0; node < numberOfNodes; node++)
    {
        if (!isVisited[node])
        {
            DFS(node, isVisited);
            numberOfConexComponents++;
        }
    }
    return numberOfConexComponents;
}

vector<vector<int> > Graph::getBiconnectedComponents(int startNode)
{
    vector<vector<int> > biconnectedComponents;
    int depth[numberOfNodes];
    int low[numberOfNodes];
    bool isVisited[numberOfNodes];
    stack<int> visitedNodes;
    for (int node = 0; node < numberOfNodes; node++)
    {
        isVisited[node] = false;
    }

    // Call recursive function with startNode as root
    findBiconnectedComponents(startNode, 0, NO_PARENT_NODE, biconnectedComponents, depth, low, isVisited, visitedNodes);

    return biconnectedComponents;
}

vector<vector<int> > Graph::getStronglyConnectedComponents()
{
    vector<vector<int> > stronglyConnectedComponents;
    int depth[numberOfNodes];
    int low[numberOfNodes];
    bool isVisited[numberOfNodes];
    bool inStack[numberOfNodes];
    stack<int> visitedNodes;
    int index = 0;
    for (int node = 0; node < numberOfNodes; node++)
    {
        isVisited[node] = false;
        inStack[node] = false;
    }

    for (int node = 0; node < numberOfNodes; node++)
    {
        if (!isVisited[node])
        {
            // Call recursive function with node as root
            findStronglyConnectedComponents(node, index, NO_PARENT_NODE, stronglyConnectedComponents, depth, low, isVisited, inStack, visitedNodes);
        }
    }

    return stronglyConnectedComponents;
}

vector<vector<int> > Graph::getCriticalEdges()
{
    vector<vector<int> > criticalEdges;
    int depth[numberOfNodes];
    int low[numberOfNodes];
    bool isVisited[numberOfNodes];
    int index = 0;
    for (int node = 0; node < numberOfNodes; node++)
    {
        isVisited[node] = false;
    }

    for (int node = 0; node < numberOfNodes; node++)
    {
        if (!isVisited[node])
        {
            // Call recursive function with node as root
            findCriticalEdges(node, index, NO_PARENT_NODE, criticalEdges, depth, low, isVisited);
        }
    }

    return criticalEdges;
}

bool Graph::isGraph(vector<int> nodeDegrees)
{
    while (!nodeDegrees.empty())
    {
        // Sort degrees in descending order
        sort(nodeDegrees.begin(), nodeDegrees.end(), std::greater<int>());

        // Erase nodes with degree 0
        while (!nodeDegrees.empty() && nodeDegrees[nodeDegrees.size() - 1] == 0)
        {
            nodeDegrees.pop_back();
        }
        if (nodeDegrees.empty())
        {
            // All degrees are 0 so it is a graph
            return true;
        }

        // Take the highest remaining degree
        int edges = nodeDegrees[0];
        if (edges > nodeDegrees.size() - 1)
        {
            // More edges then remaining nodes
            return false;
        }

        // Consume the edges
        nodeDegrees[0] = 0;
        for (int targetNode = 1; targetNode <= edges; targetNode++)
        {
            nodeDegrees[targetNode]--;
        }
    }
    return true;
}

vector<int> Graph::getNodesInTopologicalOrder()
{
    vector<int> topologicalOrder;
    bool isVisited[numberOfNodes];
    for (int node = 0; node < numberOfNodes; node++)
    {
        isVisited[node] = false;
    }

    for (int node = 0; node < numberOfNodes; node++)
    {
        if (!isVisited[node])
        {
            // Call recursive function with node as root
            findTopologicalOrder(node, topologicalOrder, isVisited);
        }
    }

    return topologicalOrder;
}

int Graph::getRootUpdatePath(int node, int root[])
{
    // Find root node
    if (root[node] <= NO_PARENT_NODE)
    {
        return node;
    }

    root[node] = getRootUpdatePath(root[node], root);
    return root[node];
}

void Graph::uniteSets(int node1Root, int node2Root, int root[], int height[])
{
    if (node1Root != node2Root)
    {
        if (height[node1Root] > height[node2Root])
        {
            height[node1Root] += height[node2Root];
            root[node2Root] = node1Root;
        }
        else
        {
            height[node2Root] += height[node1Root];
            root[node1Root] = node2Root;
        }
    }
}

vector<string> Graph::solveDisjointSetsTasks(vector<pair<int, pair<int, int> > > tasks)
{
    int root[numberOfNodes], height[numberOfNodes];
    for (int node = 0; node < numberOfNodes; node++)
    {
        root[node] = NO_PARENT_NODE;
        height[node] = 1;
    }

    vector<string> answers;
    for (int i = 0; i < tasks.size(); i++)
    {
        int task_number = tasks[i].first;
        int node1 = tasks[i].second.first;
        int node2 = tasks[i].second.second;

        int node1Root = getRootUpdatePath(node1, root);
        int node2Root = getRootUpdatePath(node2, root);

        if (task_number == TASK_NUMBER_UNITE_SETS)
        {
            uniteSets(node1Root, node2Root, root, height);
        }

        if (task_number == TASK_NUMBER_QUERY_SAME_SET)
        {
            if (node1Root == node2Root)
            {
                answers.push_back("DA");
            }
            else
            {
                answers.push_back("NU");
            }
        }
    }
    return answers;
}

int Graph::getTreeDiameter(int rootNode)
{
    // Find furthest node from root
    vector<int> minimumDistances = getMinimumDistances(rootNode);
    int furthestNode = rootNode;
    int maximumDistance = 0;
    for (int i = 0; i < minimumDistances.size(); i++)
    {
        if (minimumDistances[i] > maximumDistance)
        {
            maximumDistance = minimumDistances[i];
            furthestNode = i;
        }
    }

    // Find furthest node from the previous furthest node
    minimumDistances = getMinimumDistances(furthestNode);
    maximumDistance = 0;
    for (int i = 0; i < minimumDistances.size(); i++)
    {
        if (minimumDistances[i] > maximumDistance)
        {
            maximumDistance = minimumDistances[i];
        }
    }

    return maximumDistance + 1;
}

#pragma endregion EndGrapClassImplementation

class Solution
{
public:
    vector<vector<int> > criticalConnections(int numberOfNodes, vector<vector<int> > &connections)
    {
        Graph graph(numberOfNodes, false);
        graph.setEdges(connections);
        return graph.getCriticalEdges();
    }
};

class WeightedGraph : public Graph
{
private:
    // Maps the edges to their weights
    map<pair<int, int>, int> weightMap;

    /**
     * @brief Get a vector of all edges sorted by weight
     * 
     * @return vector<pair<int, pair<int, int> > > sorted edges {cost, {node, targetNode}}
     */
    vector<pair<int, pair<int, int> > > getSortedEdges();

public:
    /**
     * @brief Construct a new Weighted Graph object
     * 
     * @param numberOfNodes 
     * @param isOriented 
     */
    WeightedGraph(int numberOfNodes, bool isOriented)
        : Graph(numberOfNodes, isOriented) {}

    /**
     * @brief Read edges from a input stream
     * 
     * @param in 
     * @param numberOfEdges 
     * @param isZeroBased 
     */
    void readEdges(istream &in, int numberOfEdges, bool isZeroBased);

    /**
     * @brief Read the edges from a stream (matrix format)
     * 
     * @param in 
     */
    void readEdgesFromMatrix(istream &in);

    /**
     * @brief Get the Adjacency Matrix of the Graph
     * 
     * @return vector<vector<int> > 
     */
    vector<vector<int> > getAdjacencyMatrix();

    /**
     * @brief Prints a Matrix in a stream
     * 
     * @param out 
     * @param matrix
     */
    static void printMatrix(ostream &out, vector<vector<int> > matrix);

    /**
     * @brief Get the minimum distances from startNode to all nodes. 
     * (Dijkstra Algorithm)
     * WARNING: Does not work for negative weights!
     * 
     * @param startNode base node from which the distances are calculated
     * @return the distances from the startNode to all other nodes
     */
    vector<int> getMinimumDistances(int startNode);

    /**
     * @brief Get the minimum distances from startNode to all nodes. 
     * (Bellman Ford Algorithm)
     * Throws error when there is a negative cycle
     * 
     * @param startNode base node from which the distances are calculated
     * @return the distances from the startNode to all other nodes
     */
    vector<int> getMinimumDistancesNegativeWeights(int startNode);

    // Faster function
    vector<int> getMinimumDistancesNegativeWeights(int startNode, vector<vector<pair<int, int> > > weightedEdges);

    /**
     * @brief Get the Minimum Spanning Tree of the Graph.
     * (Kruskal Algorithm: 
     * Sort all edges
     * Iterate through the edges. If the 2 nodes are not in the same component add the edge
     * Stop when we have N-1 edges)
     * 
     * @param totalCost The cost will be stored here
     * @return Graph 
     */
    Graph getMinimumSpanningTree(int &totalCost);
};

#pragma region WeightedGraphClassImplementation
void WeightedGraph::readEdges(istream &in, int numberOfEdges, bool isZeroBased)
{

    for (int i = 0; i < numberOfEdges; i++)
    {
        int baseNode, targetNode, weight;
        in >> baseNode >> targetNode >> weight;

        // Make nodes zero-based
        if (!isZeroBased)
        {
            baseNode--;
            targetNode--;
        }

        // Add edges
        edges[baseNode].push_back(targetNode);
        if (weightMap.find(make_pair(baseNode, targetNode)) == weightMap.end())
        {
            weightMap[make_pair(baseNode, targetNode)] = weight;
        }
        else
        {
            int minWeight = min(weight, weightMap[make_pair(baseNode, targetNode)]);
            weightMap[make_pair(baseNode, targetNode)] = minWeight;
        }
    }
}

void WeightedGraph::readEdgesFromMatrix(istream &in)
{
    for (int i = 0; i < numberOfNodes; i++)
    {
        for (int j = 0; j < numberOfNodes; j++)
        {
            int weight;
            fin >> weight;
            edges[i].push_back(j);
            weightMap[make_pair(i, j)] = weight;
        }
    }
}

vector<vector<int> > WeightedGraph::getAdjacencyMatrix()
{
    vector<vector<int> > adjacencyMatrix(numberOfNodes);
    for (int i = 0; i < numberOfNodes; i++)
    {
        adjacencyMatrix[i] = vector<int>(numberOfNodes);
        for (int j = 0; j < numberOfNodes; j++)
        {
            if (weightMap.find(make_pair(i, j)) == weightMap.end())
            {
                adjacencyMatrix[i][j] = NO_EDGE;
            }
            else
            {
                adjacencyMatrix[i][j] = weightMap[make_pair(i, j)];
            }
        }
    }
    return adjacencyMatrix;
}

void WeightedGraph::printMatrix(ostream &out, vector<vector<int> > matrix)
{
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++)
        {
            out << matrix[i][j] << " ";
        }
        out << "\n";
    }
}

vector<pair<int, pair<int, int> > > WeightedGraph::getSortedEdges()
{
    vector<pair<int, pair<int, int> > > sortedEdges;
    for (int node = 0; node < numberOfNodes; node++)
    {
        for (int i = 0; i < edges[node].size(); i++)
        {
            int targetNode = edges[node][i];
            if (weightMap.find(make_pair(node, targetNode)) == weightMap.end())
            {
                throw "No weight found for this edge!";
            }
            int cost = weightMap[make_pair(node, targetNode)];
            sortedEdges.push_back(make_pair(cost, make_pair(node, targetNode)));
        }
    }
    sort(sortedEdges.begin(), sortedEdges.end());
    return sortedEdges;
}

Graph WeightedGraph::getMinimumSpanningTree(int &totalCost)
{
    Graph minimumSpanningTree(numberOfNodes, true);

    vector<pair<int, pair<int, int> > > sortedEdges = getSortedEdges();

    int root[numberOfNodes];
    for (int node = 0; node < numberOfNodes; node++)
    {
        root[node] = NO_PARENT_NODE;
    }

    totalCost = 0;
    int numberOfEdgesInTree = 0;
    for (int i = 0; i < sortedEdges.size(); i++)
    {
        // Check if all nodes are in the tree
        if (numberOfEdgesInTree == numberOfNodes - 1)
        {
            break;
        }

        int cost = sortedEdges[i].first;
        int node = sortedEdges[i].second.first;
        int targetNode = sortedEdges[i].second.second;

        int nodeRoot = getRootUpdatePath(node, root);
        int targetNodeRoot = getRootUpdatePath(targetNode, root);

        if (nodeRoot != targetNodeRoot)
        {
            // Add to solution
            numberOfEdgesInTree++;
            totalCost += cost;
            minimumSpanningTree.addEdge(node, targetNode);

            // Make both trees have the same root
            root[nodeRoot] = targetNodeRoot;
        }
    }
    return minimumSpanningTree;
}

vector<int> WeightedGraph::getMinimumDistancesNegativeWeights(int startNode, vector<vector<pair<int, int> > > weightedEdges)
{
    vector<int> minimumDistance(numberOfNodes);
    int frequency[numberOfNodes];
    bool inQueue[numberOfNodes];
    queue<int> nodesQueue;
    for (int i = 0; i < numberOfNodes; i++)
    {
        minimumDistance[i] = MAX_DISTANCE;
        frequency[i] = 0;
        inQueue[i] = false;
    }

    nodesQueue.push(startNode);
    minimumDistance[startNode] = 0;
    inQueue[startNode] = true;

    while (!nodesQueue.empty())
    {
        int node = nodesQueue.front();
        frequency[node]++;
        if (frequency[node] > numberOfNodes)
        {
            string error("Ciclu negativ!");
            throw error;
        }
        for (int i = 0; i < weightedEdges[node].size(); i++)
        {
            int targetNode = weightedEdges[node][i].first;
            int cost = weightedEdges[node][i].second;
            if (minimumDistance[node] + cost < minimumDistance[targetNode])
            {
                minimumDistance[targetNode] = minimumDistance[node] + cost;
                if (!inQueue[targetNode])
                {
                    nodesQueue.push(targetNode);
                    inQueue[targetNode] = true;
                }
            }
        }
        nodesQueue.pop();
        inQueue[node] = false;
    }
    return minimumDistance;
}

vector<int> WeightedGraph::getMinimumDistancesNegativeWeights(int startNode)
{
    vector<int> minimumDistance(numberOfNodes);
    int frequency[numberOfNodes];
    bool inQueue[numberOfNodes];
    queue<int> nodesQueue;
    for (int i = 0; i < numberOfNodes; i++)
    {
        minimumDistance[i] = MAX_DISTANCE;
        frequency[i] = 0;
        inQueue[i] = false;
    }

    nodesQueue.push(startNode);
    minimumDistance[startNode] = 0;
    inQueue[startNode] = true;

    while (!nodesQueue.empty())
    {
        int node = nodesQueue.front();
        frequency[node]++;
        if (frequency[node] > numberOfNodes)
        {
            string error("Ciclu negativ!");
            throw error;
        }
        for (int i = 0; i < edges[node].size(); i++)
        {
            int targetNode = edges[node][i];
            int cost = weightMap[make_pair(node, targetNode)];
            if (minimumDistance[node] + cost < minimumDistance[targetNode])
            {
                minimumDistance[targetNode] = minimumDistance[node] + cost;
                if (!inQueue[targetNode])
                {
                    nodesQueue.push(targetNode);
                    inQueue[targetNode] = true;
                }
            }
        }
        nodesQueue.pop();
        inQueue[node] = false;
    }
    return minimumDistance;
}

vector<int> WeightedGraph::getMinimumDistances(int startNode)
{
    vector<int> minimumDistance(numberOfNodes);
    set<pair<int, int> > nodeDistanceSet;
    for (int i = 0; i < numberOfNodes; i++)
    {
        minimumDistance[i] = MAX_DISTANCE;
    }

    nodeDistanceSet.insert(make_pair(0, startNode));
    minimumDistance[startNode] = 0;

    while (!nodeDistanceSet.empty())
    {
        int currentNode = nodeDistanceSet.begin()->second;
        int currentCost = minimumDistance[currentNode];
        nodeDistanceSet.erase(nodeDistanceSet.begin());

        for (int i = 0; i < edges[currentNode].size(); i++)
        {
            int targetNode = edges[currentNode][i];
            int targetCost = weightMap[make_pair(currentNode, targetNode)];
            if (currentCost + targetCost < minimumDistance[targetNode])
            {
                nodeDistanceSet.erase(make_pair(currentNode, currentCost));
                minimumDistance[targetNode] = currentCost + targetCost;
                nodeDistanceSet.insert(make_pair(currentCost + targetCost, targetNode));
            }
        }
    }
    for (int node = 0; node < numberOfNodes; node++)
    {
        if (minimumDistance[node] == MAX_DISTANCE)
        {
            minimumDistance[node] = 0;
        }
    }
    return minimumDistance;
}
#pragma endregion EndWeightedGraphClassImplementation

int main()
{
    int numberOfNodes;
    fin >> numberOfNodes;

    WeightedGraph graph(numberOfNodes, false);
    graph.readEdgesFromMatrix(fin);

    graph.printMatrix(fout, graph.getAdjacencyMatrix());
}