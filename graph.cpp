#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;
ifstream fin("sortaret.in");
ofstream fout("sortaret.out");

#define NO_PATH -1
#define NO_PARENT_NODE -1

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
    }

    /**
     * @brief Set the edges
     * 
     * @param connections 
     */
    void setEdges(vector<vector<int> > connections);

    /**
     * @brief Read the edges from a stream
     * 
     * @param in 
     * @param numberOfEdges 
     * @param isZeroBased 
     */
    void readEdges(istream &in, int numberOfEdges, bool isZeroBased);

    /**
     * @brief Get the minimum distances from startNode to all nodes
     * 
     * @param startNode base node from which the distances are calculated
     */
    vector<int> getMinimumDistances(int startNode);

    /**
     * @brief Get the number of conex components of the Graph
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
     * 
     * @return vector<int> nodes in topological order
     */
    vector<int> getNodesInTopologicalOrder();
};

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
    // Create vectors for every node
    for (int i = 0; i < numberOfNodes; i++)
    {
        vector<int> targetNodes;
        edges.push_back(targetNodes);
    }

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

void Graph::readEdges(istream &in, int numberOfEdges, bool isZeroBased)
{
    // Create vectors for every node
    for (int i = 0; i < numberOfNodes; i++)
    {
        vector<int> targetNodes;
        edges.push_back(targetNodes);
    }

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

class WeightedGraph : Graph
{
private:
    // Maps the edges to their weights
    map<pair<int, int>, int> weightMap;

    // Sorted edges
    vector<pair<int, pair<int, int> > > sortedEdges;

    // Compare two edges based on weight
    bool compareEdges(pair<int, int> edge1, pair<int, int>edge2)
    {
        return weightMap[edge1] < weightMap[edge2];
    }

public:
    /**
     * @brief Construct a new weighted Graph object
     * 
     * @param numberOfNodes 
     * @param isOriented 
     */
    WeightedGraph(int numberOfNodes, bool isOriented) 
    : Graph(numberOfNodes, isOriented){}

    /**
     * @brief Read the edges from a stream
     * 
     * @param in 
     * @param numberOfEdges 
     * @param isZeroBased 
     * @param sortedEdges
     */
    void readEdges(istream &in, int numberOfEdges, bool isZeroBased, bool holdSortedEdges);

};

void WeightedGraph::readEdges(istream &in, int numberOfEdges, bool isZeroBased, bool holdSortedEdges)
{
    // Create vectors for every node
    for (int i = 0; i < numberOfNodes; i++)
    {
        vector<int> targetNodes;
        edges.push_back(targetNodes);
    }

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
        weightMap.insert(make_pair(make_pair(baseNode, targetNode), weight));
        sortedEdges.push_back(make_pair(baseNode, targetNode));
        if (!isOriented)
        {
            edges[targetNode].push_back(baseNode);
            weightMap.insert(make_pair(make_pair(targetNode, baseNode), weight));
            sortedEdges.push_back(make_pair(targetNode, baseNode));
        }
    }
    sort(sortedEdges.begin(), sortedEdges.end(), compareEdges);
}


int main()
{
    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    WeightedGraph graph(numberOfNodes, true);
    graph.readEdges(fin, numberOfEdges, false, true);

}