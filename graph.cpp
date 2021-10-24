#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <iostream>
using namespace std;
ifstream fin("ctc.in");
ofstream fout("ctc.out");

#define NO_PATH -1
#define NO_PARENT_NODE -1

/**
 * @brief Graph class that implements algorithms
 * 
 */
class Graph
{
private:
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
    void DFS(int node, bool isVisited[])
    {
        isVisited[node] = 1;
        for (int i = 0; i < edges[node].size(); i++)
        {
            int targetNode = edges[node][i];
            if (!isVisited[targetNode])
            {
                DFS(targetNode, isVisited);
            }
        }
    }

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
     * @brief Read the edges from a stream
     * 
     * @param in 
     * @param numberOfEdges 
     * @param isZeroBased 
     */
    void readEdges(istream &in, int numberOfEdges, bool isZeroBased)
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

    /**
     * @brief Get the minimum distances from startNode to all nodes
     * 
     * @param startNode base node from which the distances are calculated
     */
    vector<int> getMinimumDistances(int startNode)
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

    /**
     * @brief Get the number of conex components of the Graph
     * 
     * @return int number of conex components
     */
    int getNumberOfConexComponents()
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

    /**
     * @brief Get the biconnected components of the Graph
     * 
     * @param startNode node from which to start looking
     * @return vector<vector<int> > vector of biconnected components (list of nodes)
     */
    vector<vector<int> > getBiconnectedComponents(int startNode)
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

    /**
     * @brief Get the strongly connected components of the Graph
     * 
     * @return vector<vector<int> > vector of strongly connected components (list of nodes)
     */
    vector<vector<int> > getStronglyConnectedComponents()
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
};

int main()
{
    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph graph(numberOfNodes, true);
    graph.readEdges(fin, numberOfEdges, false);

    vector<vector<int> > stronglyConnectedComponents = graph.getStronglyConnectedComponents();

    fout << stronglyConnectedComponents.size() << "\n";
    for (int i = 0; i < stronglyConnectedComponents.size(); i++)
    {
        vector<int> component = stronglyConnectedComponents[i];
        for (int j = 0; j < component.size(); j++)
        {
            int node = component[j] + 1;
            fout << node << " ";
        }
        fout << "\n";
    }
}