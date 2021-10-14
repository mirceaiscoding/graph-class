#include <vector>
#include <queue>
#include <fstream>
#include <iostream>
using namespace std;
ifstream fin("dfs.in");
ofstream fout("dfs.out");

#define NO_PATH -1

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
     * @brief Get the minimum distances
     * 
     * @param startNode 
     * @param distances array where the distances will be placed. It should have [number_of_nodes] elements
     */
    void getMinimumDistances(int startNode, int distances[])
    {
        queue<int> bfsNodesQueue;
        bool isVisited[numberOfNodes];
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
};

int main()
{
    int numberOfNodes, numberOfEdges;
    fin >> numberOfNodes >> numberOfEdges;

    Graph graph(numberOfNodes, false);
    graph.readEdges(fin, numberOfEdges, false);

    fout << graph.getNumberOfConexComponents();
}