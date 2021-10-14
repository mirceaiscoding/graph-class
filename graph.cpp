#include <vector>
#include <queue>
#include <fstream>
#include <iostream>
using namespace std;
ifstream fin("bfs.in");
ofstream fout("bfs.out");

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
};

int main()
{
    int numberOfNodes, numberOfEdges, startNode;
    fin >> numberOfNodes >> numberOfEdges >> startNode;
    startNode--; // make it zero-based

    Graph graph(numberOfNodes, true);
    graph.readEdges(fin, numberOfEdges, false);

    int distances[numberOfNodes];
    graph.getMinimumDistances(startNode, distances);
    for (int i = 0; i < numberOfNodes; i++)
    {
        fout << distances[i] << " ";
    }
}