#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <random>
#include <queue>
#include <vector>
#include <unordered_map>
#include <chrono>

// node object
struct keyValue {
int key;
double value;

keyValue(int keyinput, double valueinput) {
    key = keyinput;
    value = valueinput;
}
};

// make a heap with a vector of key-value pairs.
class BinaryHeap {
public:
    std::vector<keyValue> nodes;
    // this map is to make change O(log n).
    // it creates a mapping key -> index of the key-value pair in nodes
    std::unordered_map<int, int> keyToIndexMap;

    BinaryHeap() {
        // empty constructor
        // heap constructed as empty. elements inserted later.
    }
    
    bool empty() {
        return nodes.empty();
    }

    // Insert a new element with double value
    void insert(int key, double value) {
        nodes.push_back(keyValue(key, value));
        keyToIndexMap[key] = nodes.size()-1;
        int curindex = nodes.size() - 1;        
        bubbleup(curindex);
    }

    // note that this should be O(log n).
    // finding the key uses the hashing in unordered_map so O(1)
    // bubbling up or down is O(log n)
    void change(int key, double newValue) {
        if (keyToIndexMap.find(key) == keyToIndexMap.end()) {
            // if key is not found, simply insert
            insert(key, newValue);
            return;
        }
    
        int curindex = keyToIndexMap[key];
        int previndex = (curindex - 1) / 2;
        nodes[curindex].value = newValue;
    
        if (previndex < 0) {
            return;
        }
        bubbleup(curindex);
        bubbledown(curindex);
    }

    int deletemin() {
        if (nodes.empty()) {
            printf("Heap is empty. Cannot delete minimum.");
            return -1; // indicating an error or an empty heap.
        }
        int returnkey = nodes[0].key;
        keyToIndexMap.erase(returnkey);
        nodes[0] = nodes.back(); // brings the last element forward
        nodes.pop_back(); // deletes the last element.
        bubbledown(0);
        return returnkey;
    }

    void bubbleup(int curindex) {
        int previndex = (curindex - 1) / 2;
        if (previndex < 0) {
            return;
        }
        while (nodes[curindex].value < nodes[previndex].value) {
            if (previndex < 0) {
                break;
            }
            // the swap happens here.
            std::swap(nodes[curindex], nodes[previndex]);
            // makes sure the keyToIndexMap is well preserved.
            keyToIndexMap[nodes[previndex].key] = previndex;
            keyToIndexMap[nodes[curindex].key] = curindex;
        
            curindex = previndex;
            previndex = (curindex - 1) / 2;
        }
    }

    void bubbledown(int curindex) {
        while (true) {
            int nextindex = curindex * 2 + 1; // left child index.
            int rightindex = curindex * 2 + 2; // right child index.
            // find the smaller of the two children.
            if (rightindex < nodes.size() && nodes[rightindex].value < nodes[nextindex].value) {
                nextindex = rightindex;
            }
            // if no children or the heap property is not violated, break.
            if (nextindex >= nodes.size() || nodes[curindex].value <= nodes[nextindex].value) {
                break;
            }
            // swap with the smaller child
            std::swap(nodes[nextindex], nodes[curindex]);
            keyToIndexMap[nodes[curindex].key] = curindex;
            keyToIndexMap[nodes[nextindex].key] = nextindex;
            // move down the heap
            curindex = nextindex;
        }
    }
};

struct Edge {
    int target;
    double weight;


    Edge(int targ, double wgt) {
        target = targ;
        weight = wgt;
    }
};

class Graph {
    public:
    // a graph has a collection of vertices
    int vertices;
    std::vector<Edge>* edges; // vector of edges

    Graph(int num) {
        // construct a graph with n vertices and no edges
        vertices = num;
        edges = new std::vector<Edge>[num];
    }

    void addEdge(int src, int target, double weight) {
        edges[src].push_back(Edge(target, weight));
    }
};

// The maxedge cutoff was tested for, as reported in our writeup.
// We return the cutoffs here given n, dimension.
double maxedge (int n, int dimension) {
    if (dimension == 0) {
        return 16 / (float) n;
    }
    else if (dimension == 2) {
        return 3 / sqrt(n);
    }
    else if (dimension == 3) {
        return 2.2 / std::pow(n, 1/3.);
    }
    else if (dimension == 4) {
        return 2.2 / std::pow(n, 1/4.);
    }
    return -1;
}

Graph zerodimGraph(int n, int dimension) {
    double cutoff = maxedge(n, dimension);
    Graph G = Graph(n);
    // Set up random number generation
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double temp = dist(e2);
            if(temp < cutoff){
                G.addEdge(i,j,temp);
                G.addEdge(j,i,temp);
            }
        }
    }
    return G;
}

Graph euclidGraph(int n, int dimension) {
    double cutoff = maxedge(n, dimension);
    // printf("%f\n", cutoff);
    Graph G = Graph(n);

    // Set up random number generation
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    // Fill coordinate array with random coordinates
    double** coordinates = new double*[n];

    for (int i = 0; i < n; ++i) {
        coordinates[i] = new double[dimension];  // Allocate memory for each row
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dimension; j++) {
            coordinates[i][j] = dist(e2);
        }
    }

    // Fill graph with euclidean distances between the coordinates
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double distance = 0;
            for (int k = 0; k < dimension; k++) {
                distance += pow((coordinates[i][k] - coordinates[j][k]),2);
            }
            distance = sqrt(distance);
            if(distance < cutoff){
                G.addEdge(i,j,distance);
                G.addEdge(j,i,distance);
            }
        }
    }

    // Delete coordinate array
    for (int i = 0; i < n; ++i) {
        delete[] coordinates[i];
    }
    delete[] coordinates;

    return G;
}

// Graph construction
Graph constructor (int n, int dimension) {
    if (dimension == 0) {
        Graph G = zerodimGraph(n, dimension);
        return G;
    }
    else {
        Graph G = euclidGraph(n, dimension);
        return G;
    }
}

// Prim's algorithm
std::vector<double> prim(Graph G, int n) {
    // Prim's algorithm -- implicitly, vertices are labeled 0 through n-1
    int s = 0; // starting key value
    std::vector<double> d(n, 2); // distances
    int Prev[n]; // parents
    BinaryHeap H; // heap
    int inHeap[n]; // keeps track of heap contents -- gives 1 in index i if vertex i in heap, otherwise 0
    int inS[n]; // keeps track of S contents

    // Fill these with values - we take s=0 here
    d[s] = 0;
    Prev[0] = -1;
    for (int i = 1; i < n; i++) {
        d[i] = 2;
        inHeap[i] = 0;
        inS[i] = 0;
    }

    H.insert(s,0);

    // Run the for loop portion of Prim
    while (!H.empty()) {
        int u = H.deletemin();
        inS[u] = 1;
        inHeap[u] = 0;

        // we examine all edges starting from u
        for (int i = 0; i < G.edges[u].size(); i++){
            int v = G.edges[u][i].target;
            if (inS[v] == 0 && d[v] > G.edges[u][i].weight) {
                d[v] = G.edges[u][i].weight;
                Prev[v] = u;
                // and we now update the heap; this could be a change or an insert
                    if (inHeap[v] == 0) {
                        H.insert(v, d[v]);
                        inHeap[v] = 1; // v is now in the heap
                    } else {
                        H.change(v, d[v]); // otherwise we simply update
                    }
            }
        }
    }
    return d;
}

double max(std::vector<double> distances) {
    double max_edge = 0;
    for (int i = 0; i < distances.size(); i++){
        if (distances[i] > max_edge) {
            max_edge = distances[i];
        }
    }
    return max_edge;
}

double min_tree_length (std::vector<double> distances) {
    // calculate length of the minimum tree
    double min_tree_length = 0;
    for (int i = 0; i < distances.size(); i++){
        min_tree_length += distances[i]; // this should give tree length; double check
    }
    return min_tree_length;
}

double average (std::vector<double> values) {
    double sum = 0;
    for (int i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    return sum / (float) values.size();
}

// This collects data for some n.
// In particular it stores the tree_lengths for individual trials.
struct dataCollector {
    int n;
    std::vector<double> trials;
    double average;
};

// makes the desired MSTs for a certain number of trials
dataCollector makeTrees (int n, int numtrials, int dimension){
    dataCollector myTrees;
    myTrees.n = n;
    for (int i = 0; i < numtrials; i++) {
        Graph G = constructor(n, dimension);
        std::vector<double> distances = prim(G, n);
        myTrees.trials.push_back(min_tree_length(distances));
    }
    myTrees.average = average(myTrees.trials);
    return myTrees;
}

void outputs(std::vector<int> n, int numtrials, int dimension){
    for (int i = 0; i < n.size(); i++) {
        dataCollector myTrees = makeTrees(n[i], numtrials, dimension);
        // In case you want individual trials
        // for (int i = 0; i < myTrees.trials.size(); i++) {
        //     printf("%f\n", myTrees.trials[i]);
        // }
        printf("%f %i %i %i \n", myTrees.average, n[i], numtrials, dimension);
    }
}


int main(int argc, char *argv[]) {
    // Initialize empty graph
    if (argc == 5) {
        int flag = std::stoi(argv[1]);
        int numpoints = std::stoi(argv[2]);
        int numtrials = std::stoi(argv[3]);
        int dimension = std::stoi(argv[4]);
        std::vector<int> n; // use a vector here so we could test multiple n
        n.push_back(numpoints);
        outputs(n, numtrials, dimension);
    }
    else {
        printf("Wrong input. Please try again. ");
    }
}

