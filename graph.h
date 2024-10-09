/*
graph.h

Visible structs and functions for graph construction and manipulation.

*/

/* Definition of a graph. */
struct graph;

enum problemPart;

struct solution;

/* Moved class decleration to this header file so I can pass
   struct edge array to function */
struct edge;

/* A particular solution to a graph problem. */
#ifndef SOLUTION_STRUCT
#define SOLUTION_STRUCT
struct solution {
  int connectedSubnets;
  int largestSubnet;
  int *largestSubnetSIDs;
  int postOutageDiameter;
  int postOutageDiameterCount;
  int *postOutageDiameterSIDs;
  int criticalServerCount;
  int *criticalServerSIDs;
};
#endif

/* Which part the program should find a solution for. */
#ifndef PART_ENUM
#define PART_ENUM
enum problemPart {
  TASK_2=0,
  TASK_3=1,
  TASK_4=2,
  TASK_7=3
};
#endif

/* Creates an undirected graph with the given numVertices and no edges and
returns a pointer to it. NumEdges is the number of expected edges. */
struct graph *newGraph(int numVertices);

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end);

/* Finds:
  - Number of connected subnetworks (before outage) (Task 2)
  - Number of servers in largest subnetwork (before outage) (Task 3)
  - SIDs of servers in largest subnetwork (before outage) (Task 3)
  - Diameter of largest subnetworks (after outage) (Task 4)
  - Number of servers in path with largest diameter - should be one more than
    Diameter if a path exists (after outage) (Task 4)
  - SIDs in largest subnetwork (after outage) (Task 4)
  - Number of critical servers (before outage) (Task 7)
  - SIDs of critical servers (before outage) (Task 7)
 */
struct solution *graphSolve(struct graph *g, enum problemPart part,
  int numServers, int numOutages, int *outages);

/* Frees all memory used by graph. */
void freeGraph(struct graph *g);

/* Sets all values to initial values so free can work for all tasks without change. */
void initaliseSolution(struct solution *solution);

/* Frees all data used by solution. */
void freeSolution(struct solution *solution);

/* Initialises visited array, sets all servers as unvisited */
void initVisited(int[], int);

/* Performs DFS search to find number of connected subnetworks */
void dfsExplore(enum problemPart,
                int, struct edge **, int, int[],
                int*);

/* Records SIDs of servers in a subnetwork */
void recordSIDs(enum problemPart, int, struct edge **, int, int[],
                int *, int *);

/* Finds all critial servers in a network */
void dfsExploreCritical(int, int, struct edge **, int, int[],
                        int[], int[], int[], int[], int);

/* Performs Dijkstra's algorithm to find shortest paths */
void dijkstra(int[], int, struct edge **, int, int[], int[], int[],
              int*, int*, int);

/* Write longest shortest path to solution array */
void retracePath(int[], int *, int[], int*, int);

/* Checks if all post-outage servers in a subnetwork have been
   visited */
int subnetAllVisited(int[], int[], int);

/* Extracts current minimum priority server */
int getMin(int[], int[], int, int[]);

/* Checks for adjacency between two nodes */
int adjacency(int, int, struct edge **, int);

/* Sorts an array */
void insertionSort(int *, int);

/* Swaps two sequential elements in an array */
void swap(int *, int);

/* Returns the maximum value between two ints */
int max(int, int);

/* Returns the minimum value between two ints */
int min(int, int);

/* Checks if an element is in an array */
int inArr(int *, int, int);

/* DEBUGGING FUNCTIONS */
void printArr(int *, int);
