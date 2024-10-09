/*
graph.c

Set of vertices and edges implementation.

Implementations for helper functions for graph construction and manipulation.

*/
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include "graph.h"
#include "utils.h"
#include "pq.h"

#define INITIALEDGES 32
#define UNVISITED 0
#define VISITED 1
#define TRUE_INT 1
#define FALSE_INT 0
#define NULL_INT -1

//struct edge;

/* Definition of a graph. */
struct graph {
  int numVertices;
  int numEdges;
  int allocedEdges;
  struct edge **edgeList;
};

/* Definition of an edge. */
struct edge {
  int start;
  int end;
};

struct graph *newGraph(int numVertices){
  struct graph *g = (struct graph *) malloc(sizeof(struct graph));
  assert(g);
  /* Initialise edges. */
  g->numVertices = numVertices;
  g->numEdges = 0;
  g->allocedEdges = 0;
  g->edgeList = NULL;
  return g;
}

/* Adds an edge to the given graph. */
void addEdge(struct graph *g, int start, int end){
  assert(g);
  struct edge *newEdge = NULL;
  /* Check we have enough space for the new edge. */
  if((g->numEdges + 1) > g->allocedEdges){
    if(g->allocedEdges == 0){
      g->allocedEdges = INITIALEDGES;
    } else {
      (g->allocedEdges) *= 2;
    }
    g->edgeList = (struct edge **) realloc(g->edgeList,
      sizeof(struct edge *) * g->allocedEdges);
    assert(g->edgeList);
  }

  /* Create the edge */
  newEdge = (struct edge *) malloc(sizeof(struct edge));
  assert(newEdge);
  newEdge->start = start;
  newEdge->end = end;

  /* Add the edge to the list of edges. */
  g->edgeList[g->numEdges] = newEdge;
  (g->numEdges)++;
}

/* Frees all memory used by graph. */
void freeGraph(struct graph *g){
  int i;
  for(i = 0; i < g->numEdges; i++){
    free((g->edgeList)[i]);
  }
  if(g->edgeList){
    free(g->edgeList);
  }
  free(g);
}

/* Finds:
  - Number of connected subnetworks (before outage) (Task 2)
  - Number of servers in largest subnetwork (before outage) (Task 3)
  - SIDs of servers in largest subnetwork (before outage) (Task 3)
  - Diameter of largest subnetworks (after outage) (Task 4)
  - Number of servers in path with largest diameter - should be one more than
    Diameter if a path exists (after outage) (Task 4)
  - SIDs in path with largest diameter (after outage) (Task 4)
  - Number of critical servers (before outage) (Task 7)
  - SIDs of critical servers (before outage) (Task 7)
*/
struct solution *graphSolve(struct graph *g, enum problemPart part,
  int numServers, int numOutages, int *outages){

  struct solution *solution = (struct solution *)
    malloc(sizeof(struct solution));
  assert(solution);

  /* Initialise solution values */
  initaliseSolution(solution);

  if(part == TASK_2){
      
      int v; // server loop counter
      
      /* Stores number of connected subnetworks */
      int connectedSubnets = 0;
      
      /* Stores number of servers in current subnetwork (REQ.TASK3) */
      int currServerCount = 0;
      
      /* Keeps track of which servers have been visited */
      int visited[numServers];
      initVisited(visited, numServers);
      
      /* Activate DFS search on network */
      for (v=0; v<numServers; v++) {
          /* Found a new subnetwork! */
          if (!visited[v]) {
              dfsExplore(part, v, g->edgeList, g->numEdges, visited,
                         &currServerCount);
              connectedSubnets++;
          }
      }
      solution->connectedSubnets = connectedSubnets;

  } else if(part == TASK_3) {
    
      int v; // server (vertice) loop counter
      
      /* Stores number of servers in current subnetwork */
      int currServerCount = 0;
      
      /* Stores number of servers in largest subnetwork */
      int maxServerCount = 0;
      
      /* Keeps track of which servers have been visited */
      int visited[numServers];
      initVisited(visited, numServers);
      int reVisited[numServers];
      initVisited(reVisited, numServers);
      
      /* Initialize solution array */
      solution->largestSubnetSIDs =
      malloc(sizeof(solution->largestSubnetSIDs) * numServers);
      assert(solution->largestSubnetSIDs);
      
      /* Activate DFS search on network */
      for (v=0; v<numServers; v++) {
          /* Found a new subnetwork! */
          if (!visited[v]) {
              dfsExplore(part, v, g->edgeList, g->numEdges, visited,
                         &currServerCount);
              
              /* Found new largest subnet, record its SIDs */
              if (currServerCount > maxServerCount) {
                  int insertIndex = 0;
                  recordSIDs(part, v, g->edgeList, g->numEdges,
                             reVisited, solution->largestSubnetSIDs,
                             &insertIndex);
              }
              /* Check for new max server count */
              maxServerCount = max(currServerCount, maxServerCount);
              currServerCount = 0;
          }
      }
      solution->largestSubnet = maxServerCount;
      
      /* Sort the solution array */
      insertionSort(solution->largestSubnetSIDs,
                    solution->largestSubnet);
      
  } else if(part == TASK_4) {
    
      /* Initialize solution array */
      solution->postOutageDiameterSIDs =
      malloc(sizeof(solution->postOutageDiameterSIDs) * numServers);
      assert(solution->postOutageDiameterSIDs);
      
      /* Initialization */
      int newServerCount=0;
      int currDiameter=0;
      int subnetServers[numServers];
      int dfsVisited[numServers];
      int visited[numServers];
      int dist[numServers];
      int prev[numServers];
      
      int v, i;
      for (v=0; v<numServers; v++) {
          dfsVisited[v] = UNVISITED;
          visited[v] = UNVISITED;
          dist[v] = INT_MAX;
          prev[v] = NULL_INT;
      }
      for (i=0; i<numOutages; i++) {
          dfsVisited[outages[i]] = VISITED;
          visited[outages[i]] = NULL_INT;
      }
      /* Find diameter of graph */
      for (v=0; v<numServers; v++) {
          /* Found new subnetwork */
          if (!dfsVisited[v]) {
              /* Find servers in subnetwork post-outage */
              int insertIndex = 0;
              dfsExplore(part, v, g->edgeList, g->numEdges, dfsVisited,
                         &newServerCount);
              
              initVisited(dfsVisited, numServers);
              for (i=0; i<numOutages; i++) {
                  dfsVisited[outages[i]] = VISITED;
              }
              recordSIDs(part, v, g->edgeList, g->numEdges,
                         dfsVisited, subnetServers, &insertIndex);
              
              /* Now that we know which servers are available post-
                 outage, we can find the longest shortest path by
                 activating Dijkstra's algorithm on the subnetwork */
              dist[v]=0;
              insertionSort(subnetServers, newServerCount);
              dijkstra(subnetServers, newServerCount, g->edgeList,
                       g->numEdges, visited, dist, prev,
                       &currDiameter,
                       solution->postOutageDiameterSIDs, numServers);
              newServerCount=0;
          }
      }
      solution->postOutageDiameter = currDiameter;
      solution->postOutageDiameterCount = currDiameter+1;
      
  } else if(part == TASK_7) {
    
      /* Initialize solution array */
      solution->criticalServerSIDs =
      malloc(sizeof(solution->criticalServerSIDs) * numServers);
      assert(solution->criticalServerSIDs);
      
      int currPushOrder = 0;
      int visited[numServers];
      int parent[numServers];
      int pushOrder[numServers];
      int hraPushOrder[numServers];
      int criticalServers[numServers];
      
      /* Initialization */
      int v;
      for (v=0; v<numServers; v++) {
          visited[v] = UNVISITED;
          parent[v] = NULL_INT;
          pushOrder[v] = 0;
          hraPushOrder[v] = INT_MAX;
          criticalServers[v] = FALSE_INT;
      }
      /* Activate DFS traversal on graph to find critical servers
         in each subnetwork */
      int u;
      for (u=0; u<numServers; u++) {
          if (!visited[u]) {
              dfsExploreCritical(u, numServers, g->edgeList,
                                 g->numEdges, visited, parent,
                                 pushOrder, hraPushOrder,
                                 criticalServers, currPushOrder);
          }
      }
      /* Write to solution array and compute critical server count */
      int i, j=0, count=0;
      for (i=0; i<numServers; i++) {
          if (criticalServers[i] == TRUE_INT) {
              (solution->criticalServerSIDs)[j] = i;
              j++;
              count++;
          }
      }
      solution->criticalServerCount = count;
  }
  return solution;
}

void dijkstra(int subnetServers[], int newServerCount,
              struct edge **edgeList, int numEdges, int visited[],
              int dist[], int prev[], int *currDiameter, int* sol,
              int numServers) {
    int u, e;
    while (!subnetAllVisited(visited, subnetServers, newServerCount)) {
        u = getMin(dist, subnetServers, newServerCount, visited);
        visited[u] = VISITED;
        for (e=0; e<numEdges; e++) {
            if (edgeList[e]->start == u && !visited[edgeList[e]->end]
                && dist[u]+1 < dist[edgeList[e]->end]) {
                
                dist[edgeList[e]->end] = dist[u]+1;
                prev[edgeList[e]->end] = u;
                
                /* Keeps track of graph diameter */
                if (*currDiameter < dist[edgeList[e]->end]) {
                    *currDiameter = dist[edgeList[e]->end];
                    retracePath(dist, sol, prev, currDiameter,
                                numServers);
                }
            }
            if (edgeList[e]->end == u && !visited[edgeList[e]->start]
                && dist[u]+1 < dist[edgeList[e]->start]) {
                
                dist[edgeList[e]->start] = dist[u]+1;
                prev[edgeList[e]->start] = u;
                
                /* Keeps track of graph diameter */
                if (*currDiameter < dist[edgeList[e]->start]) {
                    *currDiameter = dist[edgeList[e]->start];
                    retracePath(dist, sol, prev, currDiameter,
                                numServers);
                }
            }
        }
    }
}

void retracePath(int dist[], int *sol, int prev[], int *currDiameter,
                 int numServers) {
    int j, i, dec = *currDiameter;
    for (j=0; j<numServers; j++) {
        if (dist[j] == *currDiameter) {
            i = j;
        }
    }
    for (j=dec; j>=0; j--) {
        if (j == *currDiameter) {
            sol[j] = i;
        } else {
            sol[j] = prev[i];
            i = prev[i];
        }
    }
}

int subnetAllVisited(int visited[], int subnetServers[],
                     int newServerCount) {
    int v;
    for (v=0; v<newServerCount; v++) {
        if (!visited[subnetServers[v]]) {
            return UNVISITED;
        }
    }
    return VISITED;
}

int getMin(int dist[], int subnetServers[], int newServerCount,
           int visited[]) {
    int v;
    int min=0;
    int minPriority = INT_MAX;
    for (v=0; v<newServerCount; v++) {
        if (dist[subnetServers[v]] < minPriority &&
            !visited[subnetServers[v]]) {
            
            min = v;
            minPriority = dist[subnetServers[v]];
        }
    }
    return subnetServers[min];
}

void dfsExploreCritical(int u, int numServers, struct edge **edgeList,
                        int numEdges, int visited[], int parent[],
                        int pushOrder[], int hraPushOrder[],
                        int criticalServers[], int currPushOrder) {
    
    visited[u] = VISITED;
    pushOrder[u] = currPushOrder+1;
    hraPushOrder[u] = pushOrder[u];
    int children = 0;
    
    int v;
    for (v=0; v<numServers; v++) {
        if (adjacency(u, v, edgeList, numEdges)) {
            if (!visited[v]) {
                children++;
                parent[v] = u;
                dfsExploreCritical(v, numServers, edgeList, numEdges,
                                   visited, parent, pushOrder,
                                   hraPushOrder, criticalServers,
                                   currPushOrder+1);
                hraPushOrder[u] = min(hraPushOrder[u],
                                      hraPushOrder[v]);
                if (parent[u] == NULL_INT && children > 1) {
                    criticalServers[u] = TRUE_INT;
                } else if (parent[u] != NULL_INT &&
                           hraPushOrder[v] >= pushOrder[u]) {
                    criticalServers[u] = TRUE_INT;
                }
            } else if (parent[u] != v) {
                hraPushOrder[u] = min(hraPushOrder[u], pushOrder[v]);
            }
        }
    }
}

int adjacency(int u, int v, struct edge **edgeList, int numEdges) {
    int e;
    for (e=0; e<numEdges; e++) {
        if (edgeList[e]->start == u && edgeList[e]->end == v) {
            return TRUE_INT;
        } else if (edgeList[e]->start == v && edgeList[e]->end == u) {
            return TRUE_INT;
        }
    }
    return FALSE_INT;
}

void initVisited(int visited[], int numServers) {
    int i;
    for (i=0; i<numServers; i++) {
        visited[i] = UNVISITED;
    }
}

void dfsExplore(enum problemPart part,
                int v,
                struct edge **edgeList, int numEdges, int visited[],
                int *currServerCount) {
    
    int e; // edge loop counter
    
    /* Mark current server as visited */
    visited[v] = VISITED;
    
    /* Increment subnetwork's server count */
    if (part == TASK_3 || part == TASK_4) {
        /* Increment server count */
        (*currServerCount)++;
    }
    
    for (e=0; e<numEdges; e++) {
        
        /* if v is the start edge AND destination (end) server
           is unvisited */
        if (edgeList[e]->start == v
            && !visited[edgeList[e]->end]) {
            
            dfsExplore(part, edgeList[e]->end, edgeList,
                       numEdges, visited, &(*currServerCount));
        }
        /* if v is the end edge AND destination (start) server
           is unvisited */
        if (edgeList[e]->end == v
            && !visited[edgeList[e]->start]) {
            
            dfsExplore(part, edgeList[e]->start, edgeList,
                       numEdges, visited, &(*currServerCount));
        }
    }
}

void recordSIDs(enum problemPart part, int v, struct edge **edgeList,
                int numEdges, int reVisited[],
                int *insertArr, int *insertIndex) {
    
    int e; // edge loop counter
    
    /* Mark current server as visited */
    reVisited[v] = VISITED;
    
    /* Insert SID into solution array */
    if (part == TASK_3 || part == TASK_4) {
        insertArr[(*insertIndex)] = v;
        (*insertIndex)++;
    }
    
    for (e=0; e<numEdges; e++) {
        
        /* if v is the start edge AND destination (end) server
           is unvisited */
        if (edgeList[e]->start == v
            && !reVisited[edgeList[e]->end]) {
            
            recordSIDs(part, edgeList[e]->end, edgeList,
                       numEdges, reVisited,
                       insertArr, insertIndex);
        }
        /* if v is the end edge AND destination (start) server
           is unvisited */
        if (edgeList[e]->end == v
            && !reVisited[edgeList[e]->start]) {
    
            recordSIDs(part, edgeList[e]->start, edgeList,
                       numEdges, reVisited,
                       insertArr, insertIndex);
        }
    }
}

void insertionSort(int *arr, int size) {
    int i, j;
    for (i=1; i<size; i++) {
        for (j=i-1; j>=0 && arr[j+1]<arr[j]; j--) {
            swap(arr, j);
        }
    }
}

void swap(int *arr, int j) {
    int temp;
    temp = arr[j];
    arr[j] = arr[j+1];
    arr[j+1] = temp;
}

int max(int x, int y) {
    if (x > y) {
        return x;
    }
    return y;
}

int min(int x, int y) {
    if (x > y) {
        return y;
    }
    return x;
}

int inArr(int *arr, int size, int element) {
    int i;
    for (i=0; i<size; i++) {
        if (arr[i] == element) {
            return TRUE_INT;
        }
    }
    return FALSE_INT;
}

/* DEBUGGING FUNCTIONS */
void printArr(int *arr, int size) {
    int i;
    for (i=0; i<size; i++) {
        if (arr[i] != INT_MAX) {
            printf("%d ", arr[i]);
        } else {
            printf("_ ");
        }
    }
    printf("\n");
}
