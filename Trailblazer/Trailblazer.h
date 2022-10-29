/******************************************************************************
 * File: Trailblazer.h
 *
 * Exports functions that use Dijkstra's algorithm, A* search, and Kruskal's
 * algorithm as specified in the assignment handout.
 */

#ifndef Trailblazer_Included
#define Trailblazer_Included

#include "TrailblazerTypes.h"
#include "set.h"
#include "grid.h"
#include "TrailblazerPQueue.h"
#include "hashmap.h"
#include "hashset.h"

/* Function: shortestPath
 * 
 * Finds the shortest path between the locations given by start and end in the
 * specified world.	 The cost of moving from one edge to the next is specified
 * by the given cost function.	The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.	If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this function prototype and the
 * implementation inside of Trailblazer.cpp.
 */
Vector<Loc>
shortestPath(Loc start,
             Loc end,
             Grid<double>& world,
             double costFn(Loc from, Loc to, Grid<double>& world),
			double heuristic(Loc start, Loc end, Grid<double>& world));

/* Function: createMaze
 * 
 * Creates a maze of the specified dimensions using a randomized version of
 * Kruskal's algorithm, then returns a set of all of the edges in the maze.
 *
 * As specified in the assignment handout, the edges you should return here
 * represent the connections between locations in the graph that are passable.
 * Our provided starter code will then use these edges to build up a Grid
 * representation of the maze.
 */
Set<Edge> createMaze(int numRows, int numCols);

struct Node
{
	Loc loc;		//location of this node in grid
	Color color;	//color of this node
	double distance;	//distance from start node to this node
	Loc preNodeLoc;	//location of previous node in grid
};

void enqueueFirstNodeLoc(Grid<double>& world, Grid<Node>& worldNodes,
							Loc startLoc, Loc endLoc,
							TrailblazerPQueue<Loc>& pQueue, 
							double heuristic(Loc start, Loc end, Grid<double>& world));

void fillWorldNodes(Grid<double>& world, Grid<Node>& worldNodes);

void checkAllNeighbors(Grid<double>& world, Grid<Node>& worldNodes, 
						TrailblazerPQueue<Loc>& pQueue, 
						Set<Loc>& greenNodesLoc, Loc nodeLoc, Loc endLoc,
						double costFn(Loc from, Loc to, Grid<double>& world),
						double heuristic(Loc start, Loc end, Grid<double>& world));

void enqueueNewYellowNodeLoc(Grid<double>& world, Grid<Node>& worldNodes, 
								TrailblazerPQueue<Loc>& pQueue,
								Loc nodeLoc, Loc neighborNodeLoc, Loc endLoc,
								double costFn(Loc from, Loc to, Grid<double>& world),
								double heuristic(Loc start, Loc end, Grid<double>& world));

void updatePriority(Grid<double>& world, Grid<Node>& worldNodes,
					TrailblazerPQueue<Loc>& pQueue,
					Loc nodeLoc, Loc neighborNodeLoc, Loc endLoc,
					double costFn(Loc from, Loc to, Grid<double>& world),
					double heuristic(Loc start, Loc end, Grid<double>& world));

void recCreatePath(Vector<Loc>& vector, Grid<Node>& worldNodes, Loc start, Loc end);

void createRandomGridGraph(Grid<Loc>& gridGraph, HashMap<Loc, HashSet<Loc> >& clusters,
							TrailblazerPQueue<Edge>& edgesQueue);

void recMerging(HashMap<Loc, HashSet<Loc> >& clusters, Loc startLoc, Loc endLoc);

#endif
