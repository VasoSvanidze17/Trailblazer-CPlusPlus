/******************************************************************************
 * File: Trailblazer.cpp
 *
 * Implementation of the graph algorithms that comprise the Trailblazer
 * assignment.
 */

#include "Trailblazer.h"
#include "TrailblazerGraphics.h"
#include "TrailblazerTypes.h"
#include "TrailblazerPQueue.h"
#include "set.h"
#include "grid.h"
#include "vector.h"
#include "console.h"
#include "random.h"
#include "foreach.h"
using namespace std;

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
 * search.  Make sure to update both this implementation prototype and the
 * function prototype in Trailblazer.h.
 */
Vector<Loc>
shortestPath(Loc start,  Loc end,
             Grid<double>& world,
             double costFn(Loc from, Loc to, Grid<double>& world),
			double heuristic(Loc start, Loc end, Grid<double>& world))
{
	Grid<Node> worldNodes;	//declare of Grid
	fillWorldNodes(world, worldNodes);	//create new world with node in grid


	TrailblazerPQueue<Loc> pQueue;	//this queue stores location of Nodes

	//this function gives relevant values first node and add in priority queue
	enqueueFirstNodeLoc(world, worldNodes, start, end, pQueue, heuristic);

	Set<Loc> greenNodesLoc;
	Loc currNodeLoc; 

	while (!pQueue.isEmpty())
	{
		currNodeLoc = pQueue.dequeueMin();	//get first min node from priority queue
		worldNodes[currNodeLoc.row][currNodeLoc.col].color = GREEN;	//fill this node with GREEN
		greenNodesLoc.add(currNodeLoc);
		colorCell(world, currNodeLoc, GREEN);

		if (currNodeLoc == end) break; //if current Node is end Node

		checkAllNeighbors(world, worldNodes, pQueue, greenNodesLoc, 
							currNodeLoc, end, costFn, heuristic);
	}

	Vector<Loc> result;
	recCreatePath(result, worldNodes, start, end);

	return result;
}

Set<Edge> createMaze(int numRows, int numCols)
{	
	Set<Edge> result;
	Grid<Loc> gridGraph(numRows, numCols);
	HashMap<Loc, HashSet<Loc> > clusters;
	TrailblazerPQueue<Edge> edgesQueue;
	createRandomGridGraph(gridGraph, clusters, edgesQueue);

	Edge currEdge;
	int clustersNum = clusters.size();

	while (true)
	{
		currEdge = edgesQueue.dequeueMin();

		if (!clusters[currEdge.start].contains(currEdge.end))
		{	
			//recMerging(clusters, currEdge.start, currEdge.end);
			/*
				recMerging this function doesn't properly. if this function will be 
				implemented correctly then next 3 line must be deleted for correct working
			*/
			clusters[currEdge.start].add(currEdge.end);
			clusters.remove(currEdge.end);
			clustersNum--;
		}

		if (clustersNum == 1) break;
	}

	return result;
}



/*private functions implementation*/

/*
	function fills Grid with Nodes,
	also function paints all Nodes with GRAY and 
	gives them relevant location
*/
void fillWorldNodes(Grid<double>& world, Grid<Node>& worldNodes)
{
	worldNodes.resize(world.nRows, world.nCols);

	Node node;
	node.color = GRAY;
	node.distance = 0;

	for (int i = 0; i < worldNodes.nRows; i++)
	{
		for (int j = 0; j < worldNodes.nCols; j++)
		{
			node.loc = makeLoc(i, j);
			worldNodes[i][j] = node;
		}
	}
}

/*
	function gives relevant values first node and add in priority queue
*/
void enqueueFirstNodeLoc(Grid<double>& world, Grid<Node>& worldNodes,
	Loc startLoc, Loc endLoc,
	TrailblazerPQueue<Loc>& pQueue,
	double heuristic(Loc start, Loc end, Grid<double>& world))
{
	worldNodes[startLoc.row][startLoc.col].color = YELLOW;	//paint first Node with YELLOW
	worldNodes[startLoc.row][startLoc.col].distance =
		heuristic(startLoc, endLoc, world);

	//add in priority queue
	pQueue.enqueue(startLoc, worldNodes[startLoc.row][startLoc.col].distance);
}

/*
	function checks all neighbors node and
	if neighbor node is gray then this node will be yellow and
	function adds in priority queue and if neighbor node is yellow then
	function compares current distance and strored distance and
	if current node length is less then function updates priority of node
*/
void checkAllNeighbors(Grid<double>& world, Grid<Node>& worldNodes,
						TrailblazerPQueue<Loc>& pQueue,
						Set<Loc>& greenNodesLoc, Loc nodeLoc, Loc endLoc,
						double costFn(Loc from, Loc to, Grid<double>& world),
						double heuristic(Loc start, Loc end, Grid<double>& world))
{
	for (int i = nodeLoc.row - 1; i <= nodeLoc.row + 1; i++)
	{
		for (int j = nodeLoc.col - 1; j <= nodeLoc.col + 1; j++)
		{
			if (worldNodes.inBounds(i, j) && !greenNodesLoc.contains(makeLoc(i, j)))
			{
				if (worldNodes[i][j].color == GRAY)
				{
					enqueueNewYellowNodeLoc(world, worldNodes, pQueue, nodeLoc,
											 makeLoc(i, j), endLoc, costFn, heuristic);
				}

				if (worldNodes[i][j].color == YELLOW &&
					worldNodes[i][j].distance > worldNodes[nodeLoc.row][nodeLoc.col].distance +
					costFn(nodeLoc, makeLoc(i, j), world))
				{
					updatePriority(world, worldNodes, pQueue, nodeLoc,
									makeLoc(i, j), endLoc, costFn, heuristic);
				}
			}
		}
	}
}

/*
	function paints gray nodes with yellow and adds in priority queue
*/
void enqueueNewYellowNodeLoc(Grid<double>& world, Grid<Node>& worldNodes,
								TrailblazerPQueue<Loc>& pQueue,
								Loc nodeLoc, Loc neighborNodeLoc, Loc endLoc,
								double costFn(Loc from, Loc to, Grid<double>& world),
								double heuristic(Loc start, Loc end, Grid<double>& world))
{
	//paint node with yellow
	worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].color = YELLOW; 
	colorCell(world, makeLoc(neighborNodeLoc.row, neighborNodeLoc.col), YELLOW);

	//calculate distance value and add in pQueue
	worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].distance =
		worldNodes[nodeLoc.row][nodeLoc.col].distance +
		costFn(nodeLoc, neighborNodeLoc, world) +
		heuristic(neighborNodeLoc, endLoc, world);
	pQueue.enqueue(neighborNodeLoc, worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].distance);

	//save previous node location
	worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].preNodeLoc = nodeLoc;
}

/*
	function updates priority in queue and also update
	previous node and distance
*/
void updatePriority(Grid<double>& world, Grid<Node>& worldNodes,
					TrailblazerPQueue<Loc>& pQueue,
					Loc nodeLoc, Loc neighborNodeLoc, Loc endLoc,
					double costFn(Loc from, Loc to, Grid<double>& world),
					double heuristic(Loc start, Loc end, Grid<double>& world))
{
	double tempdistance = worldNodes[nodeLoc.row][nodeLoc.col].distance +
		costFn(nodeLoc, neighborNodeLoc, world);
	
	if (tempdistance < worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].distance)
	{
		/*
		update value of distance from start node to current node and
		also update this value in pQeueue with decreaseKey function
		*/
		worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].distance = tempdistance;

		pQueue.decreaseKey(neighborNodeLoc, 
			worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].distance +
			heuristic(neighborNodeLoc, endLoc, world));

		//update previous location of node
		worldNodes[neighborNodeLoc.row][neighborNodeLoc.col].preNodeLoc = nodeLoc;
	}
}

/*
	function adds location of nodes in vector so creates path recursively  
*/
void recCreatePath(Vector<Loc>& vector, Grid<Node>& worldNodes, Loc start, Loc end)
{
	if (end == start)
	{
		vector.add(end); //add start location of node
		return;
	}

	recCreatePath(vector, worldNodes, start, 
				worldNodes[end.row][end.col].preNodeLoc); //add previous location of nodes
	vector.add(end);	//add last location of node
}

void createRandomGridGraph(Grid<Loc>& gridGraph, HashMap<Loc, HashSet<Loc> >& clusters,
							TrailblazerPQueue<Edge>& edgesQueue)
{
	Loc startLocEdge;
	HashSet<Loc> emptyClusters;

	Edge currEdge;
	double edgeWeight;

	for (int i = 0; i < gridGraph.nRows; i++)
	{
		for (int j = 0; j < gridGraph.nCols; j++)
		{
			startLocEdge = makeLoc(i, j);
			gridGraph[i][j] = startLocEdge;
			clusters.put(startLocEdge, emptyClusters);

			if (gridGraph.inBounds(i, j + 1))
			{
				currEdge = makeEdge(startLocEdge, makeLoc(i, j + 1));
				edgeWeight = randomReal(0, 1);
				edgesQueue.enqueue(currEdge, edgeWeight);
			}

			if (gridGraph.inBounds(i + 1, j))
			{
				currEdge = makeEdge(startLocEdge, makeLoc(i + 1, j));
				edgeWeight = randomReal(0, 1);
				edgesQueue.enqueue(currEdge, edgeWeight);
			}
		}
	}
}

void recMerging(HashMap<Loc, HashSet<Loc> >& clusters, Loc startLoc, Loc endLoc)
{
	if (!clusters.containsKey(endLoc)) return;
	if (clusters[endLoc].size() == 0)
	{
		clusters[startLoc].add(endLoc);
		clusters.remove(endLoc);
		return;
	}

	foreach (Loc currLoc in clusters[endLoc])
	{
		recMerging(clusters, startLoc, currLoc);
	}

	clusters[startLoc].add(endLoc);
	clusters.remove(endLoc);
}