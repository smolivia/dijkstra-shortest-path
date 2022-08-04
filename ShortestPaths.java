package dijkstra-shortest-path;

package spath;

import java.util.LinkedList;
import java.util.HashMap;

// heap-related structures from Lab 3
import heaps.PQEntry;
import heaps.MinHeap;

// directed graph structure
import spath.graphs.DirectedGraph;
import spath.graphs.Edge;
import spath.graphs.Vertex;

import timing.Ticker;


public class ShortestPaths {
    // "infinity" value for path lengths
    private final static Integer inf = Integer.MAX_VALUE;
    
    // a directed graph, and a weighting function on its edges
    private final DirectedGraph g;
    private HashMap<Edge, Integer> weights;	
    
    // starting vertex for shortest path computation
    private Vertex startVertex;
    
    // map from vertices to their handles in the priority queue
    private HashMap<Vertex, PQEntry<VertexAndDist, Integer>> handles;
    
    // map from vertices to their parent edges in the shortest-path tree
    private HashMap<Vertex, Edge> parentEdges;
    
    
    //
    // constructor
    //
    public ShortestPaths(DirectedGraph g, HashMap<Edge,Integer> weights, 
			 Vertex startVertex) {
    	this.g           = g;
    	this.weights     = weights;

    	this.startVertex = startVertex;	
	
    	this.handles     = new HashMap<Vertex, PQEntry<VertexAndDist, Integer>>();
    	this.parentEdges = new HashMap<Vertex, Edge>();
    }

    
    //
    // run() 
    //
    // Given a weighted digraph stored in g/weights, compute a
    // shortest-path tree of parent edges back to a given starting
    // vertex.
    //
    public void run() {
    	Ticker ticker = new Ticker(); // heap requires a ticker
	
    	MinHeap<VertexAndDist, Integer> pq = 
    			new MinHeap<VertexAndDist, Integer>(g.getNumVertices(), ticker);
	
    	//
    	// Put all vertices into the heap, infinitely far from start.
    	// Record handle to each inserted vertex, and initialize
    	// parent edge of each to null (since we have as yet found 
    	// no path to it.)
    	//
    	for (Vertex v : g.vertices()) {
    		PQEntry<VertexAndDist, Integer> d = pq.insert(new VertexAndDist(v, inf), inf);
    		handles.put(v, d);
    		parentEdges.put(v, null);
    	}
	
    	//
    	// Relax the starting vertex's distance to 0.
    	//   - get the handle to the vertex from the heap
    	//   - extract the vertex + distance object from the handle
    	//   - update the distance in the vertex + distance object
    	//   - update the heap through the vertex's handle
    	//
    	PQEntry<VertexAndDist, Integer> startHandle = handles.get(startVertex);
    	VertexAndDist vd = startHandle.getElement();
    	vd.setDistance(0);
    	startHandle.updatePriority(0);
	
    	//CODE BEGINS
    	while(!pq.isEmpty()) {
    		VertexAndDist current = pq.extractMin().getElement(); //get vertex with minimum distance from priority queue
    		Iterable<Edge> outgoing = current.getVertex().edgesFrom();
    		for(Edge e : outgoing) { //iterate through outgoing edges of minimum vertex
    			Vertex adj = e.to; //get vertex adjacent to current edge

    			PQEntry<VertexAndDist, Integer> adjHandle = handles.get(adj);
    			VertexAndDist adjAndDist = adjHandle.getElement(); 
    			int adjDist = adjAndDist.getDistance(); // get adjacent vertex's distance through handle

    			int newDist = current.getDistance() + weights.get(e); //calculate current vertex's distance + weight of edge going to the adjacent vertex
    			if(newDist < adjDist) { //if the newly calculated distance is less than the adjacent vertex's current distance,
										//update the adjacent vertex's distance to new distance and parent edge to current edge
    				adjAndDist.setDistance(newDist);
    				parentEdges.put(adj, e);
    				adjHandle.updatePriority(newDist); //update the priority of the adjacent vertex so it can be correctly recognized by the pq
    			} 			
    		}
    	//CODE ENDS
    		
    	}
       
    }
    
    
    //
    // returnPath()
    //
    // Given an ending vertex v, compute a linked list containing every
    // edge on a shortest path from the starting vertex (stored) to v.
    // The edges should be ordered starting from the start vertex.
    //
    public LinkedList<Edge> returnPath(Vertex endVertex) {
    	LinkedList<Edge> path = new LinkedList<Edge>();
        //CODE BEGINS
    	Vertex v = endVertex; //start at the ending vertex
    	while(parentEdges.get(v) != null) { //continue to iterate through the parent edges and add them to the path
											//in reverse order until we reach the starting vertex with a null parent edge
    		path.addFirst(parentEdges.get(v));
    		v = parentEdges.get(v).from; //update the current vertex by following the parent edge
    	}
    	return path;
        //CODE ENDS
    }
    
    ////////////////////////////////////////////////////////////////
    //BELOW ARE HELPER METHODS PRE-PROVIDED FOR THE PROJECT:
    //
    // returnLength()
    // Compute the total weight of a putative shortest path
    // from the start vertex to the specified end vertex.
    // No user-serviceable parts inside.
    //
    public int returnLength(Vertex endVertex) {
    	LinkedList<Edge> path = returnPath(endVertex);
	
    	int pathLength = 0;
    	for(Edge e : path) {
    		pathLength += weights.get(e);
    	}
	
    	return pathLength;
    }

    //
    // returnLengthDirect()
    // Expose the current-best distance estimate stored at a vertex.
    // Useful for comparing to ground-truth shortest-path distance
    //   in the absence of parent pointers.
    public int returnLengthDirect(Vertex endVertex) {
    	PQEntry<VertexAndDist, Integer> endHandle = handles.get(endVertex);
    	return endHandle.getElement().getDistance();
    }
}
