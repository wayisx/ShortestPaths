package de.unikiel.npr.thorup.algs;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Hashtable;
import java.util.ListIterator;

import de.unikiel.npr.thorup.ds.ArrayPriorityQueue;
import de.unikiel.npr.thorup.ds.PriorityQueue;
import de.unikiel.npr.thorup.ds.PriorityQueueItem;
import de.unikiel.npr.thorup.ds.graph.Edge;
import de.unikiel.npr.thorup.ds.graph.Graph;
import de.unikiel.npr.thorup.ds.graph.WeightedEdge;

/**
 * An implementation of <i>Dijkstra</i>'s single-source shortest paths
 * algorithm.
 * 
 * @author
 * 		<a href="mailto:npr@informatik.uni-kiel.de">Nick Pr&uuml;hs</a>
 * @version
 * 		1.0, 09/17/09
 */
public class NewKColorDijkstra {
	/**
	 * The predecessors of all vertices of the checked graph on their way to
	 * the source vertex.
	 */
	private int[] predecessors;
	
	/**
	 * The distances of all vertices of the checked graph from the source
	 * vertex.
	 */
	private int[] distances;
	
	Hashtable EtS;
	Hashtable CEt;
	int[] colors;
	PriorityQueueItem<Integer>[] items;

	
	
	public void initialize(Graph<? extends Edge> g, int u,  
			 List lengthDist, PriorityQueue q) throws IllegalArgumentException {
		
		//Create tables to store distinct edges(EtS) and CurrentEdges
				EtS = new Hashtable(lengthDist.size()*2);
			//	CEt = new Hashtable(lengthDist.size()*2);
				
				
				for(int i=0; i<lengthDist.size(); i++){
					int keyweights = (Integer) lengthDist.get(i);
					ArrayList edges = new ArrayList();
					EtS.put(keyweights, edges);
					//CEt.put(keyweights, edges);
				}
		// check arguments
				if (g == null) {
					String errorMessage = "The passed graph musn't be null.";
					
					throw new IllegalArgumentException(errorMessage);
				}
				
				if (u < 0 || u >= g.getNumberOfVertices()){
					String errorMessage = "The vertex with index " + u +
						" is not within the passed graph.";
					throw new IllegalArgumentException(errorMessage);
				}
				
				
				// get the number of vertices of the passed graph
				int n = g.getNumberOfVertices();
				int blackV = 0;
				// initialize help arrays
				colors = new int[n];
				predecessors = new int[n];
				distances = new int[n];
				//Create priority queue for vertices
				items = new PriorityQueueItem[n];
				
				int white = 0;//unvisited
				int gray = 1; //visited, tentative distance.
				int black = 2;//visited, distance confirmed.
				
				for (int v = 0; v < n; v++) {
					colors[v] = white;
					predecessors[v] = -1;
					distances[v] = -1;
				}
				
				// mark the start vertex as visited
				colors[u] = black;
				blackV++;
				distances[u] = 0;
				//Get neighbors
				int[] ne = g.getArrayOfAdjacentVertices(u);
				int[] we = g.getArrayOfIncidentEdgeWeights(u);
			
				for (int i = 0; i < ne.length; i++) {
					//get neighbor
					int neighbor =  ne[i];
					int weights = we[i];
					ArrayList edges = (ArrayList) EtS.get(weights);
					edges.add(new WeightedEdge(u,neighbor,weights));
					EtS.put(weights, edges);
					/*
					ArrayList CEtedges = (ArrayList) CEt.get(weights);
					if(CEtedges.isEmpty()){
						CEtedges.add(new WeightedEdge(u,neighbor,weights));
						CEt.put(weights, CEtedges);
					}*/
				}
				
				items = new PriorityQueueItem[lengthDist.size()];
				
				
				for(int i=0; i<lengthDist.size(); i++){
					int keyweights = (Integer) lengthDist.get(i);
					PriorityQueueItem<Integer> item = q.insert(keyweights, Integer.MAX_VALUE);
					items[i] = item;
					update(i, keyweights, lengthDist, q);
				}
				
				
				//New-Dijkstra
				while (blackV<n){
					//System.out.println(blackV);
					PriorityQueueItem ft = q.findMin();
					int r = (Integer) ft.getItem();
					ArrayList edges = (ArrayList) EtS.get(r);
					WeightedEdge edge = (WeightedEdge) edges.get(0);
					int i = edge.getSource();
					int j = edge.getTarget();
					distances[j] = distances[i] + r;
					predecessors[j] = i;
					colors[j] = 2;
					blackV++;
					
					ne = g.getArrayOfAdjacentVertices(j);
					we = g.getArrayOfIncidentEdgeWeights(j);
				
					for (int nei = 0; nei < ne.length; nei++) {
						//get neighbor
						int neighbor =  ne[nei];
						if(colors[neighbor]!=2){
							int weights = we[nei];
							ArrayList Nedges = (ArrayList) EtS.get(weights);
							Nedges.add(new WeightedEdge(j,neighbor,weights));
							EtS.put(weights, Nedges);
						}
					}
					q = new ArrayPriorityQueue<Integer>(n);
					for(int t=0; t<lengthDist.size(); t++){
						int keyweights = (Integer) lengthDist.get(t);
						//PriorityQueueItem<Integer> item = q.insert(keyweights, Integer.MAX_VALUE);
						//items[t] = item;
						update(t, keyweights, lengthDist, q);
					}
				}
				
				
	}
	
	
	public void update(int pos, int keyweights,  List lengthDist, PriorityQueue q) throws IllegalArgumentException {
	
		ArrayList EtSedges = (ArrayList) EtS.get(keyweights);
		Iterator k = EtSedges.iterator();
	
		if(!k.hasNext()){
			return;
		}
		
		WeightedEdge we = (WeightedEdge) k.next();
		int i = we.getSource();
		int j = we.getTarget();
		if(colors[j] != 2){
			int ft = distances[i] + keyweights;
			//q.updateKeyTo(items[pos] , ft);
			q.insert(keyweights, ft);
			return;
		}
		
	
		while(colors[j]==2){
		//	k.remove(we);
			k.remove();
			if(k.hasNext()){
				we = (WeightedEdge)k.next();
				j = we.getTarget();
			}else{
				break;
			}
			
		}
		
		if(colors[j] != 2){
			i = we.getSource();
			int ft = distances[i] + keyweights;
			q.insert(keyweights, ft);
		}
	}
	

	/**
	 * Iterates the passed weighted graph, computing the paths from the passed
	 * source vertex to all others, including their corresponding distances.
	 * It is imperative that u is a vertex of g. The results of the algorithm
	 * can be accessed using {@link #getDistances()} and
	 * {@link #getPredecessors()} as soon as the algorithm has finished.
	 * 
	 * @param g
	 * 		the graph to run the algorithm on
	 * @param u
	 * 		the source vertex
	 * @param q 
	 * 		the priority queue used for maintaining the order the vertices are
	 * 		visited in
	 * @throws IllegalArgumentException
	 * 		if the passed graph is <code>null</code>
	 * @throws IllegalArgumentException
	 * 		if <code>u</code> is not a vertex of <code>g</code>
	 * @see #getDistances()
	 * @see #getPredecessors()
	 */
	/*
	@SuppressWarnings("unchecked")
	public void findShortestPaths(Graph<? extends Edge> g, int u,  
			 List lengthDist, PriorityQueue q) throws IllegalArgumentException {
		
		
		
		
		// check arguments
		if (g == null) {
			String errorMessage = "The passed graph musn't be null.";
			
			throw new IllegalArgumentException(errorMessage);
		}
		
		if (u < 0 || u >= g.getNumberOfVertices()){
			String errorMessage = "The vertex with index " + u +
				" is not within the passed graph.";
			throw new IllegalArgumentException(errorMessage);
		}
		//Create tables to store distinct edges(EtS) and CurrentEdges
		
		
		
		// get the number of vertices of the passed graph
		int n = g.getNumberOfVertices();
		
		// initialize help arrays
		int[] colors = new int[n];
		predecessors = new int[n];
		distances = new int[n];
		
		int white = 0;//unvisited
		int gray = 1; //visited, tentative distance.
		int black = 2;//visited, distance confirmed.
		
		for (int v = 0; v < n; v++) {
			colors[v] = white;
			predecessors[v] = -1;
			distances[v] = -1;
		}
		
		// mark the start vertex as visited
		colors[u] = gray;
		
		// initialize the priority queue
		PriorityQueueItem<Integer>[] items = new PriorityQueueItem[n];
		
		PriorityQueueItem<Integer> item = q.insert(u, 0);
		items[u] = item;
		
		// extend distance tree until border is empty
		while (!q.isEmpty()) {
			// get next vertex for the distance tree and set its distance
			item = q.deleteMin();
			int d = (int)item.getKey();
			int v = (int)item.getItem();
			distances[v] = d;
			
			// mark that vertex as visited
			colors[v] = black;
			
			// update border and border approximation
			int[] a = g.getArrayOfAdjacentVertices(v);
			int[] da = g.getArrayOfIncidentEdgeWeights(v);
			
			for (int i = 0; i < a.length; i++) {
				// get next neighbor
				int w = a[i];
				
				// update approximation
				int dw = d + da[i];
				
				// update entry if necessary
				if (colors[w] == white ||
					(colors[w] == gray && items[w].getKey() > dw)) {
					colors[w] = gray;
					predecessors[w] = v;
					
					if (items[w] == null) {
						items[w] = q.insert(w, dw);
					} else {
						q.decreaseKeyTo(items[w], dw);
					}
				}
			}
		}
	}
	*/
	
	
	/**
	 * Gets the predecessors of all vertices of the checked graph on their way
	 * to the source vertex.
	 * 
	 * @return
	 * 		the predecessors of all vertices of the checked graph on their way
	 * 		to the source vertex
	 */
	public int[] getPredecessors() {
		return predecessors;
	}

	/**
	 * Gets the distances of all vertices of the checked graph from the source
	 * vertex.
	 * 
	 * @return
	 * 		the distances of all vertices of the checked graph from the source
	 * 		vertex
	 */
	public int[] getDistances() {
		return distances;
	}
}
