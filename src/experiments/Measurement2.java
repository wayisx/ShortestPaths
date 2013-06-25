package experiments;

import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Scanner;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.io.ParseException;
import com.vividsolutions.jts.io.WKTReader;

import de.unikiel.npr.thorup.algs.Dijkstra;
import de.unikiel.npr.thorup.algs.Kruskal;
import de.unikiel.npr.thorup.algs.Thorup2;
import de.unikiel.npr.thorup.ds.ArrayPriorityQueue;
import de.unikiel.npr.thorup.ds.FibonacciHeap;
import de.unikiel.npr.thorup.ds.SplitFindminStructureGabow;
import de.unikiel.npr.thorup.ds.UnionFindStructureTarjan;
import de.unikiel.npr.thorup.ds.graph.AdjacencyListWeightedDirectedGraph;
import de.unikiel.npr.thorup.ds.graph.WeightedEdge;

/**
 * An abstract series of measurent for running the algorithms by <i>Dijkstra</i>
 * and <i>Thorup</i> several times, taking the average over all passes, and
 * increasing a single value step by step.
 * 
 * @author
 * 		<a href="mailto:npr@informatik.uni-kiel.de">Nick Pr&uuml;hs</a>
 * @version
 * 		1.0, 09/17/09
 */
public abstract class Measurement2 {
	/**
	 * The number of passes to take the average of.
	 */
	int numberOfPasses = 1;
	
	/**
	 * The maximum number of vertices of this series of measurement. 
	 */
	int numberOfVertices;
	
	/**
	 * The number of vertices of the graph the current performance test runs on.
	 */
	int numberOfVerticesCurrent;
	
	/**
	 * The maximum number of edges per vertex of this series of measurement.
	 */
	int numberOfEdgesPerVertex;
	
	/**
	 * The maximum edge weight of this series of measurement.
	 */
	int maximumEdgeWeight;
	
	/**
	 * The value of the current performance test which is increased step by
	 * step.
	 */
	int currentStepValue;
	
	/**
	 * The maximum value of all performance tests.
	 */
	int maximumStepValue;
	
	/**
	 * All values that have been used so far for the performance tests.
	 */
	long[] allStepValues;
	
	/**
	 * The number of the current iteration of this series of measurement.
	 */
	int currentStep;
	
	/**
	 * The total number of the iterations of this series of measurement.
	 */
	int numberOfSteps;
	
	/**
	 * The difference, measured in milliseconds, between the start time of the
	 * current performance test and midnight, January 1, 1970 UTC.
	 */
	double start;
	
	/**
	 * The difference, measured in milliseconds, between the stop time of the
	 * current performance test and midnight, January 1, 1970 UTC.
	 */
	double stop;

	/**
	 * The running times of the most recent performance tests to compute the
	 * average of.
	 */
	double[] timesToComputeTheAverageOf;
	
	/**
	 * The running times of <i>Dijkstra</i>'s algorithm using an array heap
	 * in the most recent performance tests.
	 */
	double[] timesDijkstraArrayHeap;
	
	/**
	 * The running times of <i>Dijkstra</i>'s algorithm using a Fibonacci heap
	 * in the most recent performance tests.
	 */
	double[] timesDijkstraFibHeap;
	
	/**
	 * The times required for constructing the <i>msb</i>-minimum spanning trees
	 * for <i>Thorup</i>'s algorithm in the most recent performance tests.
	 */
	double[] timesThorupMST;
	
	/**
	 * The times required for constructing the component tree, initializing the
	 * bucket structure, and preparing the unvisited data structure for
	 * <i>Thorup</i>'s algorithm in the most recent performance tests.
	 */
	double[] timesThorupDS;
	
	/**
	 * The running times of the visiting part of <i>Thorup</i>'s algorithm in
	 * the most recent performance tests.
	 */
	double[] timesThorupVisit;
	
	/**
	 * The scanner used to read the user input.
	 */
	Scanner in;
	
	/**
	 * Whether all values and running times are written to the console again
	 * after this series of measurement.
	 */
	boolean writeTableColumns = true;
	
	DecimalFormat df = new DecimalFormat("0.0000"); 
	/**
	 * Runs the algorithms by <i>Dijkstra</i> and <i>Thorup</i> several times,
	 * taking the average over all passes, and increasing a single value step
	 * by step.
	 * @throws Exception 
	 */
	public void execute() throws Exception{
	//	FileReader fRdr = new FileReader("D:\\Down\\GIS Data\\Shapefiles\\AFG_rds\\AFG_roads.WKT");    
		FileReader fRdr = new FileReader("D:\\Down\\GIS Data\\Shapefiles\\gm-jpn-trans_u_2\\5\\roadl_jpn.WKT");    
		
		WKTReader wktRdr = new WKTReader();
		Geometry geom = wktRdr.read(fRdr);
		fRdr.close();		
		MultiLineString mls=(MultiLineString)geom;
		Hashtable hb=new Hashtable(mls.getNumGeometries()*2);
		Hashtable hbw=new Hashtable(mls.getNumGeometries()*2);
		int vertexNum=0;
		List edges= new ArrayList<WeightedEdge>();
		
		   	for (int i = 0; i < mls.getNumGeometries(); i++) {

	   //	for (int i = 0; i < 1000; i++) {
		   		 LineString ls = (LineString) mls.getGeometryN(i);
		   		 String sp = ls.getStartPoint().toString();
		   		 String ep = ls.getEndPoint().toString();
		   		 int weight =  (int) (ls.getLength()*1000000);
		   		 
		   		 if(weight==0){
		   			 System.out.println(weight);
		   		 }
		   		 int sEdge=0;
		   		 int eEdge=0;
		   		 Object s = hb.get(sp);
		   		 Object e = hb.get(ep);
		   		 if(s==null){
		   			 hb.put(sp, vertexNum);
		   			 sEdge = vertexNum;
				     vertexNum++; 
		   		 }else{
		   			sEdge =  Integer.parseInt(s.toString());
		   		 }
		   		 if(e==null){
		   			 hb.put(ep, vertexNum);
		   			 eEdge = vertexNum;
				     vertexNum++; 
		   		 }else{
		   			 eEdge =  Integer.parseInt(e.toString()); 
		   		 }
		   		//WeightedEdge we =  new WeightedEdge(sEdge, eEdge);
		   		 String we = sEdge+"-"+eEdge;
		   		if(hbw.get(we)==null){
		   		 hbw.put(we, weight);
		   		
		   		 edges.add(new WeightedEdge(sEdge, eEdge));
		   		}else{ //There are more than one way between two vertexes.
		   			int w = Integer.parseInt(hbw.get(we).toString());
		   			if(w>weight){
		   				hbw.remove(we);
		   				hbw.put(we, weight);
		   			}
		   		}
			}
		   
		   	AdjacencyListWeightedDirectedGraph<WeightedEdge> graph =
				new AdjacencyListWeightedDirectedGraph<WeightedEdge>(vertexNum);
		   	for(int i = 0; i<edges.size(); i++){
		   		WeightedEdge we = (WeightedEdge) edges.get(i);

		   		 int sEdge=we.getSource();
		   		 int eEdge=we.getTarget();
		   		 int weight = Integer.parseInt(hbw.get(new String(sEdge+"-"+eEdge)).toString());
		   		graph.addEdge(new WeightedEdge(sEdge,eEdge,weight));
		   		graph.addEdge(new WeightedEdge(eEdge,sEdge,weight));
		   	}
		/*
		
		   	AdjacencyListWeightedDirectedGraph<WeightedEdge> graph =
				generateGraph();
		   int	vertexNum=10000;*/
		System.out.println(" Graph has " +graph.getNumberOfEdges() + " edges,"+ vertexNum+" vertexes.");
		System.out.println("");
		
		numberOfPasses =1;
		//maximumStepValue and  currentStepValue controls the times of loop
		currentStepValue = 1;
		maximumStepValue = 30;
		numberOfVerticesCurrent = vertexNum;
		// initialize all variables
		numberOfSteps = maximumStepValue / currentStepValue;
		
		timesToComputeTheAverageOf = new double[numberOfPasses];
		
		timesDijkstraArrayHeap = new double[numberOfSteps];
		timesDijkstraFibHeap = new double[numberOfSteps];
		timesThorupMST = new double[numberOfSteps];
		timesThorupDS = new double[numberOfSteps];
		timesThorupVisit = new double[numberOfSteps];
		
		allStepValues = new long[numberOfSteps];
		int sourceID = 0;
		System.out.println();
		
		// run this series of measurement until the maximum step value
	
		
		
	
		
		
		
		while (currentStepValue <= maximumStepValue) {
			numberOfVerticesCurrent = vertexNum;
			// remember the current step value for evaluation purposes
			allStepValues[currentStep] = currentStepValue;
			
			// generate the next test instance
		//	System.out.println();
		//	System.out.print("Generating a random graph with " +
				//	numberOfVerticesCurrent + " vertices...");
			
			
			
			
			Thorup2 thorup = new Thorup2();
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				thorup.constructMinimumSpanningTree(graph,
						new Kruskal(new UnionFindStructureTarjan<Integer>()));
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}

			// compute the average
			timesThorupMST[currentStep] =
				getAverage(timesToComputeTheAverageOf);
					
			// show the result
		//	System.out.print(" Took " + timesThorupMST[currentStep] + " ms for constructing the MST, ");
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				thorup.constructOtherDataStructures
					(new UnionFindStructureTarjan<Integer>());
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}

			// compute the average
			timesThorupDS[currentStep] =
				getAverage(timesToComputeTheAverageOf);
					
			// show the result
		//	System.out.print(timesThorupDS[currentStep] +" ms for constructing the other data structures,");
			
			//System.out.println();
			
			
			
			
			
			
			
			
			
			
			// run Dijkstra's algorithm with an array heap and take the time
			if(currentStepValue==1){
				System.out.print("Running Dijkstra with an array priority " +
				"queue...");
			}
			
			//sourceID=0;
			Dijkstra dijkstra = new Dijkstra();
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				dijkstra.findShortestPaths(graph, sourceID,
						new ArrayPriorityQueue<Integer>
							(numberOfVerticesCurrent));
				stop = System.nanoTime();
				//sourceID++;
				timesToComputeTheAverageOf[pass] = stop - start;
			}
	
			// compute the average
			timesDijkstraArrayHeap[currentStep] =
				getAverage(timesToComputeTheAverageOf);
					
			// show the result
			if(currentStepValue==1){
			System.out.println(" took " + timesDijkstraArrayHeap[currentStep] +
					" ms (average of " + numberOfPasses + " passes).");
			// run Dijkstra's algorithm with a Fibonacci heap and take the time
			System.out.print("Running Dijkstra with a Fibonacci heap...");
			}
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				dijkstra.findShortestPaths(graph, sourceID,
						new FibonacciHeap<Integer>());
				stop = System.nanoTime();
				//sourceID++;
				timesToComputeTheAverageOf[pass] = stop - start;
			}
			int [] dd = dijkstra.getDistances();
			/*
			for(int i = 0; i<numberOfVerticesCurrent;i++){
				
				System.out.println(i+"-> "+dd[i]);
			}
			*/
	
			// compute the average
			timesDijkstraFibHeap[currentStep] =
				getAverage(timesToComputeTheAverageOf);
					
			// show the result
			if(currentStepValue==1){
			System.out.println(" took " + timesDijkstraFibHeap[currentStep] +
					" ms (average of " + numberOfPasses + " passes).");
	
			
			/*
			 * construct the msb-minimum spanning tree for Thorup's algorithm
			 * and take the time
			 */
			System.out.print("Running Thorup...");
			
			}
			
			
			/*
			 * construct all other data structures for Thorup's algorithm and
			 * take the time
			 */
			
			
			
			// run Thorup's algorithm and take the time
			int[] td = new int[numberOfVerticesCurrent];
			for (int pass = 0; pass < numberOfPasses; pass++) {
				
				if (currentStep > 1) {
					thorup.cleanUpBetweenQueries();
				}

				start = System.nanoTime();
				 td = thorup.findShortestPaths(sourceID);
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}
			
			/*
			for(int i = 0; i<td.length;i++){
				
				System.out.println(i+"-> "+td[i]);
			}*/
			// compute the average
			timesThorupVisit[currentStep] =
				getAverage(timesToComputeTheAverageOf);
		
					
			// show the result
			if(currentStepValue==1){
			System.out.println(" and " + timesThorupVisit[currentStep] +
					" ms for visiting all vertices (average of " +
					numberOfPasses + " passes).");
			
			}
			// prepare next iteration
			currentStepValue += allStepValues[0];
			currentStep++;
			//sourceID++;
			customUpdate();
		}
		
		if (writeTableColumns) {
			// write all values and running times to the console
			System.out.println();
			
			System.out.println("All step values:");
			//writeTableColumn(allStepValues);
			//writeLatexTableRow(allStepValues);
			System.out.println("Loop: "+maximumStepValue+"times.");
			System.out.println();
			
			System.out.println("Avarage times of Dijkstra with array heap:");
			//writeTableColumn(timesDijkstraArrayHeap);
			//writeLatexTableRow(timesDijkstraArrayHeap);
			System.out.println(new java.text.DecimalFormat("0.0000").format(getHalfAverage(timesDijkstraArrayHeap)*1000000));
			System.out.println();
			
			System.out.println("Avarage times of Dijkstra with Fibonacci heap:");
			//writeTableColumn(timesDijkstraFibHeap);
		//	writeLatexTableRow(timesDijkstraFibHeap);
			System.out.println(getHalfAverage(timesDijkstraFibHeap)*1000000);
			System.out.println();
			
			System.out.println("Average times of Thorup (Total)");
			//writeTableColumn(timesThorupVisit);
			//writeLatexTableRow(timesThorupVisit);
			System.out.println(new java.text.DecimalFormat("0.0000").format((getHalfAverage(timesThorupVisit)+getAverage(timesThorupMST)+getAverage(timesThorupDS))*1000000));
			
			System.out.println();
			
			/*
			System.out.println("Average times of Thorup (construct MST):");
		//	writeTableColumn(timesThorupMST);
			System.out.println(getAverage(timesThorupMST)*1000000);
		//	writeLatexTableRow(timesThorupMST);
			System.out.println();
			
			System.out.println("Average times of Thorup (construct other DS):");
			//writeTableColumn(timesThorupDS);
			//writeLatexTableRow(timesThorupDS);
			System.out.println(getAverage(timesThorupDS)*1000000);
			System.out.println();
			*/
			
			System.out.println("Average times of Thorup (construct structures):");
			//	writeTableColumn(timesThorupMST);
				System.out.println(new java.text.DecimalFormat("0.0000").format((getHalfAverage(timesThorupMST)+getAverage(timesThorupDS))*1000000));
			//	writeLatexTableRow(timesThorupMST);
				System.out.println();
				
			
			System.out.println("Average times of Thorup (visit)");
			//writeTableColumn(timesThorupVisit);
			//writeLatexTableRow(timesThorupVisit);
			System.out.println(new java.text.DecimalFormat("0.0000").format(getHalfAverage(timesThorupVisit)*1000000));
			System.out.println();
		}
	}
	
	/**
	 * Generates the connected, weighted, undirected graph for the next
	 * iteration of this series of measurement.
	 * 
	 * @return
	 * 		the connected, weighted, undirected graph for the next iteration of
	 * 		this series of measurement
	 */
	public abstract AdjacencyListWeightedDirectedGraph<WeightedEdge>
		generateGraph();
	
	/**
	 * Reads some values from user input, depending on the type of this series
	 * of measurement.
	 */
	public abstract void readValuesFromUserCustom();
	
	/**
	 * Allows subclasses to update some values in every iteration, depending on
	 * the type series of measurement.
	 */
	public abstract void customUpdate();
	
	/**
	 * Reads the parameters for this series of measurement from user input.
	 */
	public void readValuesFromUser() {
		System.out.println("This application allows measuring the " +
				"performance of Thorup's SSSP algorithm compared to " +
				"the one of the algorithm by Dijkstra.");
		
		System.out.println();
		
		in = new Scanner(System.in);
		
		System.out.print("Please specify the number of passes for each " +
				"experiment to take the average of: ");
		numberOfPasses = in.nextInt();
		//numberOfPasses = 100;
		System.out.print("Please specify the maximum number of vertices: ");
		numberOfVertices = in.nextInt();
		numberOfVerticesCurrent = numberOfVertices;
		
		System.out.print("Please specify the maximum number of edges per" +
				" vertex: ");
		numberOfEdgesPerVertex = in.nextInt();
		
		System.out.print("Please specify the the maximum edge weight: ");
		maximumEdgeWeight = in.nextInt();

		readValuesFromUserCustom();
		
		System.out.print("Do you want to write the table columns to the " +
				"console (Y/N)? ");
		in.nextLine();
		writeTableColumns = (in.nextLine()).equals("Y");
		writeTableColumns = true;
	}
	
	/**
	 * Computes and returns the average of the values in the passed
	 * <code>long</code> array.
	 * 
	 * @param values
	 * 		the values to compute the average of
	 * @return
	 * 		the average of the values in the passed <code>long</code> array
	 */
	public double getAverage(double[] values) {
		double avg = values[0];
		
		for (int i = 1; i < values.length; i++) {
			avg += values[i];
		}
		avg /= values.length;
		avg /=1000000;
		return avg;
	}
	
	public double getHalfAverage(double[] values) {
		double avg = 0;
		int count=0;
		for (int i = 0; i < values.length; i++) {
			if(i>=values.length/2){
				avg += values[i];
				count++;
			}
		}
		avg /= count;
		avg /=1000000;
		
		return avg;
	}
	/**
	 * Writes a column with the values of the passed <code>int</code> array
	 * to the console.
	 * 
	 * @param values
	 * 		the values to write to the console
	 */
	public void writeTableColumn(double[] values) {
		double total=0;
		for (int i = 0; i < values.length; i++) {
			total+=values[i];
		}
		System.out.println(df.format(total));
	}
	
	/**
	 * Writes a column with the values of the passed <code>long</code> array
	 * to the console.
	 * 
	 * @param values
	 * 		the values to write to the console
	 */
	public void writeTableColumn(long[] values) {
		for (int i = 0; i < values.length; i++) {
			System.out.println(values[i]);
		}
	}
	
	/**
	 * Writes a latex table row with every second value of the passed
	 * <code>int</code> array to the console.
	 * 
	 * @param values
	 * 		the values to write to the console
	 */
	public void writeLatexTableRow(double[] values) {
		for (int i = 1; i < values.length; i+= 2) {
			System.out.print(values[i]);
			
			if (i < values.length - 1) {
				System.out.print(" & ");
			} else {
				System.out.println("\\\\");
			}
		}
	}
	
	/**
	 * Writes a latex table row with every second value of the passed
	 * <code>long</code> array to the console.
	 * 
	 * @param values
	 * 		the values to write to the console
	 */
	public void writeLatexTableRow(long[] values) {
		for (int i = 1; i < values.length; i+= 2) {
			System.out.print(values[i]);
			
			if (i < values.length - 1) {
				System.out.print(" & ");
			} else {
				System.out.println("\\\\");
			}
		}
	}
}
