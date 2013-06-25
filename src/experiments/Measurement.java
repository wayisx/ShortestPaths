package experiments;

import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.ListIterator;
import java.util.Scanner;




import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.io.ParseException;
import com.vividsolutions.jts.io.WKTReader;

import de.unikiel.npr.thorup.algs.Dijkstra;
import de.unikiel.npr.thorup.algs.Kruskal;
import de.unikiel.npr.thorup.algs.NewKColorDijkstra;
import de.unikiel.npr.thorup.algs.Thorup;
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
public abstract class Measurement {
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
	
	double[] timesNewKColorDijkstra;
	
	/**
	 * The times required for constructing the <i>msb</i>-minimum spanning trees
	 * for <i>Thorup</i>'s algorithm in the most recent performance tests.
	 */
	double[] timesThorupMST;
	double[] timesThorup2MST;
	
	/**
	 * The times required for constructing the component tree, initializing the
	 * bucket structure, and preparing the unvisited data structure for
	 * <i>Thorup</i>'s algorithm in the most recent performance tests.
	 */
	double[] timesThorupDS;
	double[] timesThorup2DS;
	
	/**
	 * The running times of the visiting part of <i>Thorup</i>'s algorithm in
	 * the most recent performance tests.
	 */
	double[] timesThorupVisit;
	double[] timesThorup2Visit;
	
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
		FileReader fRdr = new FileReader("D:\\Down\\GIS Data\\Shapefiles\\gm-jpn-trans_u_2\\1\\roadl_jpn.WKT");    
		
		WKTReader wktRdr = new WKTReader();
		Geometry geom = wktRdr.read(fRdr);
		
		
		MultiLineString mls=(MultiLineString)geom;
		Hashtable hb=new Hashtable(mls.getNumGeometries()*2);
		Hashtable hbw=new Hashtable(mls.getNumGeometries()*2);
		
		int vertexNum=0;
		List edges= new ArrayList<WeightedEdge>();
		List lengthDist = new ArrayList();
		int sum =0;
		   	for (int i = 0; i < mls.getNumGeometries(); i++) {

	   //	for (int i = 0; i < 1000; i++) {
		   		 LineString ls = (LineString) mls.getGeometryN(i);
		   		 String sp = ls.getStartPoint().toString();
		   		 String ep = ls.getEndPoint().toString();

		   		// int weight =  (int) (ls.getLength()*1000000);
		   		
		   		 /*Calculate distance in meters*/
		   		 int weight = (int)(length(ls.getCoordinateSequence()));
		   		// sum = sum+ weight;
		   		 if(weight<1){
		   			weight = 1;
		   			// System.out.println(ls);
		   		 }
		   
		   		//System.out.println(sums);
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
		   		 String rwe = eEdge+"-"+sEdge;
		   		if(hbw.get(we)==null&&hbw.get(rwe)==null){
		   		 hbw.put(we, weight);
		   		 edges.add(new WeightedEdge(sEdge, eEdge));
		   		}else if(hbw.get(we)==null){ //There are more than one way between two vertexes.
		   			int w = Integer.parseInt(hbw.get(rwe).toString());
		   			if(w>weight){
		   				hbw.remove(rwe);
		   				hbw.put(rwe, weight);
		   			}
		   		}else{
		   			int w = Integer.parseInt(hbw.get(we).toString());
		   			if(w>weight){
		   				hbw.remove(we);
		   				hbw.put(we, weight);
		   			}
		   		}
			}
		   	fRdr.close();		
		   	AdjacencyListWeightedDirectedGraph<WeightedEdge> graph =
				new AdjacencyListWeightedDirectedGraph<WeightedEdge>(vertexNum);
		   	for(int i = 0; i<edges.size(); i++){
		   		WeightedEdge we = (WeightedEdge) edges.get(i);

		   		 int sEdge=we.getSource();
		   		 int eEdge=we.getTarget();
		   		 int weight = Integer.parseInt(hbw.get(new String(sEdge+"-"+eEdge)).toString());
		   		graph.addEdge(new WeightedEdge(sEdge,eEdge,weight));
		   		graph.addEdge(new WeightedEdge(eEdge,sEdge,weight));
		   		
		   	//	weight =1;
		   		/* Add distinct length to list*/
		   		if(!lengthDist.contains(weight)){
		   			lengthDist.add(weight);
		   		}
		   	}
		   	
		    /*Sort list*/
		    //Collections.sort(lengthDist);
		   	/*
		   	lengthDist.add(3);
			lengthDist.add(2);
			lengthDist.add(5);
			lengthDist.add(8);
			lengthDist.add(4);
			lengthDist.add(9);
			lengthDist.add(1);
			lengthDist.add(10);
			lengthDist.add(18);*/
		    Object[] a = lengthDist.toArray();
		  	quicksort(a, 0, a.length-1);
		  	ListIterator i = lengthDist.listIterator();
		  	for (int j=0; j<a.length; j++) {
		  	    i.next();
		  	    i.set(a[j]);
		  	}
		   	
		/*
		
		   	AdjacencyListWeightedDirectedGraph<WeightedEdge> graph =
				generateGraph();
		   int	vertexNum=10000;*/
		System.out.println(" Graph has " +graph.getNumberOfEdges() + " edges,"+ vertexNum+" vertexes.");
		System.out.println("");
		numberOfPasses =1;
		currentStepValue = 1;
		maximumStepValue =50;
		numberOfVerticesCurrent = vertexNum;
		// initialize all variables
		numberOfSteps = maximumStepValue / currentStepValue;
		
		timesToComputeTheAverageOf = new double[numberOfPasses];
		
		timesDijkstraArrayHeap = new double[numberOfSteps];
		timesDijkstraFibHeap = new double[numberOfSteps];
		timesThorupMST = new double[numberOfSteps];
		timesThorupDS = new double[numberOfSteps];
		timesThorupVisit = new double[numberOfSteps];
		timesThorup2MST = new double[numberOfSteps];
		timesThorup2DS = new double[numberOfSteps];
		timesThorup2Visit = new double[numberOfSteps];
		timesNewKColorDijkstra = new double[numberOfSteps];
		
		
		
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
			
			
			
			
			
			/**
			 * Preparing for New K color Dijkstra
			 * 
			 * */
			
			NewKColorDijkstra NewKColorDijkstra = new NewKColorDijkstra();
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				NewKColorDijkstra.initialize(graph, sourceID, lengthDist,
						new ArrayPriorityQueue<Integer>
							(numberOfVerticesCurrent));
				stop = System.nanoTime();
				/*
				int pre[] = NewKColorDijkstra.getDistances();
				for(int v = 0; v<vertexNum; v++){
					System.out.println(pre[v]);
					
				}*/
				timesToComputeTheAverageOf[pass] = stop - start;
			}
			timesNewKColorDijkstra[currentStep] = getAverage(timesToComputeTheAverageOf);
			

			
			//Preparing for Thorup
			Thorup thorup = new Thorup();
			
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
			//System.out.print(" Took " + timesThorupMST[currentStep] + " ms for constructing the MST, ");
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				thorup.constructOtherDataStructures
					(new UnionFindStructureTarjan<Integer>(),
					 new SplitFindminStructureGabow<Integer>
						(numberOfVerticesCurrent));
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}

			// compute the average
			timesThorupDS[currentStep] =
				getAverage(timesToComputeTheAverageOf);
					
			// show the result
			//System.out.print(timesThorupDS[currentStep] +" ms for constructing the other data structures,");
			
		//	System.out.println();
			
			
			Thorup2 thorup2 = new Thorup2();
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				thorup2.constructMinimumSpanningTree(graph,
						new Kruskal(new UnionFindStructureTarjan<Integer>()));
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}

			// compute the average
			timesThorup2MST[currentStep] =
				getAverage(timesToComputeTheAverageOf);
					
			// show the result
		//	System.out.print(" Took " + timesThorupMST[currentStep] + " ms for constructing the MST, ");
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				thorup2.constructOtherDataStructures
					(new UnionFindStructureTarjan<Integer>());
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}

			// compute the average
			timesThorup2DS[currentStep] =
				getAverage(timesToComputeTheAverageOf);
			
			
			
			// run Dijkstra's algorithm with an array heap and take the time
			if(currentStepValue==1){
				System.out.print("Running Dijkstra with an array priority " +
				"queue...");
			}
			
			sourceID=0;
			Dijkstra dijkstra = new Dijkstra();
			
			for (int pass = 0; pass < numberOfPasses; pass++) {
				start = System.nanoTime();
				dijkstra.findShortestPaths(graph, sourceID,
						new ArrayPriorityQueue<Integer>
							(numberOfVerticesCurrent));
				stop = System.nanoTime();
			/*	
				int pre[] = dijkstra.getPredecessors();
				for(int v = 0; v<vertexNum; v++){
					System.out.println(pre[v]);
					
				}*/
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
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}
	
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
			for (int pass = 0; pass < numberOfPasses; pass++) {
				if (currentStep > 0) {
					thorup.cleanUpBetweenQueries
						(new SplitFindminStructureGabow<Integer>
							(graph.getNumberOfVertices()));
				}

				start = System.nanoTime();
				thorup.findShortestPaths(sourceID);
				stop = System.nanoTime();
			
				timesToComputeTheAverageOf[pass] = stop - start;
			}
	
			// compute the average
			timesThorupVisit[currentStep] =
				getAverage(timesToComputeTheAverageOf);
		
					
			// show the result
			if(currentStepValue==1){
			System.out.println(" and " + timesThorupVisit[currentStep] +
					" ms for visiting all vertices (average of " +
					numberOfPasses + " passes).");
			
			}
			


			
			
//==================
			int[] td = new int[numberOfVerticesCurrent];
			for (int pass = 0; pass < numberOfPasses; pass++) {
				
				if (currentStep > 1) {
					thorup2.cleanUpBetweenQueries();
				}

				start = System.nanoTime();
				 td = thorup2.findShortestPaths(sourceID);
				stop = System.nanoTime();
				
				timesToComputeTheAverageOf[pass] = stop - start;
			}
			
			/*
			for(int i = 0; i<td.length;i++){
				
				System.out.println(i+"-> "+td[i]);
			}*/
			// compute the average
			timesThorup2Visit[currentStep] =
				getAverage(timesToComputeTheAverageOf);
		
			
		
			
			
			
			
			
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
			
			
			
			System.out.println("Avarage times of New K Color Dijkstra with array heap:");
			//writeTableColumn(timesDijkstraArrayHeap);
			//writeLatexTableRow(timesDijkstraArrayHeap);
			System.out.println(new java.text.DecimalFormat("0.0000").format(getHalfAverage(timesNewKColorDijkstra)*1000000));
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
			
			/*
			System.out.println("Average times of Thorup (Total)");
			//writeTableColumn(timesThorupVisit);
			//writeLatexTableRow(timesThorupVisit);
			System.out.println(new java.text.DecimalFormat("0.0000").format((getHalfAverage(timesThorupVisit)+getAverage(timesThorupMST)+getAverage(timesThorupDS))*1000000));
			
			System.out.println();*/
			
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
			
			
			/*
			System.out.println("Average times of Thorup2 (Total)");
			//writeTableColumn(timesThorupVisit);
			//writeLatexTableRow(timesThorupVisit);
			System.out.println(new java.text.DecimalFormat("0.0000").format((getHalfAverage(timesThorup2Visit)+getAverage(timesThorup2MST)+getAverage(timesThorup2DS))*1000000));
			
			System.out.println();
			*/
			
			System.out.println("Average times of Thorup2 (construct structures):");
			//	writeTableColumn(timesThorupMST);
				System.out.println(new java.text.DecimalFormat("0.0000").format((getHalfAverage(timesThorup2MST)+getAverage(timesThorup2DS))*1000000));
			//	writeLatexTableRow(timesThorupMST);
				System.out.println();
				
			
			System.out.println("Average times of Thorup2 (visit)");
			//writeTableColumn(timesThorupVisit);
			//writeLatexTableRow(timesThorupVisit);
			System.out.println(new java.text.DecimalFormat("0.0000").format(getHalfAverage(timesThorup2Visit)*1000000));
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
	
	/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	/*::                                                                         :*/
	/*::  This routine calculates the distance between two points (given the     :*/
	/*::  latitude/longitude of those points). It is being used to calculate     :*/
	/*::  the distance between two locations using GeoDataSource (TM) prodducts  :*/
	/*::                                                                         :*/
	/*::  Definitions:                                                           :*/
	/*::    South latitudes are negative, east longitudes are positive           :*/
	/*::                                                                         :*/
	/*::  Passed to function:                                                    :*/
	/*::    lat1, lon1 = Latitude and Longitude of point 1 (in decimal degrees)  :*/
	/*::    lat2, lon2 = Latitude and Longitude of point 2 (in decimal degrees)  :*/
	/*::    unit = the unit you desire for results                               :*/
	/*::           where: 'M' is statute miles                                   :*/
	/*::                  'K' is kilometers (default)                            :*/
	/*::                  'N' is nautical miles                                  :*/
	/*::  Worldwide cities and other features databases with latitude longitude  :*/
	/*::  are available at http://www.geodatasource.com                          :*/
	/*::                                                                         :*/
	/*::  For enquiries, please contact sales@geodatasource.com                  :*/
	/*::                                                                         :*/
	/*::  Official Web site: http://www.geodatasource.com                        :*/
	/*::                                                                         :*/
	/*::           GeoDataSource.com (C) All Rights Reserved 2013                :*/
	/*::                                                                         :*/
	/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	private double distance(double lat1, double lon1, double lat2, double lon2, String unit) {
	  double theta = lon1 - lon2;
	  double dist = Math.sin(deg2rad(lat1)) * Math.sin(deg2rad(lat2)) + Math.cos(deg2rad(lat1)) * Math.cos(deg2rad(lat2)) * Math.cos(deg2rad(theta));
	  dist = Math.acos(dist);
	  dist = rad2deg(dist);
	  dist = dist * 60 * 1.1515;
	  if (unit == "K") {
	    dist = dist * 1.609344;
	  } else if (unit == "N") {
	  	dist = dist * 0.8684;
	    }
	  return (dist);
	}

	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	/*::  This function converts decimal degrees to radians             :*/
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	private double deg2rad(double deg) {
	  return (deg * Math.PI / 180.0);
	}

	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	/*::  This function converts radians to decimal degrees             :*/
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	private double rad2deg(double rad) {
	  return (rad * 180 / Math.PI);
	}
	
  public  double length(CoordinateSequence pts) {
      if (pts.size() < 1) return 0.0;
      double sum = 0.0;
      for (int i = 1; i < pts.size(); i++) {
        sum += distance(pts.getCoordinate(i).y, pts.getCoordinate(i).x,pts.getCoordinate(i - 1).y, pts.getCoordinate(i - 1).x, "K");
      }
      return sum;
  }
  
  
  
  public static void quicksort (Object[] a, int s, int e){
	  if(s<e){
		  int i = s-1;
		  int key = (Integer) a[e];
		  int temp = 0;
		  for (int j=s;j<e;j++){
			  if((Integer) a[j]<key){
				  i++;
				  temp = (Integer) a[i];
				  a[i] = a[j];
				  a[j] = temp;
			  }
		  }
		  temp = (Integer) a[i+1];
		  a[i+1] = key;
		  a[e] = temp;
		  quicksort(a, s, i);
		  quicksort(a, i+2, e);
	  }
  	
  }
  
}
