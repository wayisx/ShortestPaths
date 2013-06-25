package experiments;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.MultiPoint;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.io.WKTReader;

import de.unikiel.npr.thorup.ds.graph.AdjacencyListWeightedDirectedGraph;
import de.unikiel.npr.thorup.ds.graph.WeightedEdge;
import de.unikiel.npr.thorup.util.RandomGraphGenerator;

/**
 * A series of measurent for running the algorithms by <i>Dijkstra</i>
 * and <i>Thorup</i> several times, taking the average over all passes, and
 * increasing the number of vertices step by step.
 * 
 * @author
 * 		<a href="mailto:npr@informatik.uni-kiel.de">Nick Pr&uuml;hs</a>
 * @version
 * 		1.0, 09/17/09
 */
public class Test extends Measurement {
	/**
	 * Generates the connected, weighted, undirected graph with the next number
	 * of vertices.
	 * 
	 * @return
	 * 		the connected, weighted, undirected graph with the next number
	 * 		of vertices
	 */
	public AdjacencyListWeightedDirectedGraph<WeightedEdge> generateGraph() {
		return RandomGraphGenerator.generateConnectedWeightedUndirectedGraph
			(currentStepValue, maximumEdgeWeight, numberOfEdgesPerVertex);
	}
	
	/**
	 * Reads the size of the steps the number vertices are increased in each
	 * iteration.
	 */
	public void readValuesFromUserCustom() {
		System.out.print("Please specify the steps the number of vertices " +
		"of the instances should be increased (step size): ");
		
		currentStepValue = in.nextInt();
		//currentStepValue = 1;
		maximumStepValue = numberOfVertices;
		
		numberOfVerticesCurrent = currentStepValue;
	}
	
	/**
	 * Increases the number of vertices for the next iteration of this series
	 * of measurement.
	 */
	public void customUpdate() {
		numberOfVerticesCurrent = currentStepValue;
	}
	
	
	/**
	 * Reads all required parameters from user input and runs a series of
	 * measurent for running the algorithms by <i>Dijkstra</i> and
	 * <i>Thorup</i> several times, taking the average over all passes, and
	 * increasing the number vertices step by step.
	 * 
	 * @param args
	 * 		ignored
	 */
	public static void main(String[] args) {
		Test m = new Test();
		
		//String j = null;
	//	System.out.println(j!=null);
		
		
		//m.readValuesFromUser();
		try {
			m.execute();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
