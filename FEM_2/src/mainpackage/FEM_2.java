package mainpackage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import Jama.Matrix;


/*
 * Jama download: http://math.nist.gov/javanumerics/jama/Jama-1.0.3.jar
 */
public class FEM_2 {	
		
	// for .csv files
	private static final String SEPARATOR = ",", FILE_NAME = "result.txt";
	
	private List<Node> nodeList;
	private List<Element> elementList;
			
	private double[][] fMatrix, globalMatrix;
	
	// constructor
	public FEM_2() {
		nodeList = new ArrayList<>();
		elementList = new ArrayList<>();
	}
	
	public void runProgram() {		
		String fileName;
		fileName = "src/data.txt";
		
		loadDataFromFile(fileName);		
		printBasicData();
		
		Element.initFmatrix(nodeList.size());
		for (Element e : elementList)
			e.initLocalMatrix();
		
		initGlobalMatrix();
		printMatrixes();	
		computeTemperatures();
//		try {
//			//printAndSaveToFile(solved);
//		} catch (IOException e1) {
//			e1.printStackTrace();
//			System.out.println("Error. Cold not save result to the file.");
//		}
	}
	
	private void loadDataFromFile(final String FILE_NAME) {
		BufferedReader fr = null; // file reader
		try {
			fr = new BufferedReader(new FileReader(FILE_NAME));
			String line = null;
			
			// set up basic input
			Element.setAlpha(Integer.parseInt(fr.readLine()));
			Element.setTemperatureOfEnnv(Integer.parseInt(fr.readLine()));
			Element.setTime(Integer.parseInt(fr.readLine()));
			
			Node n1, n2;
			Element elem;
			double ev[] = new double[3];
			String[] val = fr.readLine().split(SEPARATOR);
			if (val.length != 3) {
				System.out.println("Bad input file.");
				System.exit(1);
			}
			n1 = new Node(Double.parseDouble(val[0]), Double.parseDouble(val[1]), Integer.parseInt(val[2]));
			nodeList.add(n1);
			while ((line = fr.readLine()) != null)
			{
				// loading element parameters
				val = line.split(SEPARATOR);			
				if (val.length != 3) break; // bad input line size
				for (int i=0;i<3;i++)
					ev[i] = Double.parseDouble(val[i]);
								
				// loading 2nd node parameters
				val = fr.readLine().split(SEPARATOR);		
				if (val.length != 3) break;
				n2 = new Node(Double.parseDouble(val[0]), Double.parseDouble(val[1]), Integer.parseInt(val[2]));
				nodeList.add(n2);
				
				elem = new Element(n1, n2, ev[0], ev[1], ev[2]);
				elementList.add(elem);
				
				n1 = n2;
			}
			fMatrix = new double[1][nodeList.size()];
			if(fr != null) fr.close();
		} catch (IOException e) {
			e.printStackTrace();			
			System.out.println("Error. Could not load data from file.");
			System.exit(1);
		} 
	}
	
	private void printBasicData() {
		System.out.printf("Alpha: %d\nTemp. of environment: %d\nTime: %d\n",
				Element.getAlpha(), Element.getTempOfEnv(), Element.getTime());
		
		System.out.println("Nodes coordinates and temperatures:");
		for (Node n : nodeList) {
			System.out.printf("r: %.2f; t: %.2f\n", n.getR_COORDINATE(), n.getTEMP_BEGIN());
		}
		
		System.out.println("\nElements values:");
		for (Element e : elementList) {
			System.out.printf("ro: %.2f; c: %.2f, k: %.2f; length: %.2f\n", e.getRO(), e.getC(), e.getK(), e.getL());
		}
		System.out.print("\tElement thickness (max R): " + Element.getTotalWidth());
		System.out.printf("\n\tNodes total: %d; Elements total: %d\n", nodeList.size(), elementList.size());
	}
	
	private void printMatrixes() {
		final int SIZE = globalMatrix.length;
		System.out.printf("\nGlobal matrix [%d x %d]:\n", SIZE, SIZE);
		for (int j,i=0; i<SIZE; i++) {
			for (j=0; j<SIZE; j++)
				System.out.printf("%.2f ", globalMatrix[i][j]);
			System.out.println("");
		}
		/*System.out.println("\nF matrix:");
		double[][] fMatrix = Element.getFmatrix();
		for (int i=0;i<SIZE;i++) {
			System.out.printf("%.2f ",fMatrix[0][i]);
		} System.out.println("");*/
	}
	
	private void initGlobalMatrix() {
		final int SIZE = nodeList.size(); // = number of elements + 1
		globalMatrix = new double[SIZE][SIZE];
		
		int i, j; i=j=0;
		for (Element e : elementList) {
			double[][] matrix = e.getLocalMatrix();
			
			globalMatrix[i][j] += matrix[0][0];
			globalMatrix[i][j+1] += matrix[0][1];
			globalMatrix[j+1][i] += matrix[1][0];
			globalMatrix[j+1][i+1] += matrix[1][1];
			
			i++; j++;
		}		
	}
	
	private void printAndSaveToFile(double[][] solved) throws IOException {
		FileWriter fw = null; // file writer
		fw = new FileWriter(new File(FILE_NAME));
		fw.append(Calendar.getInstance().getTime().toString() + "\n\n");
		System.out.println("\n\tTemperatures:");
		for (int i=0; i<nodeList.size(); i++) {
			System.out.printf("%.2f ",solved[i][0]);
			fw.append(String.valueOf(solved[i][0]));
			fw.append("\n");
		}
		if (fw != null) fw.close();
	}
	
	private void computeTemperatures() {
		Matrix H = new Matrix(globalMatrix);
		//return H.solve(new Matrix(Element.getFmatrix()).transpose().uminus()).getArray();		
		//return H.solve(new Matrix(Element.getFmatrix()).transpose()).getArray();	
		
		if (nodeList.size() != (elementList.size()+1))
			System.exit(1);
		
		double avgDtau = .0;
		for (Element e : elementList)
			avgDtau += e.get_dTau();
		avgDtau = avgDtau / elementList.size();
		
		int nTime = (int) ((Element.getTime() / avgDtau) +1);
		System.out.printf("Iterations: %d\n", nTime);
		final int SIZE = elementList.size() +1;
		double[][] temperatureVector = new double[1][SIZE];
		int i=0;
		for (Node n : nodeList)
			temperatureVector[0][i++] = n.getTEMP_BEGIN();
		
		System.out.println("Nodes: "+i);
		Matrix temperatureMatrix = new Matrix(temperatureVector).transpose();
		Matrix fMatrix = new Matrix(Element.getFmatrix());
		fMatrix = fMatrix.transpose();
		
		for (i=0; i<nTime; i++) {			
			Matrix newH = H.times(temperatureMatrix);
			Matrix newTempV = newH.solve(fMatrix);
			
			
			temperatureMatrix = newTempV;
		}
		
	}
	

	public static void main(String[] args) {
		new FEM_2().runProgram();
		System.out.print("\n\nEND");
		System.exit(0);
	}

}
