package mainpackage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Scanner;

import Jama.Matrix;


/*
 * Jama download: http://math.nist.gov/javanumerics/jama/Jama-1.0.3.jar
 */
public class FEM_2 {	
		
	// for .csv files
	private static final String SEPARATOR = ",", FILE_NAME = "result.txt";
	
	private List<Node> nodeList;
	private List<Element> elementList;
			
	private double[][] globalMatrix;
	
	// constructor
	public FEM_2() {
		nodeList = new ArrayList<>();
		elementList = new ArrayList<>();
	}
	
	public void runProgram() {		
		String fileName = null;
		
		System.out.println("Enter file name, example: src/example.txt");
		Scanner sc = new Scanner(System.in);
		fileName = sc.nextLine();
		sc.close();
		
		//fileName = "src/example.txt";		
		loadDataFromFile(fileName);		
		printBasicData();
		
		for (Element e : elementList)
			e.initLocalMatrix();
		
		initGlobalMatrix();
		printMatrixes();	
		computeTemperatures();
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
			System.out.printf("r: %.2f; t: %.2f\n", n.getR_COORDINATE(), n.getTemp());
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
	
	private void computeTemperatures() {
		Matrix H = new Matrix(globalMatrix);	
				
		double avgDtau = .0;
		for (Element e : elementList)
			avgDtau += e.get_dTau();
		avgDtau = avgDtau / elementList.size();
		
		int nTime = (int) ((Element.getTime() / avgDtau) +1);
		System.out.printf("\n\tIterations: %d\n", nTime);
		final int SIZE = elementList.size() +1;
		
		//System.out.println("avg_dTau * nTime = " + (nTime * avgDtau));
		
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<nTime; i++)
		{			
			System.out.printf("\nStep: %d; time: %.2f", (i+1), (avgDtau * i));
			for (Element e : elementList)
				e.initLocalFmatrix();
			
			double[][] fVector = initNewF_vector(SIZE);
			
			double[][] tempVector = H.solve(new Matrix(fVector).transpose()).getArray();
			
			// print and append to the string
			sb.append("Step: "+ (i+1) + "; time: " + (int)(i*avgDtau) + "\n");
			System.out.print("\nTemperatures: ");
			for (int j=0; j<SIZE; j++) {
				System.out.printf("%.2f ", tempVector[j][0]);
				sb.append(String.format("%.2f", tempVector[j][0]) + " ");
				//sb.append(new DecimalFormat("#.##").format(tempVector[j][0]) + " ");
			} System.out.println(""); sb.append("\n");			
			updateNodesTemperature(tempVector);			
		}
		try {
			saveToFile(sb.toString());
		} catch (IOException e1) {
			e1.printStackTrace();
			System.out.println("\nError. Could not save results to the file.");
		}
	}
	
	private double[][] initNewF_vector(int size) {
		double[][] fVector = new double[1][size];
		
		int i=0;
		for (Element e : elementList)
		{
			double[][] localM = e.getLocalFmatrix();
			fVector[0][i] += localM[0][0];
			fVector[0][i+1] += localM[0][1];
			i++;
		}
		/*System.out.print("Vector F: ");
		for (i=0;i<size;i++)
			System.out.printf("%.2f", fVector[0][i]);*/
		
		return fVector;
	}
	
	private void updateNodesTemperature(double[][] tv) {
		int i=0;
		for (Node n : nodeList)
			n.setTemp(tv[i++][0]);
	}
	
	private void saveToFile(String text) throws IOException {
		FileWriter fw = null;
		try {
			fw = new FileWriter(new File(FILE_NAME));
			fw.append(Calendar.getInstance().getTime().toString() + "\n\n");
			fw.append(text);
			System.out.println("\n\tSaved to: "+FILE_NAME);
		} finally {
			if (fw != null) fw.close();
		}
	}

	public static void main(String[] args) {
		new FEM_2().runProgram();
		System.out.print("\n\nEND");
		System.exit(0);
	}

}
