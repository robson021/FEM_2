package mainpackage;

public class Element {	
	// steel density (ro) = 7,86 g/cm3	
	
	// points for Gauss numerical integration
	private static final double E_2 = 0.5773502692, E_1 = -E_2;
	private static final int W = 1; // weight of integration
	private static final double SIZE = 2.0;
	
	private final Node NODE_1, NODE_2;
	private final double RO, C, K;
	private final double L; // length of element
	private static int time;
	private static int dTime;
	private static double totalWidth = 0; // max r 	
	private static int alpha;
	private static int tempOfEnvironment;
	private double rMax;
	private static int counter = 0; // for F_matrix creation. Node counter
	private double[] rp;
	
	private double dTau;
	
	// possible shape functions
	private static final double[] N = {(1-E_1)/SIZE, (1-E_2)/SIZE, (1+E_1)/SIZE, (1+E_2)/SIZE};
	
	private double[][] localMatrix;
	private static double[][] fMatrix;
	
	public Element(Node n1, Node n2, double ro, double c, double k) {
		NODE_1 = n1; NODE_2 = n2;
		RO = ro; C = c; K = k;
		L = NODE_2.getR_COORDINATE() - NODE_1.getR_COORDINATE();
		totalWidth += L;				
		rMax = NODE_2.getR_COORDINATE();
		double a = k/(c*ro);
		dTau = L*L/(0.5*a);
		double nTime = (time/dTau) +1;
		dTau = time/nTime;
		//dTau = time;
	}
	
	public static void initFmatrix(int size) {
		fMatrix = new double[1][size];
	}
	
	public void initLocalMatrix() {
		localMatrix = new double[2][2];
		rp = new double[2];
		double Rp;
		for (int x=0;x<2;x++) {
			Rp = Ni(x) * NODE_1.getR_COORDINATE() + Nj(x) * NODE_2.getR_COORDINATE();
			rp[x] = Rp;
			double tmp = K * Rp * W/L;
			localMatrix[0][0] += tmp + C*RO*L*Rp*W * Ni(x) * Ni(x) / dTau;
			localMatrix[0][1] += -tmp + C*RO*L*Rp*W * Ni(x)*Nj(x) / dTau;			
			localMatrix[1][1] += K*Rp*W/L + C*RO*Rp*L*W * Nj(x)*Nj(x) / dTau;
			
			// bc check
			if (NODE_2.getBC() == Node.CONVECTION_CONDITION)
				localMatrix[1][1] += 2*alpha*rMax;	
			
			// fmatrix
			fMatrix[0][counter] += C*RO*L* (Ni(x)*NODE_1.getTEMP_BEGIN() + 
					Nj(x)*NODE_2.getTEMP_BEGIN()*Rp*W*Ni(x))/dTau;
			
			fMatrix[0][counter+1] += C*RO*L* (Ni(x)*NODE_1.getTEMP_BEGIN() + 
					Nj(x)*NODE_2.getTEMP_BEGIN()*Rp*W*Nj(x))/dTau;
			if (NODE_2.getBC() == Node.CONVECTION_CONDITION)
				fMatrix[0][counter+1] += 2*rMax*alpha*tempOfEnvironment; 
			
		}
		localMatrix[1][0] = localMatrix[0][1];
		counter++;
	}
	
	//public double[][] getFmatrix()

	private double Ni(int i) {
		if (i==0) return N[0];
		else return N[1];
	}
	private double Nj(int j) {
		if (j==0) return N[2];
		else return N[3];
	}
	
	public static double getTotalWidth() {
		return totalWidth;
	}
	
	public double getRO() {
		return RO;
	}
	public double getC() {
		return C;
	}
	public double getK() {
		return K;
	}
	public double getL() {
		return L;
	}
	public double[][] getLocalMatrix() {
		return localMatrix;
	}	
	public static void setTime(int t) {
		time = t;
		dTime = time / 30;
	}
	public static int getTime() {
		return time;
	}
	public static void setAlpha(int a) {
		alpha = a;
	}
	public static int getAlpha() {
		return alpha;
	}
	public static double[][] getFmatrix() {
		return fMatrix;
	}
	public static void setTemperatureOfEnnv(int t) {
		tempOfEnvironment = t;
	}
	public static int getTempOfEnv() {
		return tempOfEnvironment;
	}
	public double get_dTau() {
		return this.dTau;
	}
}
