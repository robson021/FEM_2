package mainpackage;

public class Node {
	private double temperature;
	private final double R_COORDINATE;	
	public static final int CONVECTION_CONDITION = 1;
	
	/*
	 * Boundary condition types:	 
	 * 0 = none
	 * 1 = convection
	 */
	private final int BC; 
	
	public Node(double r, double t, int bc) {
		temperature = t; R_COORDINATE = r; BC = bc;
	}

	public double getTemp() {
		return temperature;
	}
	public double getR_COORDINATE() {
		return R_COORDINATE;
	}
	public int getBC() {
		return BC;
	}
	public void setTemp(double t) {
		temperature = t;
	}
}
