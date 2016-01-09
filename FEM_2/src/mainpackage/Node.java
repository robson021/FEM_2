package mainpackage;

public class Node {
	private final double TEMP_BEGIN;
	private final int R_COORDINATE;	
	public static final int CONVECTION_CONDITION = 1;
	/*
	 * Boundary condition types:	 
	 * 0 = none
	 * 1 = convection
	 */
	private final int BC; 
	
	public Node(int r, double t, int bc) {
		TEMP_BEGIN = t; R_COORDINATE = r; BC = bc;
	}

	public double getTEMP_BEGIN() {
		return TEMP_BEGIN;
	}
	public int getR_COORDINATE() {
		return R_COORDINATE;
	}
	public int getBC() {
		return BC;
	}
}
