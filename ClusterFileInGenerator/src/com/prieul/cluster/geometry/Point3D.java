package com.prieul.cluster.geometry;

public class Point3D extends Point {
	
	private double z;
	
	public Point3D(double x, double y, double z) {
		super(x, y);
		this.setZ(z);
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

}
