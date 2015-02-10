package com.prieul.cluster.filein.generator;

import java.io.FileWriter;
import java.io.IOException;

import com.prieul.cluster.geometry.Point3D;

public class GeometricFile3DGenerator {
	
	public static int nbPoints = 0;
	public static final int dim = 3;
	public static String buffer = "";
	public static final String cheminFichier = "gen3D/";
	public static final String nomFichier = "fichier.txt";

	
	public static void main(String[] args) {
		
		//Création d'un premier rectangle "plein"
		Point3D A = new Point3D(1, 1, 0);
		Point3D B = new Point3D(2, 0, 1);
		writeRectangle(A, B, 0.1);
		
		//Création d'un deuxième rectangle "plein"
		Point3D C = new Point3D(1, 1, 3);
		Point3D D = new Point3D(2, 0, 4);
		writeRectangle(C, D, 0.1);
		
		//Ecriture de la première ligne et des points
		String firstLine = nbPoints + " " + dim + "\n";
		buffer = firstLine + buffer;
		FileWriter fw = null;
		try {
			fw = new FileWriter(cheminFichier + nomFichier, false);
			fw.write(buffer, 0, buffer.length());
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (fw != null) {
				try {
					fw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		

		
	}
	
	/**
	 * Ecris un rectange à partir de deux points (simplifié)
	 * Attention les deux points spécifiés ne doivent pas avoir une coordonnée commune
	 * @param supGauche
	 * @param infDroite
	 * @param pas
	 */
	public static void writeRectangle(Point3D supGauche, Point3D infDroite, double pas) {
		//Explicitation des deux autres points
		//Point3D infGauche = new Point3D(supGauche.getX(), infDroite.getY(), supGauche.getZ());
		//Point3D supDroite = new Point3D(infDroite.getX(), supGauche.getY(), infDroite.getZ());
		//Débuts et fins de parcours des boucles
		double debutX = supGauche.getX();
		double finX = infDroite.getX();
		double debutY = infDroite.getY();
		double finY = supGauche.getY();
		double debutZ = supGauche.getZ() < infDroite.getZ() ? supGauche.getZ() : infDroite.getZ();
		double finZ = supGauche.getZ() >= infDroite.getZ() ? supGauche.getZ() : infDroite.getZ();
		for (double i=debutX; i<finX; i+=pas) {
			for (double j=debutY; j<finY; j+=pas) {
				for (double k=debutZ; k<finZ; k+=pas) {
					buffer += i + " " + j + " " + k + "\n";
					nbPoints++;
				}
			}
		}
		//Ecriture du dernier point
		buffer += finX + " " + finY + " " + finZ + "\n";
		nbPoints++;
	}

}
