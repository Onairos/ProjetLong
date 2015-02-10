package com.prieul.cluster.filein.generator;

import java.io.FileWriter;
import java.io.IOException;

import com.prieul.cluster.geometry.Point;

public class GeometricFile2DGenerator {
	
	public static int nbPoints = 0;
	public static final int dim = 2;
	public static String buffer = "";
	public static final String cheminFichier = "gen2D/";
	public static final String nomFichier = "fichier.txt";
	
	public static void main(String[] args) {
		
		//Création d'un premier rectangle "plein"
		Point A = new Point(-1, 1);
		Point B = new Point(1, -1);
		writeRectangle(A, B, 0.1);
		
		//Création d'un deuxième rectangle "plein"
		Point C = new Point(3, 1);
		Point D = new Point(4, -1);
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
	
	public static void writeRectangle(Point supGauche, Point infDroite, double pas) {
		double debutX = supGauche.getX();
		double finX = infDroite.getX();
		double debutY = infDroite.getY();
		double finY = supGauche.getY();
		for (double i=debutX; i<finX; i+=pas) {
			for (double j=debutY; j<finY; j+=pas) {
				buffer += i + " " + j + "\n";
				nbPoints++;
			}
		}
		//Ecriture du dernier point
		buffer += finX + " " + finY + "\n";
		nbPoints++;
	}
	

}
