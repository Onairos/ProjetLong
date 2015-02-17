package com.prieul.cluster.filein.generator;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

import com.prieul.cluster.geometry.Point3D;

public class GeometricFile3DGenerator {

	public static int nbPoints = 0;
	public static final int dim = 3;
	public static String buffer = "";
	public static final String cheminFichier = "gen3D/";
	public static final String nomFichier = "fichier.txt";


	public static void main(String[] args) {

		writeNuageSpherique( new Point3D(0,0,0), 0.01, 0.5, 1000);
		writeSurfaceSpherique( new Point3D(1,0,0), 0.5, 0.144);

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
		
		System.out.println("succes ecriture avec "+nbPoints+" points");

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

	/**
	 * Ecris une sphere à partir de la donnée d'un centre et d'une distance
	 * @param centre
	 * @param rayon
	 * @param pas
	 */
	public static void writeSphere(Point3D centre, double rayon, double pas) {
		double r2 = Math.sqrt(2);
		double dist = rayon;
		double debutX = centre.getX() - dist*r2;
		double debutY = centre.getY() - dist*r2;
		double debutZ = centre.getZ() - dist*r2;
		double finX = centre.getX() + dist*r2;
		double finY = centre.getY() + dist*r2;
		double finZ = centre.getZ() + dist*r2;
		System.out.println("");
		System.out.print("writeSphere ");
		for (double i=debutX; i<finX; i+=pas) {
			System.out.print(".");
			for (double j=debutY; j<finY; j+=pas) {
				for (double k=debutZ; k<finZ; k+=pas) {
					double distX = (i-centre.getX())*(i-centre.getX());
					double distY = (j-centre.getY())*(j-centre.getY());
					double distZ = (k-centre.getZ())*(k-centre.getZ());
					if( distX+distY+distZ < dist*dist){
						buffer += i + " " + j + " " + k + "\n";
						nbPoints++;
					}
				}
			}
		}
		System.out.println("");
		System.out.println("Sphere de centre ("+centre.getX()+", "+centre.getY()+", "+centre.getZ()+") de rayon "+rayon);
	}

	/**
	 * Ecris un nuage de points à partir de la donnée d'un centre et d'un ecartType
	 * @param centre
	 * @param ecartType
	 * @param nbPoints
	 */
	public static void writeNuage(Point3D centre, double ecartType, int nbPoints) {
		Random alea = new Random();
		System.out.println("");
		System.out.print("writeNuage ");
		double teta, phi, dist;
		double i, j , k;
		for (int n=0; n<nbPoints; n++){
			teta = alea.nextDouble()*2*Math.PI;
			phi = alea.nextDouble()*Math.PI;
			dist = Math.abs( alea.nextGaussian() * ecartType );
			i = centre.getX() + dist * Math.sin(phi) * Math.cos(teta);
			j = centre.getY() + dist * Math.sin(phi) * Math.sin(teta);
			k = centre.getZ() + dist * Math.cos(phi);
			buffer += i + " " + j + " " + k + "\n";
		}
		GeometricFile3DGenerator.nbPoints += nbPoints;
		System.out.println("");
		System.out.println("nuage de centre ("+centre.getX()+", "+centre.getY()+", "+centre.getZ()+") d'ecart-type "+ecartType);
	}

	/**
	 * Ecris un nuage de points spherique à partir de la donnée d'un centre et d'un rayon significatif et d'un ecart-type
	 * @param centre
	 * @param rayon
	 * @param ecartType
	 * @param nbPoints
	 */
	public static void writeNuageSpherique(Point3D centre, double ecartType, double rayon , int nbPoints) {
		Random alea = new Random();
		System.out.println("");
		System.out.print("writeNuageSpherique ");
		double teta, phi, dist;
		double i, j , k;
		for (int n=0; n<nbPoints; n++){
			teta = alea.nextDouble()*2*Math.PI;
			phi = alea.nextDouble()*Math.PI;
			dist = Math.abs( alea.nextGaussian() * ecartType + rayon );
			i = centre.getX() + dist * Math.sin(phi) * Math.cos(teta);
			j = centre.getY() + dist * Math.sin(phi) * Math.sin(teta);
			k = centre.getZ() + dist * Math.cos(phi);
			buffer += i + " " + j + " " + k + "\n";
		}
		GeometricFile3DGenerator.nbPoints += nbPoints;
		System.out.println("");
		System.out.println("nuage spherique de centre ("+centre.getX()+", "+centre.getY()+", "+centre.getZ()+") de rayon significatif "+rayon);
	}
	
	
	/**
	 * Ecris une surface spherique à partir de la donnée d'un centre et d'un rayon
	 * @param centre
	 * @param rayon
	 * @param pasAngulaire
	 */
	public static void writeSurfaceSpherique(Point3D centre, double rayon, double pasAngulaire) {
		System.out.println("");
		System.out.print("writeSurfaceSpherique ");
		double teta, phi;
		double i, j , k;
		for( teta=0; teta<2*Math.PI; teta+=pasAngulaire){
			for( phi=0; phi<Math.PI; phi+=pasAngulaire){
				i = centre.getX() + rayon * Math.sin(phi) * Math.cos(teta);
				j = centre.getY() + rayon * Math.sin(phi) * Math.sin(teta);
				k = centre.getZ() + rayon * Math.cos(phi);
				buffer += i + " " + j + " " + k + "\n";
				nbPoints++;
			}
		}
		System.out.println("");
		System.out.println("surface spherique de centre ("+centre.getX()+", "+centre.getY()+", "+centre.getZ()+") de rayon "+rayon);
	}
}