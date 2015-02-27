import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import jxl.Sheet;
import jxl.Workbook;


public class ReorganizeParameters {
	public static String tableurFileName = "DocumentationDoxygen.xls";

	public static void main(String[] args){
		try{
			String argument = args[0];
			File f=new File(argument);
			String contenu = AutoRefactor.loadFile(f);
			ArrayList<String> fichiers = new ArrayList<String>();
			fichiers.add(contenu);

			// 0) recuperer liste methodes
			File tabFile = new File(tableurFileName);
			ArrayList<String> listeMethodes = new ArrayList<String>();
			Workbook workbook = Workbook.getWorkbook(tabFile);
			Sheet sheet1 = workbook.getSheet(1);
			for (int i=2; i<64; i++) {
				listeMethodes.add(sheet1.getCell(3, i).getContents());
				//System.out.println(sheet1.getCell(1, i).getContents());
			}
			Sheet sheet2 = workbook.getSheet(2);
			for( String methode: listeMethodes){
				ArrayList<String> listeVariablesTriees = new ArrayList<String>();
				ArrayList<String> listeVariablesNonTriees = new ArrayList<String>();
				int indexMethode=-1;
				for( String fichier: fichiers){
					indexMethode = fichier.indexOf("SUBROUTINE "+methode);
					if( indexMethode == -1 ) continue;

					// 1) recuperer bon ordre des parametres
					int indexDebut = fichier.indexOf("!#### Parameters ####",indexMethode);
					int indexFin = fichier.indexOf("!#### Variables  ####",indexDebut);
					int i = fichier.indexOf("!###########################################",indexDebut);
					if( indexFin==-1 || i<indexFin){
						indexFin = i;
					}
					String declarationParametres = fichier.substring(indexDebut, indexFin);
					String[] declarationParametresbyLine = declarationParametres.split("\\n");
					for( String line :  declarationParametresbyLine){
						int indexDebutVariable = line.lastIndexOf(":");
						if( indexDebutVariable == -1) continue;
						indexDebutVariable+=2;
						//System.out.println(((line.substring(indexDebutVariable)).split(" |!"))[0]);
						listeVariablesTriees.add(((line.substring(indexDebutVariable)).split(" |!"))[0]);
					}

					// 2) recupere ordre actuel des parametres
					int indexParentheseOuvrante = fichier.indexOf("(",indexMethode);
					int indexParentheseFermante = fichier.indexOf(")",indexMethode);
					String parametres = fichier.substring(indexParentheseOuvrante+1,indexParentheseFermante);
					String[] parametresSplit= parametres.split("\\s|\\,");
					for( String parametre : parametresSplit){
						/*
						for( String s : listeVariablesTriees ){
							if(s.equals(parametre)){
								listeVariablesNonTriees.add(parametre);
							}
						}
						/*/
						if( parametre.length()>0 && ! parametre.contains("&")){
							listeVariablesNonTriees.add( parametre);
						}
						//*/
					}
					// 3) Changer la signature
					StringBuffer fichierBuffer = new StringBuffer(fichier);
					fichierBuffer.delete(indexParentheseOuvrante+1, indexParentheseFermante);
					String resultat = "";
					boolean first = true;
					for( String parametre: listeVariablesTriees){
						if( first ){
							resultat += parametre;
							first = false;
						}
						else {
							resultat += ","+parametre;
						}
					}
					fichierBuffer.insert(indexParentheseOuvrante+1, resultat);
					fichiers.set(fichiers.indexOf( fichier ),fichierBuffer.toString());
					break;
				}
				if( indexMethode==-1) continue;

				// 4) Changer les CALL
				for( String fichier: fichiers){
					StringBuffer fichierBuffer = new StringBuffer(fichier);
					int indexCall=0;
					while( (indexCall=fichier.indexOf("CALL "+methode,indexCall+1)) != -1 ){
						int indexParentheseOuvrante = fichier.indexOf("(",indexCall);
						int indexParentheseFermante = fichier.indexOf(")",indexCall);
						String parametres = fichier.substring(indexParentheseOuvrante+1,indexParentheseFermante);
						String[] parametresSplit= parametres.split("\\s|\\,");
						ArrayList<String> listeExpressions = new ArrayList<String>();
						for( String parametre : parametresSplit){
							if( parametre.length()>0 && ! parametre.contains("&")){
								listeExpressions.add( parametre);
							}
						}
						fichierBuffer.delete(indexParentheseOuvrante+1, indexParentheseFermante);
						ArrayList<String> listeExpressionsTriees = new ArrayList<String>();
						for( String parametre: listeVariablesTriees){
							for(String ss : listeVariablesNonTriees){
								System.out.println("1: "+ss);
							}
							for(String ss : listeExpressions){
								System.out.println("2: "+ss);
							}
							for(String ss : listeExpressionsTriees){
								System.out.println("3: "+ss);
							}
							System.out.println();
							listeExpressionsTriees.add(
									listeExpressions.get(
											listeVariablesNonTriees.indexOf(parametre)));
						}
						String resultat = "";
						boolean first = true;
						for( String parametre: listeExpressionsTriees){
							if( first ){
								resultat += parametre;
								first = false;
							}
							else {
								resultat += ","+parametre;
							}
						}
						fichierBuffer.insert(indexParentheseOuvrante+1, resultat);
					}
					fichier = fichierBuffer.toString();
				}
			}
			System.out.println(fichiers.get(0));
			//AutoRefactor.writeFile(f.getName(), "", contenu);
		}
		catch( Exception e){
			e.printStackTrace();
		}
	}

}
