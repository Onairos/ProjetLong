import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;


public class ReorganizeDeclarationsLines {
	
	public final static String annonceParametre = "!====";
	public final static String annonceVariable = "!#### Variables";
	
	public static void main(String[] args) {
		String argument = args[0];

		if( argument.equals("all")){
			refactor("clusters.f90");
			refactor("gmsh2cluster.f90");
			refactor("module_calcul.f90");
			refactor("module_decoupe.f90");
			refactor("module_embed.f90");
			refactor("module_entree.f90");
			refactor("module_MPI.f90");
			refactor("module_solve.f90");
			refactor("module_sortie.f90");
			refactor("module_sparse.f90");
			refactor("module_structure.f90");
			refactor("module_teste_clusters.f90");
			refactor("module_visuclusters.f90");
			refactor("module_visuclusters_gmsh.f90");
			refactor("module_visuclusters_paraview.f90");
			refactor("module_visuclusters_structure.f90");
			refactor("teste_clusters.f90");
			refactor("visuclusters.f90");
			refactor("visudecoup.f90");
		}
		else{
			refactor( argument);
		}
		
		System.out.println("###################################");
		System.out.println("fini ;)");
		System.out.println("###################################");
		
	}
	
	public static void refactor(String fileName) {
		//Récupération du fichier à réusiner
		File f = new File(fileName);
		String oldFileContent = loadFile( f );
		String newFileContent = "";
		
		//Récupération ligne par ligne
		String[] lines = oldFileContent.split("\\n");
		
		//Recherche des débuts et fins locaux où le tri est nécessaire
		ArrayList<Integer> indicesDebut = new ArrayList<>();
		ArrayList<Integer> indicesFin = new ArrayList<>();
		for (int i=0; i<lines.length ; i++) {
			if ( lines[i].contains("::") ) {
				//on a trouvé le départ d'un tri local, c'est la ligne suivante
				indicesDebut.add(i);
				//on prend le type associé à la variable
				String type = getNomType(lines[i]);
				//il faut ensuite trouver une ligne ne contenant pas ce type (ce sera la fin de ce tri local)
				for ( int j = i+1; j<lines.length ; j++) {
					if ( !lines[j].contains( type ) ) {
						//c'est la fin!!
						indicesFin.add(j-1);
						//on change i pour faire la recherche à partir de la ligne suivante
						i = j-1; //tout de suite après on a i++ donc i=j
						break;
					}
				}
			}
		}
		
		/*System.out.println( indicesDebut );
		System.out.println("-------------------------------------");
		System.out.println( indicesFin );*/
		
		//Effectuer les tris à bulle
		for (int i=0; i<indicesDebut.size(); i++) {
			doTriBulle( lines, indicesDebut.get(i), indicesFin.get(i) );
		}
		
		//Ecrire sur le fichier
		for (String s : lines) {
			newFileContent += s;
			newFileContent += "\n";
		}
		System.out.println(newFileContent);
		writeFile(fileName,"",newFileContent);
		//System.out.println("#######################");
		/*for (String s : lines) {
			System.out.println(s);
		}*/

	}
	
	
	public static void swapLines(String[] lines, int ind1, int ind2) {
		//DEBUG
		System.out.println("######## Following lines will be swapped :");
		System.out.println(ind1+1 + lines[ind1]);
		System.out.println(ind2+1 + lines[ind2]);
		System.out.println("########");
		
		String tmp = lines[ind1];
		lines[ind1] = lines[ind2];
		lines[ind2] = tmp;
	}
	
	/**
	 * A partir d'une ligne de type déclaration F90 <nomType> :: nomVar
	 * 		récupère le nom du type
	 * Précondition : ligne de type déclaration avec normes du PL, cad une seule variable par ligne avec " :: " comme séparateur, et pas d'espace après le nom de la variable
	 * @param line la ligne
	 */
	public static String getNomType(String line) {
		String[] tab;
		tab = line.split("(\\s+)" + "::" + "(\\s+)");
		return tab[0];
	}
	
	/**
	 * A partir d'une ligne de type déclaration F90 <nomType> :: nomVar
	 * 		récupère le nom de la variable
	 * Précondition : ligne de type déclaration avec normes du PL, cad une seule variable par ligne avec " :: " comme séparateur, et pas d'espace après le nom de la variable
	 * @param line la ligne
	 */
	public static String getNomVariable(String line) {
		String[] tab;
		tab = line.split("(\\s+)" + "::" + "(\\s+)");
		//Enlever un commentaire éventuel
		String retour = "";
		try {
			retour = tab[1];
			String[] coupe = retour.split("(\\s*)" + "!");
			retour = coupe[0];
		} catch (ArrayIndexOutOfBoundsException e) {
			//System.out.println("## ERROR sur la ligne : ##");
			//System.out.println(line);
			e.printStackTrace();
			//System.out.println("## FIN ERROR##");
			//System.out.println();
		}
		return retour;
	}
	
	
	public static void doTriBulle(String[] lines, int debut, int fin) {
		if (debut!=fin) {
			//Voir Wikipedia pour l'implantation
			for (int i=fin-1; i>=debut; i--) {
				for (int j=debut; j<=i; j++) {
					//System.out.println("Boucle doTriBulle : i=" + i + " et j=" + j);
					if ( mustBeSwapped(lines, j)  ) {
						swapLines(lines, j, j+1);
					}
				}
			}
		}
	}
	
	public static boolean mustBeSwapped(String[] lines, int j) {
		boolean retour = false;
		String var1 = getNomVariable( lines[j] );
		String var2 = getNomVariable( lines[j+1] );
		if ( var1.compareToIgnoreCase(var2) > 0 ) { //sir var 1 est situé après var2 dans le dictionnaire
			retour = true;
		}
		return retour;
	}
	
	
	 /**
	    * Loads the specified file into a String representation
	    * @author Stephane Nicoll - Infonet FUNDP
	    * @version 0.1
	    */
		
	    public static String loadFile(File f) {
	        try {
	           BufferedInputStream in = new BufferedInputStream(new FileInputStream(f));
	           StringWriter out = new StringWriter();
	           int b;
	           while ((b=in.read()) != -1)
	               out.write(b);
	           out.flush();
	           out.close();
	           in.close();
	           return out.toString();
	        }
	        catch (IOException ie)
	        {
	             ie.printStackTrace(); 
	        }
			return null;
	    }
	    
	    public static void writeFile(String nomFichier, String extension, String s){
	    	FileWriter fw = null;
			try {
				fw = new FileWriter(nomFichier+extension, false);
				fw.write(s, 0, s.length());
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

}
