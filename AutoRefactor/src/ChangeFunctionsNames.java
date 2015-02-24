import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;


public class ChangeFunctionsNames {
	
	public static int numSheet = 1;
	public static int colonneOld = 0;
	public static int colonneNew = 2;
	public static int colonneValid = 3;
	public static int indiceDebutLigne = 2;
	public static int indiceFinLigne = 63;
	public static String tableurFileName = "DocumentationDoxygen.xls";

	
	public static void main(String[] args) {
		String argument = args[0];
		
		//Recup du Excel
		File tabFile = new File(tableurFileName);
		ArrayList<String> oldNames = new ArrayList<String>();
		ArrayList<String> newNames = new ArrayList<String>();
		ArrayList<String> valids = new ArrayList<String>();
		
		try {
			Workbook workbook = Workbook.getWorkbook(tabFile);
			Sheet sheet = workbook.getSheet(numSheet);
			
			for (int i=indiceDebutLigne; i<indiceFinLigne; i++) {
				oldNames.add(sheet.getCell(colonneOld, i).getContents());
				newNames.add(sheet.getCell(colonneNew, i).getContents());
				valids.add(sheet.getCell(colonneValid, i).getContents());
			}
		} catch (BiffException | IOException e) {
			e.printStackTrace();
		}
		
		
		
		if( argument.equals("all")){
			refactor("clusters.f90", oldNames, newNames, valids);
			refactor("gmsh2cluster.f90", oldNames, newNames, valids);
			refactor("module_calcul.f90", oldNames, newNames, valids);
			refactor("module_decoupe.f90", oldNames, newNames, valids);
			refactor("module_embed.f90", oldNames, newNames, valids);
			refactor("module_entree.f90", oldNames, newNames, valids);
			refactor("module_MPI.f90", oldNames, newNames, valids);
			refactor("module_solve.f90", oldNames, newNames, valids);
			refactor("module_sortie.f90", oldNames, newNames, valids);
			refactor("module_sparse.f90", oldNames, newNames, valids);
			refactor("module_structure.f90", oldNames, newNames, valids);
			refactor("module_teste_clusters.f90", oldNames, newNames, valids);
			refactor("module_visuclusters.f90", oldNames, newNames, valids);
			refactor("module_visuclusters_gmsh.f90", oldNames, newNames, valids);
			refactor("module_visuclusters_paraview.f90", oldNames, newNames, valids);
			refactor("module_visuclusters_structure.f90", oldNames, newNames, valids);
			refactor("teste_clusters.f90", oldNames, newNames, valids);
			refactor("visuclusters.f90", oldNames, newNames, valids);
			refactor("visudecoup.f90", oldNames, newNames, valids);
		}
		else{
			refactor( argument, oldNames, newNames, valids);
		}
		
		System.out.println("###################################");
		System.out.println("fini ;)");
		System.out.println("###################################");
	}
		
	public static void refactor(String arg, ArrayList<String> oldNames, ArrayList<String> newNames, ArrayList<String> valids){
		File f = new File(arg);
		String content = loadFile( f );

			
			int i = 0; //c'est moche xD
			for (String s : oldNames) {
				if ("x".equals(valids.get(i))) {
					content = content.replaceAll(" " + s + "\\(", " " + newNames.get(i) + "(" );
					content = content.replaceAll(" " + s + "\\s", " " + newNames.get(i) + "\n" );
				}
				i++;
			}
						
		writeFile(arg,".new",content);
		//System.out.println(content);
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
