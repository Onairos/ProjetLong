import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;


public class AutoRefactor {

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
	}
		
	public static void refactor(String arg){
		File f = new File(arg);
		String s = loadFile( f );
		writeFile(arg,".old",s);
		
		// mots-clefs generaux
		s = s.replaceAll("program", "PROGRAM");
		s = s.replaceAll("module ", "MODULE ");
		s = s.replaceAll("contains", "CONTAINS");
		s = s.replaceAll("subroutine", "SUBROUTINE");
		s = s.replaceAll("function", "FUNCTION");
		s = s.replaceAll("continue", "CONTINUE");
		s = s.replaceAll("return", "RETURN");
		s = s.replaceAll("call", "CALL");
		s = s.replaceAll("use ", "USE ");
		s = s.replaceAll("include", "INCLUDE");
		s = s.replaceAll("deallocate", "DEALLOCATE");
		s = s.replaceAll("allocate", "ALLOCATE");
		s = s.replaceAll("allocatable", "ALLOCATABLE");
		s = s.replaceAll("open", "OPEN");
		s = s.replaceAll("close", "CLOSE");
		s = s.replaceAll("file=", "FILE=");
		s = s.replaceAll("unit=", "UNIT=");
		s = s.replaceAll("inquire", "INQUIRE");
		s = s.replaceAll("exist ", "EXIST ");
		s = s.replaceAll("exist=", "EXIST=");
		s = s.replaceAll("exist\\(", "EXIST(");
		s = s.replaceAll("read", "READ");
		s = s.replaceAll("write", "WRITE");
		s = s.replaceAll("print", "PRINT");
		s = s.replaceAll("implicit none", "IMPLICIT NONE");
		s = s.replaceAll("implicit", "IMPLICIT");
		s = s.replaceAll("intent\\(in\\)", "INTENT(IN)");
		s = s.replaceAll("intent\\(out\\)", "INTENT(OUT)");
		s = s.replaceAll("external", "EXTERNAL");
		s = s.replaceAll("intrinsic", "INTRINSIC");
		
		// types
		s = s.replaceAll("type ", "TYPE ");
		s = s.replaceAll("kind", "KIND");
		s = s.replaceAll("!TYPE", "!type");
		s = s.replaceAll("type\\(", "TYPE(");
		s = s.replaceAll("integer", "INTEGER");
		s = s.replaceAll("real", "REAL");
		s = s.replaceAll("double precision", "DOUBLE PRECISION");
		s = s.replaceAll("Double precision", "DOUBLE PRECISION");
		s = s.replaceAll("complex", "COMPLEX");
		s = s.replaceAll("character", "CHARACTER");
		s = s.replaceAll("CHARACTERes", "caracteres");
		s = s.replaceAll("logical", "LOGICAL");
		s = s.replaceAll("parameter", "PARAMETER");
		s = s.replaceAll("common", "COMMON");
		s = s.replaceAll("pointer", "POINTER");
		s = s.replaceAll("dimension", "DIMENSION");
		s = s.replaceAll("reDIMENSION", "redimension");
		s = s.replaceAll("\\> DIMENSION", "> dimension");
		s = s.replaceAll("\\>DIMENSION", ">dimension");
		s = s.replaceAll("!DIMENSION", "!dimension");
		s = s.replaceAll("DIMENSIONs", "dimensions");
		
		// structures de controle
		s = s.replaceAll("goto", "GOTO");
		s = s.replaceAll("go to", "GOTO");
		s = s.replaceAll("end if", "ENDIF");
		s = s.replaceAll("endif", "ENDIF");
		s = s.replaceAll("#ENDIF", "#endif");
		s = s.replaceAll("end do", "ENDDO");
		s = s.replaceAll("enddo", "ENDDO");
		s = s.replaceAll("then", "THEN");
		s = s.replaceAll("else if", "ELSEIF");
		s = s.replaceAll("elseif", "ELSEIF");
		s = s.replaceAll("else", "ELSE");
		s = s.replaceAll("select", "SELECT");
		s = s.replaceAll("case default", "CASE DEFAULT");
		s = s.replaceAll("case", "CASE");
		s = s.replaceAll(" do", " DO");
		s = s.replaceAll("\\tdo", "	DO");
		s = s.replaceAll("DOm", "dom");
		s = s.replaceAll("DOn", "don");
		s = s.replaceAll("DOu", "dou");
		s = s.replaceAll(" if", " IF");
		s = s.replaceAll("\\tif", "	IF");
		s = s.replaceAll("while", "WHILE");
		
		//booleens
		s = s.replaceAll("\\.true\\.", ".TRUE.");
		s = s.replaceAll("\\.false\\.", ".FALSE.");
		s = s.replaceAll("\\.eq\\.", ".EQ.");
		s = s.replaceAll("\\.ne\\.", ".NE.");
		s = s.replaceAll("\\.gt\\.", ".GT.");
		s = s.replaceAll("\\.lt\\.", ".LT.");
		s = s.replaceAll("\\.ge\\.", ".GE.");
		s = s.replaceAll("\\.le\\.", ".LE.");
		s = s.replaceAll("\\.or\\.", ".OR.");
		s = s.replaceAll("\\.and\\.", ".AND.");
		s = s.replaceAll("\\.xor\\.", ".XOR.");
		s = s.replaceAll("\\.not\\.", ".NOT.");
		
		// fin
		s = s.replaceAll("stop", "STOP");
		s = s.replaceAll("end ", "END ");
		s = s.replaceAll("end=", "END=");
		s = s.replaceAll("#END ", "#end ");
		
		writeFile(arg,"",s);
		System.out.println(s);
		System.out.println("fini");
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
