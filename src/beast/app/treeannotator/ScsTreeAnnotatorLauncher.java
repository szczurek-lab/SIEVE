package beast.app.treeannotator;

import beast.app.beastapp.BeastLauncher;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

/**
 * Loads beast.jar and launches ScsTreeAnnotator
 * <p>
 * This class should be compiled against 1.6 and packaged by itself. The
 * remained of BEAST can be compiled against Java 1.7 or higher
 **/
public class ScsTreeAnnotatorLauncher extends BeastLauncher {

    public static void main(String[] args) throws NoSuchMethodException, SecurityException, ClassNotFoundException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, IOException {
        if (javaVersionCheck("ScsTreeAnnotator")) {
            String classpath = getPath(false, null);
            run(classpath, "beast.app.treeannotator.ScsTreeAnnotator", args);
        }
    }

}
