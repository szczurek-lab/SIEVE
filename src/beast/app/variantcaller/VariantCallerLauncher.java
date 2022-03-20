package beast.app.variantcaller;

import beast.app.beastapp.BeastLauncher;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

public class VariantCallerLauncher extends BeastLauncher {

    public static void main(String[] args) throws NoSuchMethodException, SecurityException, ClassNotFoundException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, IOException {
        if (javaVersionCheck("VariantCaller")) {
            String classpath = getPath(false, null);
            run(classpath, "beast.app.variantcaller.VariantCaller", args);
        }
    } // main

}
