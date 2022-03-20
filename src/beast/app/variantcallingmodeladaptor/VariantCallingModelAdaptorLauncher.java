package beast.app.variantcallingmodeladaptor;

import beast.app.beastapp.BeastLauncher;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

public class VariantCallingModelAdaptorLauncher extends BeastLauncher {

    public static void main(String[] args) throws NoSuchMethodException, SecurityException, ClassNotFoundException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, IOException {
        if (javaVersionCheck("VariantCaller")) {
            String classpath = getPath(false, null);
            run(classpath, "beast.app.variantcallingmodeladaptor.VariantCallingModelAdaptor", args);
        }
    } // main

}
