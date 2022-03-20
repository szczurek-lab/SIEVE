package beast.app.datacollector;

import beast.app.beastapp.BeastLauncher;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

public class DataCollectorLauncher extends BeastLauncher {

    public static void main(String[] args) throws NoSuchMethodException, SecurityException, ClassNotFoundException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, IOException {
        if (javaVersionCheck("DataCollector")) {
            String classpath = getPath(false, null);
            run(classpath, "beast.app.datacollector.DataCollector", args);
        }
    } // main

}
