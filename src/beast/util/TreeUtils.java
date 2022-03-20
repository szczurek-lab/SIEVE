package beast.util;

import beast.evolution.tree.Node;

import java.util.Set;

public class TreeUtils {

    public static void processMetaData(Node node) {
        for (Node child : node.getChildren()) {
            processMetaData(child);
        }

        Set<String> metaDataNames = node.getMetaDataNames();
        if (metaDataNames != null && !metaDataNames.isEmpty()) {
            StringBuilder metadata = new StringBuilder();

            for (String name : metaDataNames) {
                Object value = node.getMetaData(name);
                metadata.append(name).append("=");

                if (value instanceof Object[]) {
                    Object[] values = (Object[]) value;

                    metadata.append("{");
                    for (int i = 0; i < values.length; i++) {
                        metadata.append(values[i].toString());
                        if (i < values.length - 1) {
                            metadata.append(",");
                        }
                    }
                    metadata.append("}");

                } else {
                    metadata.append(value.toString());
                }
                metadata.append(",");
            }

            metadata = new StringBuilder(metadata.substring(0, metadata.length() - 1));
            node.metaDataString = metadata.toString();
        }
    } // processMetaData

}
