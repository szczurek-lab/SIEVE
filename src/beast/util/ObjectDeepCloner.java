package beast.util;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class ObjectDeepCloner {

    private ObjectDeepCloner() {}

    public static Object deepCopy(final Object o) throws Exception {
        ObjectOutputStream out = null;
        ObjectInputStream in = null;

        try {
            ByteArrayOutputStream byteOut = new ByteArrayOutputStream();

            out = new ObjectOutputStream(byteOut);
            out.writeObject(o);
            out.flush();

            ByteArrayInputStream byteIn = new ByteArrayInputStream(byteOut.toByteArray());

            in = new ObjectInputStream(byteIn);
            return in.readObject();
        } finally {
            assert out != null;
            out.close();
            assert in != null;
            in.close();
        }
    }

}
