package beast.util;

import java.io.File;

public class FileNameProcessor {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public final static String delimiter = String.valueOf(File.separatorChar);


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static String getBaseName(final String fullName) {
        if (fullName != null && !fullName.equals("")) {
            String[] tmp = fullName.split(delimiter);
            return tmp[tmp.length - 1].split("\\.")[0];
        } else {
            return "";
        }
    } // getBaseName

    public static String getSuffix(final String fullName) {
        if (fullName != null && !fullName.equals("")) {
            String[] tmp1 = fullName.split(delimiter);
            String[] tmp2 = tmp1[tmp1.length - 1].split("\\.");
            if (tmp2 != null && tmp2.length > 1) {
                return "." + tmp2[tmp2.length - 1];
            } else {
                return "";
            }
        } else {
            return "";
        }
    } // getSuffix

    public static String getRootPath(final String fullName) {
        if (fullName != null && !fullName.equals("")) {
            String[] tmp1 = fullName.split(delimiter);
            String[] tmp2 = new String[tmp1.length - 1];

            if (tmp2.length == 0) {
                return "";
            } else if (tmp2.length > 0) {
                System.arraycopy(tmp1, 0, tmp2, 0, tmp2.length);
                return String.join(delimiter, tmp2) + "/";
            }
        }

        return "";
    } // getRootPath

}
