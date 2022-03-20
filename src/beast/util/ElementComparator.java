package beast.util;

import java.util.Arrays;
import java.util.List;

public class ElementComparator {

    /**
     * check whether an array is in a list of arrays
     * the element type of array should be primitive
     * otherwise the result is incorrect
     *
     * @param group  a list of arrays to be compared
     * @param target an array to be found
     * @param <E>    primitive data type
     * @return target is in group or not
     */
    public static <E> boolean compListArrays(final List<E[]> group, final E[] target) {
        if (group == null || group.size() == 0) {
            return false;
        }

        for (E[] e : group) {
            if (Arrays.equals(e, target)) {
                return true;
            }
        }

        return false;
    } // compListArrays

    /**
     * check whether an integer array is in a list of integer arrays
     *
     * @param group  a list of arrays to be compared
     * @param target an array to be found
     * @return target is in group or not
     */
    public static boolean compIntArrays(final List<int[]> group, final int[] target) {
        if (group == null || group.size() == 0)
            return false;

        for (int[] e : group)
            if (Arrays.equals(e, target))
                return true;

        return false;
    } // compIntArrays

}
