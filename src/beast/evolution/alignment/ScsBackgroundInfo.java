package beast.evolution.alignment;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;
import com.google.common.primitives.Longs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Description("Background information of single sequence read counts data in an alignment")
public class ScsBackgroundInfo extends BEASTObject {

    private enum BackgroundInfoType {
        Universal,
        CovSup,
        FullSupsCov
    }

    final private static String COV_SUP_FORMAT = "^variant$";
    final private static String FULL_SUPS_COV_FORMAT = "^variant[1-3]$";

    final private static Pattern COV_SUP_PATTERN = Pattern.compile(COV_SUP_FORMAT);
    final private static Pattern FULL_SUPS_COV_PATTERN = Pattern.compile(FULL_SUPS_COV_FORMAT);

    /**
     * default name of background information
     */
    final private static String DEFAULT_NAME = "coverage";

    /**
     * background names to choose from
     */
    final private static List<String> availableNames = Stream.of("variant1", "variant2", "variant3", "variant", "coverage", "normal").collect(Collectors.toList());

    final public Input<String> backgroundNameInput = new Input<>(
            "backgroundName",
            "name of background information, one of " + availableNames,
            DEFAULT_NAME,
            availableNames.toArray(new String[0])
    );

    final public Input<String> dataInput = new Input<>("value", "the details of background information " +
            "of current data set, in a default form of 'reads,occurrences;' (note that semicolon is used to " +
            "separate different reads).", Input.Validate.REQUIRED);

    private BackgroundInfoType backgroundInfoType = null;

    /**
     * name of background information
     */
    private String backgroundName = null;

    /**
     * store background information
     */
    private List<long[]> backgroundList = null;

    /**
     * the number of background information points
     */
    private long backgroundInfoPoints = 0;

    /**
     * constructor for testing purpose
     */
    public ScsBackgroundInfo() {
    }

    /**
     * constructor for testing purpose
     *
     * @param infoName
     * @param info
     */
    public ScsBackgroundInfo(String infoName, String info) {
        backgroundNameInput.setValue(infoName, this);
        dataInput.setValue(info, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        if (backgroundNameInput.get() == null)
            throw new IllegalArgumentException("Error! Background name missing.");

        this.backgroundName = backgroundNameInput.get();
        this.backgroundInfoPoints = 0;

        initBackgroundInfoType();

        initBackgroundInfo();
    } // initAndValidate

    private void initBackgroundInfoType() {
        if (FULL_SUPS_COV_PATTERN.matcher(this.backgroundName).matches())
            this.backgroundInfoType = BackgroundInfoType.FullSupsCov;
        else if (COV_SUP_PATTERN.matcher(this.backgroundName).matches())
            this.backgroundInfoType = BackgroundInfoType.CovSup;
        else
            this.backgroundInfoType = BackgroundInfoType.Universal;
    } // initBackgroundInfoType

    /**
     * save the background information into an ArrayList consisting of ImmutableList
     */
    private void initBackgroundInfo() throws IllegalArgumentException {
        // remove spaces
        String data = dataInput.get().trim().replaceAll("\\s+", "");

        this.backgroundList = new ArrayList<>();

        // separate different counts
        String[] ctsInfo = data.split(";");

        List<Long> readCountsOccurrences = new ArrayList<>();
        for (String str : ctsInfo) {
            readCountsOccurrences.clear();

            String[] indiInfo = str.split(",");

            // convert string to integer
            for (String str_2 : indiInfo) {
                long integer;

                // mark non-integer data type as -1
                try {
                    integer = Long.parseLong(str_2.trim());
                } catch (NumberFormatException e) {
                    integer = -1;
                }

                if (integer < 0)
                    throw new IllegalArgumentException("Read counts and corresponding occurrences should both be " +
                            "non-negative integers. Negative integer or non-integer data type found in " + str_2 +
                            ". (" + this.getClass().getName() + ")");
                else
                    readCountsOccurrences.add(integer);

            }

            // size of readCountsOccurrences should be 2
            if (readCountsOccurrences.size() != 2)
                throw new IllegalArgumentException("Background information should be in a form of " +
                        "'the number of read count + the corresponding unique occurrences', nothing more or less. " +
                        readCountsOccurrences.toString() + "found. (" + this.getClass().getName() + ")");

            this.backgroundInfoPoints += readCountsOccurrences.get(1);

            this.backgroundList.add(Longs.toArray(readCountsOccurrences));
        }

        // sort and check if current counts has appeared before
        backgroundList.sort(new CombiComparator());
    } // initBackgroundInfo

    public boolean sanityCheckDataType(ReadCounts dataType) {
        return (this.backgroundInfoType == BackgroundInfoType.Universal) ||
                (dataType instanceof CovSup && this.backgroundInfoType == BackgroundInfoType.CovSup) ||
                (dataType instanceof FullSupsCov && this.backgroundInfoType == BackgroundInfoType.FullSupsCov);
    } // sanityCheckDataType

    @Override
    public String toString() {
        if (this.backgroundName != null) {
            return dataInput.get();
        } else {
            return null;
        }
    } // toString

    public static List<String> getOrderedNamesCovSup() {
        return Stream.of("variant", "normal", "coverage").collect(Collectors.toList());
    } // getOrderedNamesCovSup

    public static List<String> getOrderedNamesFullSupsCov() {
        return Stream.of("variant1", "variant2", "variant3", "normal", "coverage").collect(Collectors.toList());
    } // getOrderedNamesFullSupsCov


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    /**
     * getter method for outside calls
     *
     * @return backgroundName
     */
    public String getBackgroundName() {
        return this.backgroundName;
    }

    /**
     * getter method for outside calls
     *
     * @return backgroundList
     */
    public List<long[]> getBackgroundList() {
        return this.backgroundList;
    }

    public long getBackgroundInfoPoints() {
        return this.backgroundInfoPoints;
    }

    static class CombiComparator implements Comparator<long[]> {
        @Override
        public int compare(long[] v1, long[] v2) {
            if (v1.length != v2.length)
                throw new IllegalArgumentException("Error! Arrays should have the same length.");

            for (int i = 0; i < v1.length; i++) {
                if (v1[i] > v2[i])
                    return 1;
                else if (v1[i] < v2[i])
                    return -1;
            }

            throw new IllegalArgumentException("Error! No duplicates allowed: " + Arrays.toString(v1) + "; " +
                    Arrays.toString(v2));
        }
    }

} // class ScsBackgroundInfo
