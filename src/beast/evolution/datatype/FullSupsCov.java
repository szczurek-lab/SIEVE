package beast.evolution.datatype;

import beast.evolution.alignment.sequence.ReadCountsDetails;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FullSupsCov extends ReadCounts.Base {

    /**
     * convert read counts from string to integer and properly store them
     *
     * @param data    data in a string format
     * @param details detailed data for later reference
     * @return integer read counts
     */
    @Override
    public List<int[]> stringToInteger(String data, @NotNull List<ReadCountsDetails> details) {
        List<int[]> sequence = new ArrayList<>();

        // variant1, variant2, variant3, ref, cov
        int[] fullReads = new int[5];

        char[] nucs = new char[3]; // first 3 items
        int[] nucReads = new int[3]; // the next 3 items

        // remove spaces
        data = data.replaceAll("\\s+", "");

        String[] lociReadCounts = data.split(";");
        for (String str_1 : lociReadCounts) {
            String[] locusReadCountsStr = str_1.split(",");

            // initialize and sanity check stateCount
            if (stateCount == 0 && locusReadCountsStr.length == 7)
                stateCount = 4;

            // make sure the stateCount is the same everywhere in the data
            if (stateCount != locusReadCountsStr.length - 3)
                throw new IllegalArgumentException("The number of states is not unanimous across the data. Expect " +
                        "to be " + stateCount + ", " + (locusReadCountsStr.length - 3) + " found.");

            for (int i = 0; i < locusReadCountsStr.length; i++) {
                locusReadCountsStr[i] = locusReadCountsStr[i].trim();

                if (i < 3) {
                    if (locusReadCountsStr[i].length() != 1 || !Character.isAlphabetic(locusReadCountsStr[i].charAt(0)))
                        throw new IllegalArgumentException("Error! Only a character of alphabetic is expected, but " +
                                "observing " + locusReadCountsStr[i]);

                    nucs[i] = locusReadCountsStr[i].charAt(0);
                } else {
                    int integer;

                    // mark non-integer data type as -1
                    try {
                        integer = Integer.parseInt(locusReadCountsStr[i]);
                    } catch (NumberFormatException e) {
                        integer = -1;
                    }

                    if (integer < 0)
                        throw new IllegalArgumentException("Read counts should be non-negative integers. Negative " +
                                "integer or non-integer data type found in " + locusReadCountsStr[i]);

                    if (i < 6)
                        nucReads[i - 3] = integer;
                    else {
                        System.arraycopy(nucReads, 0, fullReads, 0, 3);
                        fullReads[3] = integer - Arrays.stream(nucReads).sum();
                        fullReads[4] = integer;

                        details.add(
                                new ReadCountsDetails(
                                        nucs,
                                        nucReads,
                                        fullReads[3],
                                        fullReads[4]
                                )
                        );
                    }
                }
            }

            sequence.add(Arrays.copyOf(fullReads, 5));
        }

        return sequence;
    }

    @Override
    public int getSequencingCoverage(int[] counts) {
        return counts[4];
    }

    /**
     * data type description
     *
     * @return data type name
     */
    @Override
    public String getTypeDescription() {
        return "full supports-coverage";
    }

}
