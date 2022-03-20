package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.alignment.sequence.ReadCountsDetails;
import com.google.common.primitives.Ints;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;

@Description("Coverage and support data type")
public class CovSup extends ReadCounts.Base {

    /**
     * constructor for testing purpose
     */
    public CovSup() {
    }

    /**
     * constructor for testing purpose
     *
     * @param states
     */
    public CovSup(int states) {
        stateCountInput.setValue(states, this);
        initAndValidate();
    }

    /**
     * convert read counts from string to integer and properly store them
     *
     * @param data    data in a string format
     * @param details detailed data for later reference
     * @return integer read counts
     */
    @Override
    public List<int[]> stringToInteger(
            String data,
            @NotNull List<ReadCountsDetails> details
    ) throws IllegalArgumentException {
        List<int[]> sequence = new ArrayList<>();

        // remove spaces
        data = data.replaceAll("\\s+", "");

        String[] lociReadCounts = data.split(";");
        for (String str_1 : lociReadCounts) {
            String[] locusReadCountsStr = str_1.split(",");

            // initialize and sanity check stateCount
            if (stateCount == 0 && locusReadCountsStr.length == 2)
                stateCount = 2;

            // make sure the stateCount is the same everywhere in the data
            if (stateCount != locusReadCountsStr.length)
                throw new IllegalArgumentException("The number of states is not unanimous across the data. Expect " +
                        "to be " + stateCount + ", " + locusReadCountsStr.length + " found.");

            // convert string to integer
            List<Integer> locusReadCountsInt = new ArrayList<>();
            for (String str_2 : locusReadCountsStr) {
                int integer;

                // mark non-integer data type as -1
                try {
                    integer = Integer.parseInt(str_2.trim());
                } catch (NumberFormatException e) {
                    integer = -1;
                }

                if (integer < 0)
                    throw new IllegalArgumentException("Read counts should be non-negative integers. Negative " +
                            "integer or non-integer data type found in " + str_2);
                else
                    locusReadCountsInt.add(integer);
            }

            details.add(
                    new ReadCountsDetails(
                            new char[]{'N'},
                            new int[]{locusReadCountsInt.get(1)},
                            locusReadCountsInt.get(0) - locusReadCountsInt.get(1),
                            locusReadCountsInt.get(0)
                    )
            );

            sequence.add(Ints.toArray(locusReadCountsInt));
        }

        return sequence;
    } // stringToInteger

    @Override
    public int getSequencingCoverage(int[] counts) {
        return counts[0];
    }

    @Override
    public String getTypeDescription() {
        return "coverage-support";
    } // getTypeDescription

} // class CovSup
