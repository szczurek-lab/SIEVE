package beast.evolution.datatype;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.sequence.ReadCountsDetails;
import org.jetbrains.annotations.NotNull;

import java.util.List;

@Description(value = "An interface for read counts data")
public interface ReadCounts {

    /**
     * convert read counts from string to integer and properly store them
     *
     * @param data    data in a string format
     * @param details detailed data for later reference
     * @return integer read counts
     */
    List<int[]> stringToInteger(String data, @NotNull List<ReadCountsDetails> details);

    /**
     * data type description
     *
     * @return data type name
     */
    String getTypeDescription();

    int getSequencingCoverage(int[] counts);


    @Description(value = "Base implementation of ReadCounts datatype", isInheritable = false)
    abstract class Base extends BEASTObject implements ReadCounts {

        final public Input<Integer> stateCountInput = new Input<>("state", "total number of read " +
                "counts at a locus in a cell. Should be unanimous across the whole data.");

        public int stateCount = 0;

        @Override
        public void initAndValidate() {
            if (stateCountInput.get() != null && sanityCheck(stateCountInput.get()))
                stateCount = stateCountInput.get();
        } // initAndValidate

        /**
         * check state count's sanity.
         * only accept locusReadCountsInt.size() or stateCountInput.get() as argument.
         *
         * @param count number of reads data provided
         * @return boolean
         */
        public boolean sanityCheck(int count) {
            if (count != 4 && count != 2) {
                stateCount = 0;
                throw new IllegalArgumentException("Total number of states of read counts " + count +" found, " +
                        "either 2 and or 4 allowed.");
            }
            else {
                return true;
            }
        } // sanityCheck

        @Override
        public String toString() {
            return getTypeDescription();
        } // toString

    } // abstract class Base

} // interface ReadCounts
