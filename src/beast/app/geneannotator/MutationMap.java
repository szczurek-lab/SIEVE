package beast.app.geneannotator;

import beast.app.utils.ChromosomeLabel;

import java.util.Arrays;
import java.util.Comparator;

public class MutationMap {

    private final static String INFO_FORMAT_ERROR = "The information used to initialize should be of length 4.";
    private final static String POSITION_FORMAT_ERROR = "The second and third components should be numbers.";

    private final ChromosomeLabel chromosome;
    private final long startPos;
    private final long endPos;
    private final String geneName;

    public MutationMap(String[] mutationInfo) {
        if (mutationInfo.length != 4) {
            throw new RuntimeException(INFO_FORMAT_ERROR + Arrays.toString(mutationInfo));
        }

        this.chromosome = new ChromosomeLabel(mutationInfo[0]);

        try {
            this.startPos = Long.parseLong(mutationInfo[1]);
            this.endPos = Long.parseLong(mutationInfo[2]);
        } catch (NumberFormatException e) {
            throw new RuntimeException(POSITION_FORMAT_ERROR + Arrays.toString(mutationInfo));
        }

        this.geneName = mutationInfo[3];
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MutationMap that = (MutationMap) o;

        if (startPos != that.startPos) return false;
        if (endPos != that.endPos) return false;
        return this.chromosome.equals(that.chromosome);
    } // equals

    @Override
    public int hashCode() {
        int result = chromosome.hashCode();
        result = 31 * result + (int) (startPos ^ (startPos >>> 32));
        result = 31 * result + (int) (endPos ^ (endPos >>> 32));
        return result;
    } // hashCode

    @Override
    public String toString() {
        return this.chromosome.getFullLabel() + "," +
                this.startPos + "," +
                this.endPos + "," +
                this.geneName;
    } // toString

    public ChromosomeLabel getChromosome() {
        return this.chromosome;
    }

    public long getStartPos() {
        return this.startPos;
    }

    public long getEndPos() {
        return this.endPos;
    }

    public String getGeneName() {
        return this.geneName;
    }

    static class ListComparator implements Comparator<MutationMap> {

        @Override
        public int compare(MutationMap o1, MutationMap o2) throws RuntimeException {

            // chromosome label
            final int chrResult = o1.chromosome.compareTo(o2.chromosome);
            if (chrResult != 0) {
                return chrResult;
            }

            // start position
            if (o1.startPos > o2.startPos) {
                return 1;
            } else if (o1.startPos < o2.startPos) {
                return -1;
            }

            // end position
            if (o1.endPos > o2.endPos) {
                return 1;
            } else if (o1.endPos < o2.endPos) {
                return -1;
            }

            // gene name
            final int gn = o1.geneName.toUpperCase().compareTo(o2.geneName.toUpperCase());
            if (gn > 0) {
                return 1;
            } else if (gn < 0) {
                return -1;
            }

            return 0;
        }

    }

}
