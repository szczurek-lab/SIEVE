package beast.evolution.alignment.sequence;

import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.ReadCounts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CovSupSeq extends ScsSequence {

    public CovSupSeq() {
    }

    /**
     * constructor for testing purpose
     *
     * @param taxon
     * @param sequence
     */
    public CovSupSeq(String taxon, String sequence) {
        taxonInput.setValue(taxon, this);
        dataInput.setValue(sequence, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.readCountsDetails = new ArrayList<>();
    }

    @Override
    public List<int[]> getReadCounts(ReadCounts dataType) {
        if (!isCompatibleDataType(dataType))
            throw new IllegalArgumentException("Error! Incompatible datatype. Only allowed: " +
                    Arrays.toString(compatibleDataTypes));

        return dataType.stringToInteger(
                dataInput.get().trim().replaceAll("\\s+", ""),
                this.readCountsDetails
        );
    } // getReadCounts

    @Override
    protected void initializeCompatibleDataTypes() {
        compatibleDataTypes = new String[]{CovSup.class.getName()};
    }

}
