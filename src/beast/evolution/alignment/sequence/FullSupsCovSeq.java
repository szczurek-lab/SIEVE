package beast.evolution.alignment.sequence;

import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FullSupsCovSeq extends ScsSequence {

    public FullSupsCovSeq() {
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
    }

    @Override
    protected void initializeCompatibleDataTypes() {
        compatibleDataTypes = new String[]{FullSupsCov.class.getName()};
    }

}
