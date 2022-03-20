package beast.evolution.alignment.sequence;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.ReadCounts;

import java.util.List;

@Description("Single sequence read counts data in an alignment")
public abstract class ScsSequence extends BEASTObject {

    final public Input<String> taxonInput = new Input<>("taxon", "name of this cell of loci information",
            Input.Validate.REQUIRED);
    final public Input<String> dataInput = new Input<>("value", "the read counts of single-cell " +
            "sequencing data, in a default form of 'c,ma;c,ma;...', where 'c' stands for the coverage, and 'ma' " +
            "represents the support read counts of the alternative nucleotide with the greatest number of read " +
            "counts (note that the comma is used to separate read counts at one locus, and the semicolon is used to " +
            "separate different loci)", Input.Validate.REQUIRED);

    protected List<ReadCountsDetails> readCountsDetails;

    protected static String[] compatibleDataTypes;

    @Override
    public void initAndValidate() {
        initializeCompatibleDataTypes();
    } // initAndValidate

    abstract public List<int[]> getReadCounts(ReadCounts dataType);

    final public String getTaxon() {
        return taxonInput.get();
    } // getTaxon

    abstract protected void initializeCompatibleDataTypes();

    protected boolean isCompatibleDataType(ReadCounts dataType) {
        final String dataTypeDesc = dataType.getClass().getName();

        for (String i : compatibleDataTypes) {
            if (i.equals(dataTypeDesc))
                return true;
        }

        return false;
    } // isCompatibleDataType

    public int getCoverage(int index) {
        return this.readCountsDetails.get(index).getCoverage();
    } // getCoverage

    public char[] getAltNucs(int index) {
        return this.readCountsDetails.get(index).getAltNucs();
    } // getAltNucs

    public int getRefReads(int index) {
        return this.readCountsDetails.get(index).getRefReads();
    } // getRefReads

    public int getAltReads(int locusIndex, int altNucIndex) {
        return this.readCountsDetails.get(locusIndex).getAltReads(altNucIndex);
    } // getAltReads

} // class ScsSequence
