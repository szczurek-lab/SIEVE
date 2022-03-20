package beast.evolution.rawreadcountsmodel.seqcovmodel;

import beast.evolution.alignment.ScsAlignment;

import java.io.PrintStream;
import java.util.List;

public class ReducedSeqCovModel extends SeqCovModelInterface.Base {


    //***********************************************
    //*                  Variables                  *
    //***********************************************



    //**********************************************
    //*             Overridden methods             *
    //**********************************************

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (allelicSeqCovInput.get() != null ||
                allelicSeqCovRawVarInput.get() != null ||
                allelicSeqCovArrayInput.get() != null ||
                allelicSeqCovRawVarArrayInput.get() != null)
            throw new IllegalArgumentException("The reduced sequencing coverage model is used. Attributes 'allelicSeqCov', " +
                    "'allelicSeqCovRawVar', 'allelicSeqCovArr', or 'allelicSeqCovRawVarArr' should not be specified.");

        if (initSharedAllelicSeqCovAndRawVarInput.get() != null)
            throw new IllegalArgumentException("\"initSharedAllelicSeqCovAndRawVar\" is only valid for ExploredSharedAllelicSeqCovModel.");
    } // initAndValidate

    /**
     * Initialize parameters with passed arguments from ScsTreeLikelihood.
     * Must be called by every inheriting class at the beginning.
     *
     * @param alignment         apparently
     * @param matrixNum         number of matrices
     * @param modeledAlleles    alleles the model of evolution models
     * @param useLogPartials apparently
     */
    @Override
    public void deeplyInitialize(
            final ScsAlignment alignment,
            final int matrixNum,
            final int[] modeledAlleles,
            final boolean useLogPartials
    ) {
        nrOfPatterns = alignment.getPatternCount();
        nrOfTaxa = alignment.getTaxonCount();
        nrOfMatrices = 1;
        this.useLogPartials = useLogPartials;

        needToUpdate = false;
    } // deeplyInitialize

    /**
     * Duplicate some important variables entirely or within ranged patterns.
     * Used when there are more than one thread.
     * Must be called by every inheriting class at the beginning.
     *
     * @param target     target instance storing variables
     * @param targetData data corresponding to the target raw read counts model
     */
    @Override
    public void duplicate(
            SeqCovModelInterface.Base target,
            final ScsAlignment targetData
    ) {
        target.nrOfMatrices = 1;
        target.nrOfPatterns = targetData.getPatternCount();
        target.nrOfTaxa = this.nrOfTaxa;
        target.useLogPartials = this.useLogPartials;

        needToUpdate = false;
    } // duplicate

    /**
     * For the purpose of debugging.
     * Only the most recent accepted allelic sequencing coverage and raw variance should be stored.
     */
    @Override
    public void storeAllelicInfo() {
        throw new IllegalArgumentException("Unsupported function.");
    } // storeAllelicInfo

    /**
     * Compute sequencing coverage likelihood under negative-binomial distribution.
     *
     * @param seqCov       sequencing coverage for current taxa under current pattern
     * @param matrixIndex  index to matrix
     * @param patternIndex index to pattern
     * @param taxonIndex   index to taxa
     */
    public void computeSeqCovLikelihood(
            final int seqCov,
            final int matrixIndex,
            final int patternIndex,
            final int taxonIndex
    ) {
        // Nothing needs to be computed.
    } // computeSeqCovLikelihood

    /**
     * Compute allelic sequencing coverage and its raw variance.
     * Should be called at the post processing stage of computing tree likelihood.
     *
     * @param changedPatterns matrices and patterns which should be updated: (matrixIndex, patternIndex)
     */
    @Override
    public void updateAllelicSeqCovAndRawVar(List<List<Integer>> changedPatterns) {
        throw new IllegalArgumentException("Unsupported function.");
    } // updateAllelicSeqCovAndRawVar

    @Override
    public double getSeqCovLikelihood(
            final int taxonIndex,
            final int index
    ) {
        return useLogPartials ? 0 : 1;
    } // getSeqCovLikelihood

    @Override
    public void setSeqCovLikelihoodsNodeForUpdate(int nodeIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    } // setSeqCovLikelihoodsNodeForUpdate

    @Override
    public void storeLikelihoods(int taxonIndex) {
        // Nothing needs to be done.
    } // storeLikelihoods

    @Override
    public void expandVariables(int taxonIndex) {
        // Nothing needs to be done.
    } // expandVariables

    @Override
    public void addMLNrOfExistingAllelesPatterns(int index, int[] values) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public void addMLNrOfSequencedAllelesWeights(int index, int value) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public List<int[]> getMLNrOfSequencedAllelesPatterns(int index) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public int getMLNrOfSequencedAllelesWeights(int index1, int index2) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public void resetMLNrOfSequencedAllelesPatternsAndWeights(int index) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public void logCovar(int index, PrintStream out) {
        throw new IllegalArgumentException("Unsupported function.");
    } // logCovar


    //*********************************************
    //*                  Methods                  *
    //*********************************************

}
