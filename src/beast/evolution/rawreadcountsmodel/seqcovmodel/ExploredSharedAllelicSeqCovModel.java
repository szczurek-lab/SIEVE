package beast.evolution.rawreadcountsmodel.seqcovmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.ScsAlignment;

import java.io.PrintStream;
import java.util.List;

public class ExploredSharedAllelicSeqCovModel extends SeqCovModelInterface.Base {


    //***********************************************
    //*                  Variables                  *
    //***********************************************


    //**********************************************
    //*             Overridden methods             *
    //**********************************************

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (inVariantCallingMode) {
            if (allelicSeqCovInput.get() != null || allelicSeqCovRawVarInput.get() != null ||
                    allelicSeqCovArrayInput.get() == null || allelicSeqCovRawVarArrayInput.get() == null)
                throw new IllegalArgumentException("Error! In variant calling mode, 'allelicSeqCovArray' and " +
                        "'allelicSeqCovRawVarArray' should be defined but not 'allelicSeqCov' and 'allelicSeqCovRawVar'.");

            double tmp1, tmp2;
            try {
                tmp1 = Double.parseDouble(allelicSeqCovArrayInput.get());
                tmp2 = Double.parseDouble(allelicSeqCovRawVarArrayInput.get());
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Error! Values of 'allelicSeqCovArray' and " +
                        "'allelicSeqCovRawVarArray' should be non-negative (float) numbers.");
            }

            if (tmp1 < 0 || tmp2 < 0)
                throw new IllegalArgumentException("Error! Values of 'allelicSeqCovArray' and " +
                        "'allelicSeqCovRawVarArray' should be non-negative (float) numbers.");

            allelicSeqCov = new RealParameter(allelicSeqCovArrayInput.get());
            allelicSeqCovRawVar = new RealParameter(allelicSeqCovRawVarArrayInput.get());
        } else {
            if (allelicSeqCovInput.get() == null || allelicSeqCovRawVarInput.get() == null ||
                    allelicSeqCovArrayInput.get() != null || allelicSeqCovRawVarArrayInput.get() != null)
                throw new IllegalArgumentException("Error! Not in variant calling mode, 'allelicSeqCov' and " +
                        "'allelicSeqCovRawVar' should be defined but not 'allelicSeqCovArray' and 'allelicSeqCovRawVarArray'.");

            if (allelicSeqCovInput.get().getValue() < 0.0)
                throw new IllegalArgumentException("'allelicSeqCov' should be defined as no smaller than 0.");
            else
                allelicSeqCov = allelicSeqCovInput.get();

            if (allelicSeqCovRawVarInput.get().getValue() < 0.0)
                throw new IllegalArgumentException("'allelicSeqCovRawVar' should be defined as no smaller than 0.");
            else
                allelicSeqCovRawVar = allelicSeqCovRawVarInput.get();
        }
    } // initAndValidate

    /**
     * Initialize parameters with passed arguments from ScsTreeLikelihood.
     * Must be called by every inheriting class at the beginning.
     *
     * @param alignment      apparently
     * @param matrixNum      number of matrices
     * @param modeledAlleles alleles the model of evolution models
     * @param useLogPartials apparently
     */
    @Override
    public void deeplyInitialize(
            final ScsAlignment alignment,
            final int matrixNum,
            final int[] modeledAlleles,
            final boolean useLogPartials
    ) {
        super.deeplyInitialize(
                alignment,
                matrixNum,
                modeledAlleles,
                useLogPartials
        );

        sizeFactors = new double[nrOfTaxa];
        processedSeqCov = new int[nrOfPatterns * nrOfTaxa];

        computeSizeFactors(alignment);
        processedSeqCov = null;

        needToUpdate = false;
        nrOfMatrices = 1;

        currentSeqCovLikelihoodsNodeIndex = new int[nrOfTaxa];
        storedSeqCovLikelihoodsNodeIndex = new int[nrOfTaxa];

        seqCovLikelihoods = new double[2][nrOfTaxa][nrOfPatterns * modeledAllelesSize];
    } // deeplyInitialize

    @Override
    public void duplicate(Base target, ScsAlignment targetData) {
        super.duplicate(target, targetData);

        // allocate memory and copy values
        target.sizeFactors = new double[nrOfTaxa];
        System.arraycopy(this.sizeFactors, 0, target.sizeFactors, 0, nrOfTaxa);

        target.seqCovLikelihoods = new double[2][nrOfTaxa][target.nrOfPatterns * modeledAllelesSize];

        target.currentSeqCovLikelihoodsNodeIndex = new int[nrOfTaxa];
        target.storedSeqCovLikelihoodsNodeIndex = new int[nrOfTaxa];
    } // duplicate

    /**
     * Set allelic sequencing coverage and raw variance with statistics from samples for a good starting point.
     * Only valid for ExploredSharedAllelicSeqCovModel.
     *
     * @param allelicSeqCov       apparently
     * @param allelicSeqCovRawVar apparently
     */
    @Override
    protected void setAllelicSeqCovAndRawVar(double allelicSeqCov, double allelicSeqCovRawVar) {
        if (initSharedAllelicSeqCovAndRawVarInput.get() != null && !initSharedAllelicSeqCovAndRawVarInput.get()) return;

        if (this.allelicSeqCov.isEstimated()) {
            this.allelicSeqCov.setValue(allelicSeqCov);
        }

        if (this.allelicSeqCovRawVar.isEstimated() && allelicSeqCovRawVar < this.allelicSeqCovRawVar.getValue()) {
            this.allelicSeqCovRawVar.setValue(allelicSeqCovRawVar);
        }
    } // setAllelicSeqCovAndRawVar

    @Override
    public void setSeqCovLikelihoodsNodeForUpdate(int nodeIndex) {
        currentSeqCovLikelihoodsNodeIndex[nodeIndex] = 1 - currentSeqCovLikelihoodsNodeIndex[nodeIndex];
    } // setSeqCovLikelihoodsNodeForUpdate

    @Override
    public void storeLikelihoods(int taxonIndex) {
        System.arraycopy(
                seqCovLikelihoods[currentSeqCovLikelihoodsNodeIndex[taxonIndex]][taxonIndex],
                0,
                seqCovLikelihoods[1 - currentSeqCovLikelihoodsNodeIndex[taxonIndex]][taxonIndex],
                0,
                nrOfPatterns * modeledAllelesSize
        );
    } // storeLikelihoods

    /**
     * For the purpose of debugging.
     * Only the most recent accepted allelic sequencing coverage and raw variance should be stored.
     */
    @Override
    public void storeAllelicInfo() {
        throw new IllegalArgumentException("Unsupported function.");
    } // storeAllelicInfo

    @Override
    protected void store() {
        super.store();

        System.arraycopy(
                currentSeqCovLikelihoodsNodeIndex,
                0,
                storedSeqCovLikelihoodsNodeIndex,
                0,
                nrOfTaxa
        );
    } // store

    @Override
    protected void restore() {
        super.restore();

        int[] tmp = currentSeqCovLikelihoodsNodeIndex;
        currentSeqCovLikelihoodsNodeIndex = storedSeqCovLikelihoodsNodeIndex;
        storedSeqCovLikelihoodsNodeIndex = tmp;
    } // restore

    @Override
    public void computeSeqCovLikelihood(
            final int seqCov,
            final int matrixIndex,
            final int patternIndex,
            final int taxonIndex
    ) {
        final int index = patternIndex * modeledAllelesSize;

        for (int i = 0; i < modeledAllelesSize; i++)
            seqCovLikelihoods[currentSeqCovLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index + i] = computeSeqCovLikelihoodPerAllele(
                    seqCov,
                    modeledAlleles[i],
                    taxonIndex
            );
    } // computeSeqCovLikelihood

    @Override
    public double getSeqCovLikelihood(
            final int taxonIndex,
            final int index
    ) {
        return seqCovLikelihoods[currentSeqCovLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index];
    } // getSeqCovLikelihood

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

    /**
     * Compute sequencing coverage likelihood per allele under negative-binomial distribution.
     *
     * @param seqCov       sequencing coverage for current taxa under current pattern
     * @param seqAlleleNum number of sequenced alleles
     * @param taxonIndex   index of taxa
     */
    private double computeSeqCovLikelihoodPerAllele(
            final int seqCov,
            final int seqAlleleNum,
            final int taxonIndex
    ) {
        final double mean = (seqAlleleNum + numericalStabilizer) * allelicSeqCov.getValue() * sizeFactors[taxonIndex];
        final double var = mean + Math.pow(sizeFactors[taxonIndex], 2) * Math.pow(seqAlleleNum + numericalStabilizer, 2) * allelicSeqCovRawVar.getValue();
        final double p = mean / var;
        final double r = Math.pow(mean, 2) / (var - mean);

        return useLogPartials ? seqCovDensity.logDensity(seqCov, r, p) : seqCovDensity.density(seqCov, r, p);
    } // computeSeqCovLikelihoodPerAllele

}
