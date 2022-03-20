package beast.evolution.rawreadcountsmodel.seqcovmodel;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.ScsAlignment;
import beast.math.distributions.NegativeBinomial;
import beast.math.statistic.DiscreteStatistics;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static beast.math.statistic.DiscreteStatistics.median;
import static beast.math.util.MathFunctions.geometricMean;

public interface SeqCovModelInterface {

    abstract class Base extends CalculationNode implements SeqCovModelInterface {


        //**********************************************
        //*                   Inputs                   *
        //**********************************************

        public final Input<Integer> zeroCovModeInput = new Input<>(
                "zeroCovMode",
                "the way to deal with zero coverage; the default is omitting zeros (0); or increasing by 1 (1)",
                Input.Validate.OPTIONAL
        );

        public final Input<RealParameter> allelicSeqCovInput = new Input<>(
                "allelicSeqCov",
                "allelic sequencing coverage"
        );

        public final Input<String> allelicSeqCovArrayInput = new Input<>(
                "allelicSeqCovArray",
                "a single or an array of allelic sequencing coverage (comma separated for a matrix, and semicolon separated among matrices)"
        );//, Input.Validate.XOR, allelicSeqCovInput);

        public final Input<RealParameter> allelicSeqCovRawVarInput = new Input<>(
                "allelicSeqCovRawVar",
                "allelic sequencing coverage raw variance"
        );

        public final Input<String> allelicSeqCovRawVarArrayInput = new Input<>(
                "allelicSeqCovRawVarArray",
                "a single or an array of allelic sequencing coverage raw variance (comma separated for a matrix, and semicolon separated among matrices)"
        );//, Input.Validate.XOR, allelicSeqCovRawVarInput);

        public final Input<Boolean> inVariantCallingModeInput = new Input<>(
                "inVariantCallingMode",
                "whether is in variant calling mode or not",
                Input.Validate.REQUIRED
        );

        public final Input<Boolean> initSharedAllelicSeqCovAndRawVarInput = new Input<>(
                "initSharedAllelicSeqCovAndRawVar",
                "whether to use the statistics of allelic sequencing coverage and raw variance from samples to " +
                        "initialize the corresponding variables in the ExploredSharedAllelicSeqCovModel",
                Input.Validate.OPTIONAL
        );


        //***********************************************
        //*                  Variables                  *
        //***********************************************

        protected int nrOfMatrices;
        protected int nrOfPatterns;
        protected int nrOfTaxa;

        protected boolean useLogPartials;

        /**
         * Number of sequenced alleles taking into account both model of evolution and ADO.
         * Used frequently.
         * Should be initialized from ScsSubstitutionModel.Base.
         * Must be sorted when initializing.
         */
        protected int modeledAllelesSize;
        protected int[] modeledAlleles;

        protected RealParameter allelicSeqCov;
        protected RealParameter allelicSeqCovRawVar;

        protected int zeroCovMode;

        protected double numericalStabilizer;

        // #taxa
        protected double[] sizeFactors;

        /**
         * Store the processed sequencing coverage.
         * Either ignore the zeros or add all entries by one.
         * #patterns * #taxa
         */
        protected int[] processedSeqCov;

        /**
         * Used to update allelic sequencing coverage and its raw variance in the raw read counts model.
         * Considering the weights for all possible maximum likelihood number of sequenced alleles.
         */
        protected List<int[]>[] MLNrOfSequencedAllelesPatterns; // only for leaves, #matrices * #patterns
        protected List<Integer>[] MLNrOfSequencedAllelesWeights; // only for leaves, #matrices * #patterns

        /**
         * Sequencing coverage likelihood.
         * Different because of number of sequenced alleles (in ascending order).
         * [2 (or 1)] * [#taxa] * [(#matrices) * #patterns * modeledAllelesSize]
         */
        protected double[][][] seqCovLikelihoods;

        protected int[] currentSeqCovLikelihoodsNodeIndex;
        protected int[] storedSeqCovLikelihoodsNodeIndex;

        protected NegativeBinomial seqCovDensity;

        /**
         * a flag to indicate whether allelic sequencing coverage and raw variance should be updated during MCMC
         */
        protected boolean needToUpdate;

        /**
         * whether in variant calling mode or not
         */
        protected boolean inVariantCallingMode = false;

        /**
         * a flag to indicate whether the caller is in debug mode or not
         */
        protected boolean inDebugMode = false;


        //**********************************************
        //*              Abstract methods              *
        //**********************************************

        /**
         * For the purpose of debugging.
         * Only the most recent accepted allelic sequencing coverage and raw variance should be stored.
         */
        public abstract void storeAllelicInfo();

        public abstract void setSeqCovLikelihoodsNodeForUpdate(int nodeIndex);

        /**
         * During likelihood initialization, likelihoods should be copied for the sake of accurate tree likelihood computation.
         *
         * @param taxonIndex apparently
         */
        public abstract void storeLikelihoods(int taxonIndex);

        /**
         * Called during initialization of tree likelihoods for internal variables expanding / copying to other matrices.
         *
         * @param taxonIndex apparently
         */
        public abstract void expandVariables(int taxonIndex);

        /**
         * Compute sequencing coverage likelihood under negative-binomial distribution.
         *
         * @param seqCov       sequencing coverage for current taxa under current pattern
         * @param matrixIndex  index to matrix
         * @param patternIndex index to pattern
         * @param taxonIndex   index to taxa
         */
        public abstract void computeSeqCovLikelihood(
                final int seqCov,
                final int matrixIndex,
                final int patternIndex,
                final int taxonIndex
        );

        /**
         * Get the sequencing coverage likelihood at {@param index}.
         *
         * @param taxonIndex apparently
         * @param index      index to the last dimension
         * @return sequencing coverage likelihood
         */
        public abstract double getSeqCovLikelihood(
                final int taxonIndex,
                final int index
        );

        /**
         * Compute allelic sequencing coverage and its raw variance.
         * Should be called at the post processing stage of computing tree likelihood.
         *
         * @param changedPatterns matrices and patterns which should be updated: (matrixIndex, patternIndex)
         */
        public abstract void updateAllelicSeqCovAndRawVar(List<List<Integer>> changedPatterns);


        //**********************************************
        //*             Overridden methods             *
        //**********************************************

        @Override
        public void initAndValidate() {
            inVariantCallingMode = inVariantCallingModeInput.get();

            if (zeroCovModeInput.get() == null)
                zeroCovMode = 0;
            else if (zeroCovModeInput.get() >= 0 && zeroCovModeInput.get() <= 1)
                zeroCovMode = zeroCovModeInput.get();
            else
                throw new IllegalArgumentException("ERROR! Illegal value for zeroCovMode. Only 0 or 1 is supported.");

            numericalStabilizer = 1e-6;

            seqCovDensity = new NegativeBinomial();
        } // initAndValidate


        //*********************************************
        //*                  Methods                  *
        //*********************************************

        /**
         * Initialize parameters with passed arguments from ScsTreeLikelihood.
         * Must be called by every inheriting class at the beginning.
         *
         * @param alignment      apparently
         * @param matrixNum      number of matrices
         * @param modeledAlleles alleles the model of evolution models
         * @param useLogPartials apparently
         */
        public void deeplyInitialize(
                final ScsAlignment alignment,
                final int matrixNum,
                final int[] modeledAlleles,
                final boolean useLogPartials
        ) {
            nrOfPatterns = alignment.getPatternCount();
            nrOfTaxa = alignment.getTaxonCount();
            nrOfMatrices = matrixNum;
            this.useLogPartials = useLogPartials;

            modeledAllelesSize = modeledAlleles.length;
            this.modeledAlleles = new int[modeledAllelesSize];
            System.arraycopy(modeledAlleles, 0, this.modeledAlleles, 0, modeledAllelesSize);
        } // deeplyInitialize

        /**
         * Duplicate some important variables entirely or within ranged patterns.
         * Used when there are more than one thread.
         * Must be called by every inheriting class at the beginning.
         *
         * @param target     target instance storing variables
         * @param targetData data corresponding to the target raw read counts model
         */
        public void duplicate(
                SeqCovModelInterface.Base target,
                final ScsAlignment targetData
        ) {
            target.nrOfMatrices = this.nrOfMatrices;
            target.nrOfPatterns = targetData.getPatternCount();
            target.nrOfTaxa = this.nrOfTaxa;

            target.needToUpdate = this.needToUpdate;
            target.useLogPartials = this.useLogPartials;

            target.modeledAllelesSize = this.modeledAllelesSize;
            target.modeledAlleles = this.modeledAlleles;
        } // duplicate

        /**
         * Set allelic sequencing coverage and raw variance with statistics from samples for a good starting point.
         * Only valid for ExploredSharedAllelicSeqCovModel.
         *
         * @param allelicSeqCov       apparently
         * @param allelicSeqCovRawVar apparently
         */
        protected void setAllelicSeqCovAndRawVar(double allelicSeqCov, double allelicSeqCovRawVar) {} // setAllelicSeqCovAndRawVar

        /**
         * Compute size factors.
         * Should only be called once through `deeplyInitialize()`.
         *
         * @param alignment apparently
         */
        protected void computeSizeFactors(final ScsAlignment alignment) {
            double[] geoMeanCovPerPattern = new double[nrOfPatterns];

            // process zero coverage and store values
            int id = 0;
            int zeroCovCounter = 0;
            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                final int patternWeight = alignment.getPatternWeight(patternIndex);

                for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
                    final int cov = alignment.getSequencingCoverageOfPattern(taxonIndex, patternIndex);

                    if (cov == 0) {
                        zeroCovCounter += patternWeight;
                    }

                    if (zeroCovMode == 0) {
                        // ignore zero coverage
                        processedSeqCov[id + taxonIndex] = cov;
                    } else if (zeroCovMode == 1) {
                        // increase all entries by 1
                        processedSeqCov[id + taxonIndex] = cov + 1;
                    }
                }

                // compute geometric means
                geoMeanCovPerPattern[patternIndex] = geometricMean(processedSeqCov, id, id + nrOfTaxa);

                id += nrOfTaxa;
            }

            final double statTmp1 = DiscreteStatistics.mean(Doubles.toArray(Ints.asList(processedSeqCov))) / 2;
            final double statTmp2 = DiscreteStatistics.variance(Doubles.toArray(Ints.asList(processedSeqCov))) / 4;

            setAllelicSeqCovAndRawVar(statTmp1, statTmp2);

            Log.info.println("Percentage of missing data: " + String.format("%.1f", ((double) zeroCovCounter * 100.0) / (alignment.getLociNr() * nrOfTaxa)) + "%");
            Log.info.println("Mean of allelic sequencing coverage: " + String.format("%.1f", statTmp1));
            Log.info.println("Mean of allelic sequencing coverage raw variance: " + String.format("%.1f", statTmp2));

            // compute size factors
            double[] intermediateArr = new double[alignment.getLociNr()];
            for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
                int index = 0;
                for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                    if (zeroCovMode == 1 ||
                            (processedSeqCov[patternIndex * nrOfTaxa + taxonIndex] > 0 && geoMeanCovPerPattern[patternIndex] > 0)) {
                        Arrays.fill(
                                intermediateArr,
                                index,
                                index + alignment.getPatternWeight(patternIndex),
                                processedSeqCov[patternIndex * nrOfTaxa + taxonIndex] / geoMeanCovPerPattern[patternIndex]
                        );

                        index += alignment.getPatternWeight(patternIndex);
                    }
                }

                sizeFactors[taxonIndex] = median(Arrays.copyOfRange(intermediateArr, 0, index));
            }
        } // computeSizeFactors

        public final double[] getSizeFactors() {
            return sizeFactors;
        } // getSizeFactors

        public int getNrOfMatrices() {
            return nrOfMatrices;
        } // getNrOfMatrices

        public boolean updateSeqCovModel() {
            return needToUpdate;
        } // updateSeqCovModel

        public boolean isInVariantCallingMode() {
            return inVariantCallingMode;
        } // isInVariantCallingMode

        public void setInDebugMode(boolean inDebugMode) {
            this.inDebugMode = inDebugMode;
        } // setInDebugMode

        /**
         * add patterns to MLNrOfSequencedAllelesPatterns[index]
         *
         * @param index  which matrix and pattern
         * @param values apparently
         */
        public void addMLNrOfExistingAllelesPatterns(final int index, final int[] values) {
            if (this.MLNrOfSequencedAllelesPatterns[index] == null)
                this.MLNrOfSequencedAllelesPatterns[index] = new ArrayList<>();

            int[] tmp = new int[values.length];
            System.arraycopy(values, 0, tmp, 0, values.length);

            this.MLNrOfSequencedAllelesPatterns[index].add(tmp);
        } // addMLNrOfExistingAllelesPatterns

        /**
         * add weight to MLNrOfSequencedAllelesWeights[index]
         *
         * @param index which matrix and pattern
         * @param value apparently
         */
        public void addMLNrOfSequencedAllelesWeights(final int index, final int value) {
            if (this.MLNrOfSequencedAllelesWeights[index] == null)
                this.MLNrOfSequencedAllelesWeights[index] = new ArrayList<>();

            this.MLNrOfSequencedAllelesWeights[index].add(value);
        } // addMLNrOfSequencedAllelesWeights

        public List<int[]> getMLNrOfSequencedAllelesPatterns(final int index) {
            return this.MLNrOfSequencedAllelesPatterns[index];
        } // getMLNrOfSequencedAllelesPatterns

        public int getMLNrOfSequencedAllelesWeights(final int index1, final int index2) {
            return this.MLNrOfSequencedAllelesWeights[index1].get(index2);
        } // getMLNrOfSequencedAllelesWeights

        public void resetMLNrOfSequencedAllelesPatternsAndWeights(final int index) {
            if (this.MLNrOfSequencedAllelesPatterns[index] != null)
                this.MLNrOfSequencedAllelesPatterns[index].clear();

            if (this.MLNrOfSequencedAllelesWeights[index] != null)
                this.MLNrOfSequencedAllelesWeights[index].clear();
        } // resetMLNrOfSequencedAllelesPatternsAndWeights


        //******************************************
        //*                 Logger                 *
        //******************************************

        /**
         * log sampled allelic sequencing coverage and raw variance
         *
         * @param index which pattern is about to be logged?
         * @param out   apparently
         */
        public abstract void logCovar(int index, PrintStream out);

    }

}
