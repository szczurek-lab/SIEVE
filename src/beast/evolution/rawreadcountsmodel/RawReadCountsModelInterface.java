package beast.evolution.rawreadcountsmodel;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.FilteredScsAlignment;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.rawreadcountsmodel.nucreadcountsmodel.NucReadCountsModelInterface;
import beast.evolution.rawreadcountsmodel.seqcovmodel.SeqCovModelInterface;
import beast.evolution.tree.Node;
import beast.math.distributions.ScsParametricDistribution;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@Description("A raw read counts model applied to the leaves on a tree due to many kinds of uncertainties raised " +
        "during single-cell DNA sequencing, e.g., allelic dropout, sequencing error, amplification error, etc.")
public interface RawReadCountsModelInterface {

    abstract class Base extends CalculationNode implements RawReadCountsModelInterface {


        //**********************************************
        //*                   Inputs                   *
        //**********************************************

        public final Input<SeqCovModelInterface.Base> seqCovModelInput = new Input<>(
                "seqCovModel",
                "model of sequencing coverage",
                Input.Validate.REQUIRED
        );

        public final Input<NucReadCountsModelInterface.Base> nucReadCountsModelInput = new Input<>(
                "nucReadCountsModel",
                "model of nucleotide read counts",
                Input.Validate.REQUIRED
        );

        public final Input<RealParameter> adoRateInput = new Input<>(
                "adoRate",
                "allelic dropout rate per site",
                Input.Validate.REQUIRED
        );

        public final Input<Boolean> singleADOInput = new Input<>(
                "singleADO",
                "single ADO mode; the default is single ADO (true); set false to have locus dropout (under experiment)",
                Input.Validate.OPTIONAL
        );


        //***********************************************
        //*                  Variables                  *
        //***********************************************

        // Should be ONLY assigned through `deeplyInitialize` or `duplicate`.
        // Should be the SAME object (a reference) as the one owes by `ScsLikelihood` instance.
        protected ScsAlignment alignment;

        protected SeqCovModelInterface.Base seqCovModel;
        protected NucReadCountsModelInterface.Base nucReadCountsModel;

        protected RealParameter adoRate;

        /**
         * Number of sequenced alleles taking into account both model of evolution and ADO.
         * Used frequently.
         * Should be initialized from ScsSubstitutionModel.Base.
         * Must be sorted when initializing.
         */
        protected int modeledAllelesSize;
        protected int[] modeledAlleles;

        /**
         * true: default, locus dropout
         * false: single ADO
         */
        protected boolean singleADO;

        protected int nrOfMatrices;
        protected int nrOfPatterns;
        protected int nrOfTaxa;
        protected int nrOfStates;
        protected int nrOfCategories;

        protected boolean useLogPartials;

        /**
         * a flag to indicate whether this instance has been deeply initialized or not
         */
        protected boolean deeplyInitialized = false;

        /**
         * a flag to indicate whether the caller is in debug mode or not
         */
        protected boolean inDebugMode = false;


        //**********************************************
        //*              Abstract methods              *
        //**********************************************

        /**
         * Compute the mixed likelihood given the sequencing coverage likelihood and nucleotide read counts likelihood.
         *
         * @param matrixIndex   apparently
         * @param patternIndex  apparently
         * @param taxonIndex    apparently
         * @param genotypeIndex apparently
         * @return the mixed likelihood
         */
        public abstract double computeMixedLikelihood(
                final int matrixIndex,
                final int patternIndex,
                final int taxonIndex,
                final int genotypeIndex
        );

        /**
         * Compute the mixed likelihood given the sequencing coverage likelihood and nucleotide read counts likelihood.
         *
         * @param matrixIndex   apparently
         * @param patternIndex  apparently
         * @param taxonIndex    apparently
         * @param genotypeIndex apparently
         * @param MLCompIndex   index to the maximum likelihood component, which corresponding to the number of sequenced alleles
         * @return the mixed likelihood
         */
        public abstract double computeMixedLikelihood(
                final int matrixIndex,
                final int patternIndex,
                final int taxonIndex,
                final int genotypeIndex,
                List<Integer>[] MLCompIndex
        );

        /**
         * Get the log likelihoods for ADOs.
         *
         * @param matrixIndex   apparently
         * @param patternIndex  apparently
         * @param taxonIndex    apparently
         * @param genotypeIndex apparently
         * @return log likelihoods of ADO for {@param genotypeIndex}
         */
        public abstract double[] getAdoLogLikelihoods(
                final int matrixIndex,
                final int patternIndex,
                final int taxonIndex,
                final int genotypeIndex
        );


        //**************************************************
        //*               Overridden methods               *
        //**************************************************

        @Override
        public void initAndValidate() {
            seqCovModel = seqCovModelInput.get();
            nucReadCountsModel = nucReadCountsModelInput.get();

            if (adoRateInput.get().getValue() < 0.0 ||
                    (adoRateInput.get().isEstimated() &&
                            (adoRateInput.get().getLower() < 0.0 ||
                                    adoRateInput.get().getUpper() < 0.0 ||
                                    adoRateInput.get().getUpper() > 1.0)))
                throw new IllegalArgumentException("adoRate and its bounds should be defined as no smaller than 0 (" +
                        this.getClass().getName() + ")");
            else
                adoRate = adoRateInput.get();

            if (singleADOInput.get() == null)
                singleADO = true;
            else
                singleADO = singleADOInput.get();
        } // initAndValidate


        //***********************************************
        //*                   Methods                   *
        //***********************************************

        /**
         * Whether the number of sequenced alleles is varying or not due to the existence of ADO.
         *
         * @return varying or not
         */
        private boolean isVaryingSeqAlleles() {
            return adoRate.isEstimatedInput.get() ||
                    (adoRateInput.get().getValue() > 0.0 && adoRateInput.get().getValue() < 1.0) ||
                    adoRateInput.get().getValue() == 1.0;
        } // isVaryingSeqAlleles

        /**
         * Initialize parameters with passed arguments from ScsTreeLikelihood.
         *
         * @param alignment      apparently
         * @param matrixNum      number of matrices
         * @param stateNum       number of states
         * @param modeledAlleles alleles the model of evolution models
         * @param useLogPartials apparently
         */
        public void deeplyInitialize(
                final ScsAlignment alignment,
                final int matrixNum,
                final int stateNum,
                final int[] modeledAlleles,
                final boolean useLogPartials
        ) {
            this.alignment = alignment;
            nrOfPatterns = alignment.getPatternCount();
            nrOfTaxa = alignment.getTaxonCount();
            nrOfStates = stateNum;
            this.useLogPartials = useLogPartials;

            // initialize the modeled alleles by taking into account both the model of evolution and ADO
            Set<Integer> tmp = Arrays.stream(modeledAlleles).boxed().collect(Collectors.toSet());
            if (isVaryingSeqAlleles()) {
                for (int i : modeledAlleles) {
                    if (i - 1 >= 0)
                        tmp.add(i - 1);
                    if (!singleADO && i - 2 >= 0)
                        tmp.add(i - 2);
                }
            }
            modeledAllelesSize = tmp.size();
            this.modeledAlleles = tmp.stream().mapToInt(Integer::intValue).sorted().toArray();

            seqCovModel.deeplyInitialize(
                    alignment,
                    matrixNum,
                    this.modeledAlleles,
                    useLogPartials
            );
            nucReadCountsModel.deeplyInitialize(
                    alignment,
                    stateNum,
                    useLogPartials
            );

            nrOfMatrices = seqCovModel.getNrOfMatrices();
            nrOfCategories = nucReadCountsModel.getNrOfCategories();

            deeplyInitialized = true;
        } // deeplyInitialize

        /**
         * Duplicate some important variables entirely or within range patterns.
         * Used when there are more than one threads.
         *
         * @param target     target instance storing variables
         * @param targetData data corresponding to the target raw read counts model
         */
        public void duplicate(
                RawReadCountsModelInterface.Base target,
                final ScsAlignment targetData
        ) {
            assert targetData instanceof FilteredScsAlignment;
            target.alignment = targetData;

            target.nrOfMatrices = this.nrOfMatrices;
            target.nrOfPatterns = targetData.getPatternCount();
            target.nrOfTaxa = this.nrOfTaxa;
            target.nrOfStates = this.nrOfStates;
            target.nrOfCategories = this.nrOfCategories;

            target.modeledAllelesSize = this.modeledAllelesSize;
            target.modeledAlleles = this.modeledAlleles;
            target.useLogPartials = this.useLogPartials;

            this.seqCovModel.duplicate(
                    target.seqCovModelInput.get(),
                    targetData
            );
            this.nucReadCountsModel.duplicate(
                    target.nucReadCountsModelInput.get(),
                    targetData
            );

            target.setDeeplyInitialized(true);
        } // duplicate

        /**
         * @param taxon the taxon name as a string
         * @param data  the alignment
         * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
         * or -1 if the taxon is not in the alignment.
         */
        private int getTaxonIndex(String taxon, ScsAlignment data) {
            int taxonIndex = data.getTaxonIndex(taxon);
            if (taxonIndex == -1) {
                if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                    taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
                }
                if (taxonIndex == -1) {
                    throw new RuntimeException("Could not find sequence " + taxon + " in the ScsAlignment (" +
                            this.getClass().getName() + ")");
                }
            }
            return taxonIndex;
        } // getTaxonIndex

        /**
         * Initialize the leaf likelihoods.
         * This is MANDATORY during the initialization of tree likelihoods.
         *
         * @param node leaf node
         * @return the likelihoods for {@param node} across all patterns but only for one matrix.
         */
        public final double[] initializeLeafLikelihood(final Node node) {
            assert node.isLeaf();
            final int taxonIndex = getTaxonIndex(node.getID(), alignment);

            // For only one matrix because the likelihoods for other matrices are the same during initialization.
            double[] partials = new double[nrOfPatterns * nrOfStates];

            int index = 0;
            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                final int[] reads = alignment.getPattern(taxonIndex, patternIndex);

                // Compute likelihoods of sequencing coverage.
                seqCovModel.computeSeqCovLikelihood(
                        alignment.getSequencingCoverage(reads),
                        0,
                        patternIndex,
                        taxonIndex
                );

                // Compute likelihoods of nucleotide read counts.
                nucReadCountsModel.computeNucReadCountsLikelihoods(
                        reads,
                        patternIndex,
                        taxonIndex
                );

                for (int genotypeIndex = 0; genotypeIndex < nrOfStates; genotypeIndex++) {
                    partials[index] = computeMixedLikelihood(
                            0,
                            patternIndex,
                            taxonIndex,
                            genotypeIndex
                    );

                    index++;
                }
            }

            seqCovModel.expandVariables(taxonIndex);

            storeLikelihoods(taxonIndex);

            return partials;
        } // initializeLeafLikelihood

        /**
         * Initialize the leaf likelihoods.
         * This is MANDATORY during the initialization of tree likelihoods.
         *
         * @param node             leaf node
         * @param MLNrOfSeqAlleles a data structure to store the maximum likelihood number of sequenced alleles
         * @return the likelihoods for {@param node} across all patterns but only for one matrix.
         */
        public final double[] initializeLeafLikelihood(
                final Node node,
                List<Integer>[][][][] MLNrOfSeqAlleles
        ) {
            assert node.isLeaf();
            final int nodeIndex = node.getNr();
            final int taxonIndex = getTaxonIndex(node.getID(), alignment);

            // If necessary, initialize `MLNrOfSeqAlleles` for each matrix.
            final boolean initializeOtherMatrices = MLNrOfSeqAlleles[0][nodeIndex].length > nrOfPatterns * nrOfStates;

            // For only one matrix because the likelihoods for other matrices are the same during initialization.
            double[] partials = new double[nrOfPatterns * nrOfStates];

            int index = 0;
            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                final int[] reads = alignment.getPattern(taxonIndex, patternIndex);

                // Compute likelihoods of sequencing coverage.
                seqCovModel.computeSeqCovLikelihood(
                        alignment.getSequencingCoverage(reads),
                        0,
                        patternIndex,
                        taxonIndex
                );

                // Compute likelihoods of nucleotide read counts.
                nucReadCountsModel.computeNucReadCountsLikelihoods(
                        reads,
                        patternIndex,
                        taxonIndex
                );

                for (int genotypeIndex = 0; genotypeIndex < nrOfStates; genotypeIndex++) {
                    partials[index] = computeMixedLikelihood(
                            0,
                            patternIndex,
                            taxonIndex,
                            genotypeIndex,
                            MLNrOfSeqAlleles[0][nodeIndex][index]
                    );

                    if (MLNrOfSeqAlleles[1][nodeIndex][index][0] == null)
                        MLNrOfSeqAlleles[1][nodeIndex][index][0] = new ArrayList<>();
                    else
                        MLNrOfSeqAlleles[1][nodeIndex][index][0].clear();
                    for (Integer i : MLNrOfSeqAlleles[0][nodeIndex][index][0])
                        MLNrOfSeqAlleles[1][nodeIndex][index][0].add(new Integer(i.intValue()));

                    if (initializeOtherMatrices) {
                        int matrixIndex = 1;

                        while (index + matrixIndex * nrOfPatterns * nrOfStates < MLNrOfSeqAlleles[0][nodeIndex].length) {
                            if (MLNrOfSeqAlleles[0][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0] == null)
                                MLNrOfSeqAlleles[0][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0] = new ArrayList<>();
                            else
                                MLNrOfSeqAlleles[0][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0].clear();

                            if (MLNrOfSeqAlleles[1][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0] == null)
                                MLNrOfSeqAlleles[1][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0] = new ArrayList<>();
                            else
                                MLNrOfSeqAlleles[1][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0].clear();

                            for (Integer i : MLNrOfSeqAlleles[0][nodeIndex][index][0]) {
                                MLNrOfSeqAlleles[0][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0].add(new Integer(i.intValue()));
                                MLNrOfSeqAlleles[1][nodeIndex][index + matrixIndex * nrOfPatterns * nrOfStates][0].add(new Integer(i.intValue()));
                            }

                            matrixIndex++;
                        }
                    }

                    index++;
                }
            }

            seqCovModel.expandVariables(taxonIndex);

            storeLikelihoods(taxonIndex);

            return partials;
        } // initializeLeafLikelihood

        /**
         * Compute the leaf likelihoods during MCMC, including routine debugging.
         *
         * @param node leaf node
         * @return the likelihoods for {@param node} across all patterns and matrices.
         */
        public final double[] computeLeafLikelihood(final Node node) {
            assert node.isLeaf();
            final int taxonIndex = getTaxonIndex(node.getID(), alignment);

            if (seqCovModel.isDirtyCalculation())
                seqCovModel.setSeqCovLikelihoodsNodeForUpdate(taxonIndex);

            if (nucReadCountsModel.isDirtyCalculation())
                nucReadCountsModel.setNucReadCountsLikelihoodsNodeForUpdate(taxonIndex);

            double[] partials = new double[nrOfMatrices * nrOfPatterns * nrOfStates];

            int index = 0;
            for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {
                for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                    final int[] reads = alignment.getPattern(taxonIndex, patternIndex);

                    // Compute likelihoods of sequencing coverage if necessary.
                    if (seqCovModel.isDirtyCalculation() || (inDebugMode && seqCovModel.updateSeqCovModel()))
                        seqCovModel.computeSeqCovLikelihood(
                                alignment.getSequencingCoverage(reads),
                                matrixIndex,
                                patternIndex,
                                taxonIndex
                        );

                    // Compute likelihoods of nucleotide read counts if necessary.
                    if (matrixIndex == 0 && nucReadCountsModel.isDirtyCalculation())
                        nucReadCountsModel.computeNucReadCountsLikelihoods(
                                reads,
                                patternIndex,
                                taxonIndex
                        );

                    for (int genotypeIndex = 0; genotypeIndex < nrOfStates; genotypeIndex++) {
                        partials[index] = computeMixedLikelihood(
                                matrixIndex,
                                patternIndex,
                                taxonIndex,
                                genotypeIndex
                        );

                        index++;
                    }
                }
            }

            return partials;
        } // computeLeafLikelihood

        /**
         * Compute the leaf likelihoods during MCMC, including routine debugging.
         *
         * @param node             leaf node
         * @param MLNrOfSeqAlleles a data structure to store the maximum likelihood number of sequenced alleles
         * @return the likelihoods for {@param node} across all patterns and matrices.
         */
        public final double[] computeLeafLikelihood(
                final Node node,
                List<Integer>[][] MLNrOfSeqAlleles
        ) {
            assert node.isLeaf();
            final int taxonIndex = getTaxonIndex(node.getID(), alignment);

            if (seqCovModel.isDirtyCalculation())
                seqCovModel.setSeqCovLikelihoodsNodeForUpdate(taxonIndex);

            if (nucReadCountsModel.isDirtyCalculation())
                nucReadCountsModel.setNucReadCountsLikelihoodsNodeForUpdate(taxonIndex);

            double[] partials = new double[nrOfMatrices * nrOfPatterns * nrOfStates];

            int index = 0;
            for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {
                for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                    final int[] reads = alignment.getPattern(taxonIndex, patternIndex);

                    // Compute likelihoods of sequencing coverage if necessary.
                    if (seqCovModel.isDirtyCalculation() || (inDebugMode && seqCovModel.updateSeqCovModel()))
                        seqCovModel.computeSeqCovLikelihood(
                                alignment.getSequencingCoverage(reads),
                                matrixIndex,
                                patternIndex,
                                taxonIndex
                        );

                    // Compute likelihoods of nucleotide read counts if necessary.
                    if (matrixIndex == 0 && nucReadCountsModel.isDirtyCalculation())
                        nucReadCountsModel.computeNucReadCountsLikelihoods(
                                reads,
                                patternIndex,
                                taxonIndex
                        );

                    for (int genotypeIndex = 0; genotypeIndex < nrOfStates; genotypeIndex++) {
                        partials[index] = computeMixedLikelihood(
                                matrixIndex,
                                patternIndex,
                                taxonIndex,
                                genotypeIndex,
                                MLNrOfSeqAlleles[index]
                        );

                        index++;
                    }
                }
            }

            return partials;
        } // computeLeafLikelihood

        /**
         * Update leaf likelihoods for changed patterns.
         * Should ONLY be called during post processing.
         *
         * @param node             leaf node
         * @param changedPatterns  changed patterns (matrixIndex, patternIndex)
         * @param MLNrOfSeqAlleles a data structure to store the maximum likelihood number of sequenced alleles
         * @param partials         partial likelihoods to be updated
         */
        public void updatePartialLeafLikelihoods(
                final Node node,
                final int[][] changedPatterns,
                List<Integer>[][] MLNrOfSeqAlleles,
                double[] partials
        ) {
            assert node.isLeaf();
            final int taxonIndex = getTaxonIndex(node.getID(), alignment);

            for (int[] pair : changedPatterns) {
                // (matrixIndex, patternIndex)
                final int index = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

                final int[] reads = alignment.getPattern(taxonIndex, pair[1]);

                seqCovModel.computeSeqCovLikelihood(
                        alignment.getSequencingCoverage(reads),
                        pair[0],
                        pair[1],
                        taxonIndex
                );

                if (pair[0] == 0)
                    nucReadCountsModel.computeNucReadCountsLikelihoods(
                            reads,
                            pair[1],
                            taxonIndex
                    );

                for (int genotypeIndex = 0; genotypeIndex < nrOfStates; genotypeIndex++)
                    partials[index + genotypeIndex] = computeMixedLikelihood(
                            pair[0],
                            pair[1],
                            taxonIndex,
                            genotypeIndex,
                            MLNrOfSeqAlleles[index]
                    );
            }
        } // updatePartialLeafLikelihoods

        /**
         * If the some parameters in the sequencing coverage model are updated in post processing stage
         * but not during MCMC, the corresponding data structures should be stored for the sake of debugging.
         */
        public void storeSeqCovInfo() {
            seqCovModel.storeAllelicInfo();
        } // storeSeqCovInfo

        /**
         * For initialization.
         *
         * @param taxonIndex apparently
         */
        protected void storeLikelihoods(int taxonIndex) {
            seqCovModel.storeLikelihoods(taxonIndex);
            nucReadCountsModel.storeLikelihoods(taxonIndex);
        } // setIndexForUpdate

        /**
         * Compute allelic sequencing coverage and its raw variance.
         * Should be called at the post processing stage of computing tree likelihood.
         *
         * @param changedPatterns matrices and patterns which should be updated: (matrixIndex, patternIndex)
         */
        public void updateAllelicSeqCovAndRawVar(List<List<Integer>> changedPatterns) {
            seqCovModel.updateAllelicSeqCovAndRawVar(changedPatterns);
        } // updateAllelicSeqCovAndRawVar

        /**
         * Get the parameters situation 1 where the wild type (0/0) applies.
         *
         * @param out communicated arguments
         */
        public void getWildTypeNucReadCountModelParams(double[] out) {
            nucReadCountsModel.getWildTypeNucReadCountModelParams(out);
        } // getWildTypeNucReadCountModelParams


        //***********************************************
        //*              Getter and Setter              *
        //***********************************************

        public int getModeledAllelesSize() {
            return modeledAllelesSize;
        } // getNrOfAdoStates

        public int getNrOfParamsInReadCountsModel() {
            return nucReadCountsModel.getNrOfParamsInReadCountsModel();
        } // getNrOfParamsInReadCountsModel

        /**
         * check whether critical variables in raw read counts model have been computed or not
         * used in threaded likelihood computation to avoid wasting computer resources
         *
         * @return true or false
         */
        public boolean isDeeplyInitialized() {
            return deeplyInitialized;
        } // isDeeplyInitialized

        /**
         * set the deeplyInitialized flag
         *
         * @param flag initialized (true) or not (false)
         */
        public void setDeeplyInitialized(boolean flag) {
            deeplyInitialized = flag;
        } // setDeeplyInitialized

        public boolean isInDebugMode() {
            return inDebugMode;
        } // isInDebugMode

        public void setInDebugMode(boolean inDebugMode) {
            this.inDebugMode = inDebugMode;
            this.seqCovModel.setInDebugMode(inDebugMode);
        } // setInDebugMode

        public boolean updateSeqCovModel() {
            return seqCovModel.updateSeqCovModel();
        } // updateSeqCovModel

        public boolean isInVariantCallingMode() {
            return seqCovModel.isInVariantCallingMode();
        } // isInVariantCallingMode

        /**
         * add patterns to MLNrOfSequencedAllelesPatterns[index]
         *
         * @param index  which matrix and pattern
         * @param values apparently
         */
        public void addMLNrOfExistingAllelesPatterns(final int index, final int[] values) {
            seqCovModel.addMLNrOfExistingAllelesPatterns(index, values);
        } // addMLNrOfExistingAllelesPatterns

        /**
         * add weight to MLNrOfSequencedAllelesWeights[index]
         *
         * @param index which matrix and pattern
         * @param value apparently
         */
        public void addMLNrOfSequencedAllelesWeights(final int index, final int value) {
            seqCovModel.addMLNrOfSequencedAllelesWeights(index, value);
        } // addMLNrOfSequencedAllelesWeights

        public List<int[]> getMLNrOfSequencedAllelesPatterns(final int index) {
            return seqCovModel.getMLNrOfSequencedAllelesPatterns(index);
        } // getMLNrOfSequencedAllelesPatterns

        public int getMLNrOfSequencedAllelesWeights(final int index1, final int index2) {
            return seqCovModel.getMLNrOfSequencedAllelesWeights(index1, index2);
        } // getMLNrOfSequencedAllelesWeights

        public void resetMLNrOfSequencedAllelesPatternsAndWeights(final int index) {
            seqCovModel.resetMLNrOfSequencedAllelesPatternsAndWeights(index);
        } // resetMLNrOfSequencedAllelesPatternsAndWeights

        public ScsParametricDistribution getNucReadCountsModelDistribution() {
            return nucReadCountsModel.getDistribution();
        } // getNucReadCountsModelDistribution


        //******************************************
        //*                 Logger                 *
        //******************************************

        /**
         * log sampled allelic sequencing coverage and raw variance
         *
         * @param index     which pattern is about to be logged?
         * @param out       apparently
         */
        public void logCovar(int index, PrintStream out) {
            seqCovModel.logCovar(index, out);
        } // logCovar

    } // Base

} // ScsErrorModel
