package beast.evolution.rawreadcountsmodel.nucreadcountsmodel;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.ScsAlignment;
import beast.math.distributions.ScsParametricDistribution;

public interface NucReadCountsModelInterface {

    abstract class Base extends CalculationNode implements NucReadCountsModelInterface {


        //**********************************************
        //*                   Inputs                   *
        //**********************************************

        public final Input<RealParameter> effSeqErrRateInput = new Input<>(
                "effSeqErrRate",
                "effective sequencing error rate",
                Input.Validate.REQUIRED
        );

        public final Input<RealParameter> shapeCtrl1Input = new Input<>(
                "shapeCtrl1",
                "shape controller for homozygous genotypes",
                Input.Validate.REQUIRED
        );

        public final Input<RealParameter> shapeCtrl2Input = new Input<>(
                "shapeCtrl2",
                "shape controller for heterozygous genotypes",
                Input.Validate.REQUIRED
        );


        //***********************************************
        //*                  Variables                  *
        //***********************************************

        protected int nrOfPatterns;
        protected int nrOfTaxa;
        protected int nrOfStates;

        // The number of category of nucleotide read counts likelihood.
        protected int nrOfCategories;

        protected boolean useLogPartials;

        protected RealParameter effSeqErrRate;
        protected RealParameter shapeCtrl1;
        protected RealParameter shapeCtrl2;

        /**
         * Number of parameters in the model of nucleotide read counts,
         * e.g., for Beta-Binomial: 2; for Dirichlet-Multinomial: 4.
         */
        protected int nrOfNucReadCountsModelParam;

        /**
         * Likelihoods for nucleotide read counts model.
         * [2] * [#taxa] * [#patterns * #categories]
         */
        protected double[][][] nucReadCountsLikelihoods;

        protected int[] currentNucReadCountsLikelihoodsNodeIndex;
        protected int[] storedNucReadCountsLikelihoodsNodeIndex;

        protected ScsParametricDistribution variantDensity;


        //**********************************************
        //*              Abstract methods              *
        //**********************************************

        /**
         * Compute nucleotide read counts likelihoods.
         *
         * @param reads        reads for current taxa under current pattern
         * @param patternIndex index to pattern
         * @param taxonIndex   index to taxon
         */
        public abstract void computeNucReadCountsLikelihoods(
                final int[] reads,
                final int patternIndex,
                final int taxonIndex
        );

        /**
         * Get the parameters situation 1 where the wild type (0/0) applies.
         *
         * @param out communicated arguments
         */
        public abstract void getWildTypeNucReadCountModelParams(double[] out);


        //**********************************************
        //*             Overridden methods             *
        //**********************************************

        @Override
        public void initAndValidate() {
            if (effSeqErrRateInput.get().getValue() < 0.0 ||
                    (effSeqErrRateInput.get().isEstimated() &&
                            (effSeqErrRateInput.get().getLower() < 0.0 ||
                                    effSeqErrRateInput.get().getUpper() < 0.0 ||
                                    effSeqErrRateInput.get().getUpper() > 1.0)))
                throw new IllegalArgumentException("effSeqErrRate and its bounds should be defined as no smaller than 0.");
            else
                effSeqErrRate = effSeqErrRateInput.get();

            if (shapeCtrl1Input.get().getValue() < 0.0 ||
                    (shapeCtrl1Input.get().isEstimated() &&
                            (shapeCtrl1Input.get().getLower() < 0.0 ||
                                    shapeCtrl1Input.get().getUpper() < 0.0)))
                throw new IllegalArgumentException("shapeCtrl1 and its bounds should be defined as no smaller than 0 (" +
                        this.getClass().getName() + ")");
            else
                shapeCtrl1 = shapeCtrl1Input.get();

            if (shapeCtrl2Input.get().getValue() < 0.0 ||
                    (shapeCtrl2Input.get().isEstimated() &&
                            (shapeCtrl2Input.get().getLower() < 0.0 ||
                                    shapeCtrl2Input.get().getUpper() < 0.0)))
                throw new IllegalArgumentException("shapeCtrl2 and its bounds should be defined as no smaller than 0 (" +
                        this.getClass().getName() + ")");
            else
                shapeCtrl2 = shapeCtrl2Input.get();
        } // initAndValidate

        @Override
        protected void store() {
            super.store();

            System.arraycopy(
                    currentNucReadCountsLikelihoodsNodeIndex,
                    0,
                    storedNucReadCountsLikelihoodsNodeIndex,
                    0,
                    nrOfTaxa
            );
        } // store

        @Override
        protected void restore() {
            super.restore();

            int[] tmp = currentNucReadCountsLikelihoodsNodeIndex;
            currentNucReadCountsLikelihoodsNodeIndex = storedNucReadCountsLikelihoodsNodeIndex;
            storedNucReadCountsLikelihoodsNodeIndex = tmp;
        } // restore


        //*********************************************
        //*                  Methods                  *
        //*********************************************

        /**
         * Initialize parameters with passed arguments from ScsTreeLikelihood.
         * Must be called by every inheriting class at the beginning.
         *
         * @param alignment      apparently
         * @param useLogPartials apparently
         */
        public void deeplyInitialize(
                final ScsAlignment alignment,
                final int stateNum,
                final boolean useLogPartials
        ) {
            nrOfPatterns = alignment.getPatternCount();
            nrOfTaxa = alignment.getTaxonCount();
            nrOfStates = stateNum;
            this.useLogPartials = useLogPartials;
            nucReadCountsLikelihoods = new double[2][nrOfTaxa][nrOfPatterns * nrOfCategories];

            currentNucReadCountsLikelihoodsNodeIndex = new int[nrOfTaxa];
            storedNucReadCountsLikelihoodsNodeIndex = new int[nrOfTaxa];
        } // deeplyInitialize

        /**
         * Duplicate some important variables entirely or within ranged patterns.
         * Used when there are more than one thread.
         * Must be called by every inheriting class at the beginning.
         *
         * @param target         target instance storing variables
         * @param targetData     data corresponding to the target raw read counts model
         */
        public void duplicate(
                NucReadCountsModelInterface.Base target,
                final ScsAlignment targetData
        ) {
            target.nrOfPatterns = targetData.getPatternCount();
            target.nrOfTaxa = this.nrOfTaxa;
            target.nrOfStates = this.nrOfStates;
            target.useLogPartials = this.useLogPartials;
            target.nucReadCountsLikelihoods = new double[2][nrOfTaxa][target.nrOfPatterns * nrOfCategories];

            target.currentNucReadCountsLikelihoodsNodeIndex = new int[nrOfTaxa];
            target.storedNucReadCountsLikelihoodsNodeIndex = new int[nrOfTaxa];
        } // duplicate

        /**
         * During likelihood initialization, likelihoods should be copied for the sake of accurate tree likelihood computation.
         *
         * @param taxonIndex apparently
         */
        public void storeLikelihoods(int taxonIndex) {
            System.arraycopy(
                    nucReadCountsLikelihoods[currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex],
                    0,
                    nucReadCountsLikelihoods[1 - currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex],
                    0,
                    nrOfPatterns * nrOfCategories
            );
        } // storeLikelihoods

        public void setNucReadCountsLikelihoodsNodeForUpdate(int nodeIndex) {
            currentNucReadCountsLikelihoodsNodeIndex[nodeIndex] = 1 - currentNucReadCountsLikelihoodsNodeIndex[nodeIndex];
        } // setSeqCovLikelihoodsNodeForUpdate

        public int getNrOfCategories() {
            return nrOfCategories;
        } // getNrOfCategories

        /**
         * Get the nucleotide read counts likelihood at {@param index}.
         *
         * @param taxonIndex apparently
         * @param index      index to the last dimension
         * @return sequencing coverage likelihood
         */
        public double getNucReadCountsLikelihood(
                final int taxonIndex,
                final int index
        ) {
            return nucReadCountsLikelihoods[currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index];
        } // getNucReadCountsLikelihood

        public int getNrOfParamsInReadCountsModel() {
            return nrOfNucReadCountsModelParam;
        } // getNrOfParamsInReadCountsModel

        protected boolean sanityCheckLengthWildTypeSituationParams(double[] out) {
            return out != null && out.length == this.nrOfNucReadCountsModelParam;
        } // sanityCheckLengthWildTypeSituationParams

        public ScsParametricDistribution getDistribution() {
            return variantDensity;
        } // getNucReadCountsModelDistribution

    }

}
