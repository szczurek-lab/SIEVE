package beast.evolution.rawreadcountsmodel.nucreadcountsmodel;

import beast.math.distributions.DirichletMultinomial;

public class NucReadCountsModelFiniteMuDM extends NucReadCountsModelInterface.Base {


    //**********************************************
    //*             Overridden methods             *
    //**********************************************

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        this.variantDensity = new DirichletMultinomial();

        // alpha_1, alpha_2, alpha_3, alpha_4, alpha_0
        this.nrOfNucReadCountsModelParam = 5;

        // 4 categories (situations) of nucleotide read counts likelihood
        this.nrOfCategories = 4;
    }

    /**
     * Compute nucleotide read counts likelihoods.
     * @param reads        reads for current taxa under current pattern
     * @param patternIndex index to pattern
     * @param taxonIndex   index to taxon
     *
     */
    @Override
    public void computeNucReadCountsLikelihoods(
            final int[] reads,
            int patternIndex,
            int taxonIndex
    ) {
        final int index = patternIndex * nrOfCategories;

        final double w1 = shapeCtrl1.getValue();
        final double w2 = shapeCtrl2.getValue();

        final double oneThirdSeqErr = effSeqErrRate.getValue() / 3;
        final double oneMinusSeqErr = 1 - effSeqErrRate.getValue();
        final double halfMinusOneThirdSeqErr = 0.5 - effSeqErrRate.getValue() / 3;

        final double oneThirdSeqErrW1 = oneThirdSeqErr * w1;
        final double oneMinusSeqErrW1 = oneMinusSeqErr * w1;

        final double oneThirdSeqErrW2 = oneThirdSeqErr * w2;
        final double halfMinusOneThirdSeqErrW2 = halfMinusOneThirdSeqErr * w2;

        // 0/0
        this.nucReadCountsLikelihoods[currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index] = useLogPartials ?
                ((DirichletMultinomial) variantDensity).logDensity(
                        reads,
                        new double[]{oneThirdSeqErrW1, oneThirdSeqErrW1, oneThirdSeqErrW1, oneMinusSeqErrW1, w1}
                ) :
                ((DirichletMultinomial) variantDensity).density(
                        reads,
                        new double[]{oneThirdSeqErrW1, oneThirdSeqErrW1, oneThirdSeqErrW1, oneMinusSeqErrW1, w1}
                );

        // 1/1
        this.nucReadCountsLikelihoods[currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index + 1] = useLogPartials ?
                ((DirichletMultinomial) variantDensity).logDensity(
                        reads,
                        new double[]{oneMinusSeqErrW1, oneThirdSeqErrW1, oneThirdSeqErrW1, oneThirdSeqErrW1, w1}
                ) :
                ((DirichletMultinomial) variantDensity).density(
                        reads,
                        new double[]{oneMinusSeqErrW1, oneThirdSeqErrW1, oneThirdSeqErrW1, oneThirdSeqErrW1, w1}
                );

        // 1/1'
        this.nucReadCountsLikelihoods[currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index + 2] = useLogPartials ?
                ((DirichletMultinomial) variantDensity).logDensity(
                        reads,
                        new double[]{halfMinusOneThirdSeqErrW2, halfMinusOneThirdSeqErrW2, oneThirdSeqErrW2, oneThirdSeqErrW2, w2}
                ) :
                ((DirichletMultinomial) variantDensity).density(
                        reads,
                        new double[]{halfMinusOneThirdSeqErrW2, halfMinusOneThirdSeqErrW2, oneThirdSeqErrW2, oneThirdSeqErrW2, w2}
                );

        // 0/1
        this.nucReadCountsLikelihoods[currentNucReadCountsLikelihoodsNodeIndex[taxonIndex]][taxonIndex][index + 3] = useLogPartials ?
                ((DirichletMultinomial) variantDensity).logDensity(
                        reads,
                        new double[]{halfMinusOneThirdSeqErrW2, oneThirdSeqErrW2, oneThirdSeqErrW2, halfMinusOneThirdSeqErrW2, w2}
                ) :
                ((DirichletMultinomial) variantDensity).density(
                        reads,
                        new double[]{halfMinusOneThirdSeqErrW2, oneThirdSeqErrW2, oneThirdSeqErrW2, halfMinusOneThirdSeqErrW2, w2}
                );
    } // computeNucReadCountsLikelihoods

    @Override
    public void getWildTypeNucReadCountModelParams(double[] out) {
        if (!sanityCheckLengthWildTypeSituationParams(out))
            throw new IllegalArgumentException("Error! The length of parameter accepting wild type arguments " +
                    "should be " + this.nrOfNucReadCountsModelParam + ", but " + out.length + " provided.");

        final double w1 = shapeCtrl1.getValue();

        final double oneThirdSeqErr = effSeqErrRate.getValue() / 3;
        final double oneMinusSeqErr = 1 - effSeqErrRate.getValue();

        out[0] = out[1] = out[2] = oneThirdSeqErr * w1;
        out[3] = oneMinusSeqErr * w1;
        out[4] = w1;
    } // getWildTypeNucReadCountModelParams

}
