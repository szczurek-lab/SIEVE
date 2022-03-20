package beast.evolution.rawreadcountsmodel;

import beast.core.Description;

import java.util.List;

@Description("An error model compatible with ScsFiniteMuModel")
public class RawReadCountsModelFiniteMuBB extends RawReadCountsModelInterface.Base {


    @Override
    public double computeMixedLikelihood(int matrixIndex, int patternIndex, int taxonIndex, int genotypeIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public double computeMixedLikelihood(int matrixIndex, int patternIndex, int taxonIndex, int genotypeIndex, List<Integer>[] MLCompIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    }

    @Override
    public double[] getAdoLogLikelihoods(int matrixIndex, int patternIndex, int taxonIndex, int genotypeIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

//    @Override
//    public void initAndValidate() {
//        super.initAndValidate();
//
//        this.variantDensity = new BetaBinomial();
//
//        this.nrOfParamsInReadCountsModel = 2;
//
//        this.stateCount = 3;
//
//        // 3 situations
//        this.situationLikelihoods = new double[3];
//    } // initAndValidate

//    /**
//     * update situation likelihoods for three situations
//     *
//     * @param reads reads for current taxa under current pattern
//     */
//    @Override
//    public void updateSituationLikelihoods(final int[] reads) {
//        final double f1 = seqErr.getValue() / 3;
//        final double f2 = seqErr.getValue();
//        final double w1 = shapeCtrl1.getValue();
//        final double f3 = 0.5 - f1;
//        final double w2 = shapeCtrl2.getValue();
//
//        final long seqCov = reads[0];
//        final long variantReads = reads[1];
//
//        situationLikelihoods[0] = ((BetaBinomial) variantDensity).density(variantReads, seqCov, f1 * w1, w1 * (1 - f1));
//        situationLikelihoods[1] = ((BetaBinomial) variantDensity).density(seqCov - variantReads, seqCov, f2 * w1, w1 * (1 - f2));
//        situationLikelihoods[2] = ((BetaBinomial) variantDensity).density(variantReads, seqCov, f3 * w2, w2 * (1 - f3));
//    } // updateSituationLikelihoods

//    /**
//     * get the overall likelihood of raw read counts at the leaves
//     * account for sequencing coverage (by negative binomial distribution) and
//     * variant reads (by beta-binomial distribution or Dirichlet-multinomial distribution)
//     *
//     * @param genotypeIndex  apparently
//     * @param matrixIndex    apparently
//     * @param patternIndex   apparently
//     * @param taxonIndex     apparently
//     * @param useLogPartials apparently
//     * @param MLCompIndex    index to the maximum likelihood component, which corresponding to the number of sequenced alleles
//     * @return likelihoods
//     */
//    @Override
//    public double getOverallLikelihood(
//            final int genotypeIndex,
//            final int matrixIndex,
//            final int patternIndex,
//            final int taxonIndex,
//            final boolean useLogPartials,
//            List<Integer>[] MLCompIndex
//    ) {
//        double lh = 0.0;
//        final double theta = adoRate.getValue();
//        final int index = matrixIndex * patternCount * taxaCount * modeledAllelesSize + patternIndex * taxaCount * modeledAllelesSize + taxonIndex * modeledAllelesSize;
//
//        if (MLCompIndex[0] == null)
//            MLCompIndex[0] = new ArrayList<>();
//        else
//            MLCompIndex[0].clear();
//
//        switch (genotypeIndex) {
//            case 0:
//                // 0/0
//                if (singleADO) {
//                    // single ado
//
//                    if (useLogPartials) {
//                        dataLikelihoods[index] = (1 - theta) * situationLikelihoods[0] * seqCovLikelihoods[index + 1];
//                        dataLikelihoods[index + 1] = theta * situationLikelihoods[0] * seqCovLikelihoods[index];
//
//                        lh = logSumExp(new double[]{dataLikelihoods[index], dataLikelihoods[index + 1]});
//                    } else {
//                        dataLikelihoods[index] = (1 - theta) * situationLikelihoods[0] * seqCovLikelihoods[index + 1];
//                        dataLikelihoods[index + 1] = theta * situationLikelihoods[0] * seqCovLikelihoods[index];
//
//                        lh = dataLikelihoods[index] + dataLikelihoods[index + 1];
//                    }
//                } else {
//                    // locus ado
//                    lh = Math.pow(1 - theta, 2) * situationLikelihoods[0] * seqCovLikelihoods[taxonIndex][index + 2] + 2 * theta * (1 - theta) * situationLikelihoods[0] * seqCovLikelihoods[taxonIndex][index + 1] + Math.pow(theta, 2) * seqCovLikelihoods[taxonIndex][index];
//                }
//                break;
//            case 1:
//                // 0/1
//                if (singleADO) {
//                    // single ado
//                    lh = (1 - theta) * situationLikelihoods[2] * seqCovLikelihoods[taxonIndex][index + 2] + (theta / 2) * (situationLikelihoods[0] + situationLikelihoods[1]) * seqCovLikelihoods[taxonIndex][index + 1];
//                } else {
//                    // locus ado
//                    lh = Math.pow(1 - theta, 2) * situationLikelihoods[2] * seqCovLikelihoods[taxonIndex][index + 2] + theta * (1 - theta) * (situationLikelihoods[0] + situationLikelihoods[1]) * seqCovLikelihoods[taxonIndex][index + 1] + Math.pow(theta, 2) * seqCovLikelihoods[taxonIndex][index];
//                }
//                break;
//            case 2:
//                // 1/1
//                if (singleADO) {
//                    // single ado
//                    lh = ((1 - theta) / 2) * (situationLikelihoods[1] + situationLikelihoods[2]) * seqCovLikelihoods[taxonIndex][index + 2] + theta * situationLikelihoods[1] * seqCovLikelihoods[taxonIndex][index + 1];
//                } else {
//                    // locus ado
//                    lh = (Math.pow(1 - theta, 2) / 2) * (situationLikelihoods[1] + situationLikelihoods[2]) * seqCovLikelihoods[taxonIndex][index + 2] + 2 * theta * (1 - theta) * situationLikelihoods[1] * seqCovLikelihoods[taxonIndex][index + 1] + Math.pow(theta, 2) * seqCovLikelihoods[taxonIndex][index];
//                }
//                break;
//            default:
//                throw new IllegalStateException("Unexpected genotype index: " + genotypeIndex + " (" + this.getClass().getName() + ")");
//        }
//
//        return lh;
//    } // getOverallLikelihood

//    /**
//     * get the parameters situation 1 where the wild type (0/0) applies
//     *
//     * @param out communicated arguments
//     */
//    @Override
//    public void getWildTypeSituationParams(double[] out) {
//        if (!sanityCheckLengthWildTypeSituationParams(out))
//            throw new IllegalArgumentException("Error! The length of parameter accepting wild type arguments " +
//                    "should be " + this.nrOfParamsInReadCountsModel + ", but " + out.length + " provided.");
//
//        final double f1 = seqErr.getValue() / 3;
//        final double w1 = shapeCtrl1.getValue();
//
//        out[0] = f1 * w1;
//        out[1] = w1 * (1 - f1);
//    } // getWildTypeSituationParams

}
