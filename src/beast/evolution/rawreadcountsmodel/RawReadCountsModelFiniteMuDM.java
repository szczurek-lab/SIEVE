package beast.evolution.rawreadcountsmodel;

import beast.core.Description;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static beast.math.util.MathFunctions.logSumExp;
import static beast.math.util.MathFunctions.maxIndices;

@Description("Model of nucleotide read counts described by Dirichlet-multinomial distribution.")
public class RawReadCountsModelFiniteMuDM extends RawReadCountsModelInterface.Base {


    //**********************************************
    //*             Overridden methods             *
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
    public double computeMixedLikelihood(
            final int matrixIndex,
            final int patternIndex,
            final int taxonIndex,
            final int genotypeIndex
    ) {
        return computeMixedLikelihoodCore(
                genotypeIndex,
                taxonIndex,
                matrixIndex * nrOfPatterns * modeledAllelesSize + patternIndex * modeledAllelesSize,
                patternIndex * nrOfCategories,
                null
        );
    } // computeMixedLikelihood

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
    public double computeMixedLikelihood(
            final int matrixIndex,
            final int patternIndex,
            final int taxonIndex,
            final int genotypeIndex,
            List<Integer>[] MLCompIndex
    ) {
        double[] comp = new double[modeledAllelesSize];
        final double lh = computeMixedLikelihoodCore(
                genotypeIndex,
                taxonIndex,
                matrixIndex * nrOfPatterns * modeledAllelesSize + patternIndex * modeledAllelesSize,
                patternIndex * nrOfCategories,
                comp
        );

        if (MLCompIndex[0] == null)
            MLCompIndex[0] = new ArrayList<>();
        else
            MLCompIndex[0].clear();
        MLCompIndex[0].addAll(maxIndices(Arrays.stream(comp).boxed().toArray(Double[]::new)));

        return lh;
    } // computeMixedLikelihood

    /**
     * Get the log likelihoods for ADOs.
     *
     * @param matrixIndex   apparently
     * @param patternIndex  apparently
     * @param taxonIndex    apparently
     * @param genotypeIndex apparently
     * @return log likelihoods of ADO for {@param genotypeIndex}
     */
    public double[] getAdoLogLikelihoods(
            final int matrixIndex,
            final int patternIndex,
            final int taxonIndex,
            final int genotypeIndex
    ) {
        double[] comp = new double[modeledAllelesSize];
        computeMixedLikelihoodCore(
                genotypeIndex,
                taxonIndex,
                matrixIndex >= nrOfMatrices ? patternIndex * modeledAllelesSize : matrixIndex * nrOfPatterns * modeledAllelesSize + patternIndex * modeledAllelesSize,
                patternIndex * nrOfCategories,
                comp
        );

        return comp;
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * Core function to compute the mixed likelihood.
     *
     * @param genotypeIndex                  apparently
     * @param taxonIndex                     apparently
     * @param indexToSeqCovLikelihood        apparently
     * @param indexToNucReadCountsLikelihood apparently
     * @param comp                           a double array to store component likelihoods
     * @return the mixed likelihood
     */
    private double computeMixedLikelihoodCore(
            final int genotypeIndex,
            final int taxonIndex,
            final int indexToSeqCovLikelihood,
            final int indexToNucReadCountsLikelihood,
            double[] comp
    ) {
        final double theta = adoRate.getValue();

        double[] tmp = comp == null ? new double[modeledAllelesSize] : comp;

        switch (genotypeIndex) {
            case 0:
                // 0/0
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = theta * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(2) + Math.log(theta) + Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = 2 * theta * (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 1:
                // 0/1
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) - Math.log(2) + logSumExp(new double[]{nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood), nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)}) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = (theta / 2) * (nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(theta) + Math.log(1 - theta) + logSumExp(new double[]{nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood), nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)}) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = theta * (1 - theta) * (nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 2:
                // 1/1
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = theta * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(2) + Math.log(theta) + Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = 2 * theta * (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 3:
                // 1/1'
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = theta * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(2) + Math.log(theta) + Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = 2 * theta * (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            default:
                throw new IllegalStateException("Unexpected genotype index: " + genotypeIndex + " (" + this.getClass().getName() + ")");
        }
    } // computeMixedLikelihood

}
