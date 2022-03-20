package beast.evolution.rawreadcountsmodel.seqcovmodel;

import beast.evolution.alignment.FilteredScsAlignment;
import beast.evolution.alignment.ScsAlignment;

import java.io.PrintStream;
import java.util.*;

public class PrivateAllelicSeqCovModel extends SeqCovModelInterface.Base {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    protected double[] allelicSeqCovPerPattern; // #matrices * #patterns
    protected double[] allelicSeqCovRawVarPerPattern; // #matrices * #patterns

    /**
     * For the purpose of debug.
     * Most recently accepted move's allelic sequencing coverage and its raw variance.
     */
    protected double[] storedAllelicSeqCovPerPattern; // #matrices * #patterns
    protected double[] storedAllelicSeqCovRawVarPerPattern; // #matrices * #patterns


    //********************************************
    //*             Abstract methods             *
    //********************************************


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

        if (initSharedAllelicSeqCovAndRawVarInput.get() != null)
            throw new IllegalArgumentException("\"initSharedAllelicSeqCovAndRawVar\" is only valid for ExploredSharedAllelicSeqCovModel.");
    } // initAndValidate

    /**
     * Initialize parameters with passed arguments from ScsTreeLikelihood.
     * Must be called by every inheriting class at the beginning.
     *  @param alignment     apparently
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

        final int nrOfLoci = alignment.getLociNr();

        sizeFactors = new double[nrOfTaxa];
        processedSeqCov = new int[nrOfPatterns * nrOfTaxa];

        computeSizeFactors(alignment);

        if (modeledAllelesSize == 1) {
            // the combination of the number of sequenced alleles for all taxa is deterministic

            needToUpdate = false;
            nrOfMatrices = 1;
        } else if (modeledAllelesSize >= 1) {
            // combinations of the number of sequenced alleles for all taxa are dynamic

            needToUpdate = true;

            MLNrOfSequencedAllelesPatterns = new ArrayList[nrOfMatrices * nrOfPatterns];
            MLNrOfSequencedAllelesWeights = new ArrayList[nrOfMatrices * nrOfPatterns];
        } else
            throw new IllegalArgumentException("At least one allele should be modeled in the evolutionary model.");

        // 0: for likelihood computation
        // 1: for debugging
        seqCovLikelihoods = new double[2][nrOfTaxa][nrOfMatrices * nrOfPatterns * modeledAllelesSize];

        allelicSeqCovPerPattern = new double[nrOfMatrices * nrOfPatterns];
        allelicSeqCovRawVarPerPattern = new double[nrOfMatrices * nrOfPatterns];
        storedAllelicSeqCovPerPattern = new double[nrOfMatrices * nrOfPatterns];
        storedAllelicSeqCovRawVarPerPattern = new double[nrOfMatrices * nrOfPatterns];

        // if in variant calling mode
        if (this.inVariantCallingMode) {
            String[] allelicSeqCovArray = allelicSeqCovArrayInput.get().trim().split(";");
            String[] allelicSeqCovRawVarArray = allelicSeqCovRawVarArrayInput.get().trim().split(";");

            if (allelicSeqCovArray.length != nrOfMatrices || allelicSeqCovRawVarArray.length != nrOfMatrices)
                throw new IllegalArgumentException("Unmatched configuration of the number of categories. Please " +
                        "do not edit the configuration document.");

            for (int i = 0; i < nrOfMatrices; i++) {
                String[] parsedCov = allelicSeqCovArray[i].trim().split(",");
                if (parsedCov.length != nrOfLoci)
                    throw new IllegalArgumentException("Unmatched number of loci. Expected: " + nrOfLoci +
                            ". Provided in allelicSeqCovArrayInput: " + parsedCov.length + ".");

                String[] parsedVar = allelicSeqCovRawVarArray[i].trim().split(",");
                if (parsedVar.length != nrOfLoci)
                    throw new IllegalArgumentException("Unmatched number of loci. Expected: " + nrOfLoci +
                            ". Provided in allelicSeqCovRawVarArrayInput: " + parsedVar.length + ".");

                double[] cov = Arrays.stream(parsedCov).mapToDouble(Double::parseDouble).toArray();
                double[] var = Arrays.stream(parsedVar).mapToDouble(Double::parseDouble).toArray();

                Set<Integer> addedPatterns = new HashSet<>();
                for (int j = 0; j < nrOfLoci; j++) {
                    int pattern = alignment.getPatternIndex(i);
                    int index = i * nrOfPatterns + pattern;

                    if (addedPatterns.contains(pattern)) {
                        if (allelicSeqCovPerPattern[index] != cov[j] ||
                                allelicSeqCovRawVarPerPattern[index] != var[j])
                            throw new IllegalArgumentException("Ambiguous values mapped to the same pattern for locus " + j + ".");
                    } else {
                        addedPatterns.add(pattern);

                        allelicSeqCovPerPattern[index] = cov[j];
                        allelicSeqCovRawVarPerPattern[index] = var[j];
                    }
                }
            }
        } else if (needToUpdate) {
            Arrays.fill(allelicSeqCovPerPattern, allelicSeqCov.getValue());
            Arrays.fill(allelicSeqCovRawVarPerPattern, allelicSeqCovRawVar.getValue());
        } else
            computeDeterministicAllelicSeqCovAndVar();

        storeAllelicInfo();
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
            Base target,
            final ScsAlignment targetData
    ) {
        super.duplicate(target, targetData);

        final ScsAlignment originalData = ((FilteredScsAlignment) targetData).alignmentInput.get();

        final int[] filteredLoci = ((FilteredScsAlignment) targetData).indices();

        // allocate memory and copy values
        target.sizeFactors = new double[nrOfTaxa];
        System.arraycopy(this.sizeFactors, 0, target.sizeFactors, 0, nrOfTaxa);

        target.processedSeqCov = new int[target.nrOfPatterns * nrOfTaxa];

        if (target.needToUpdate) {
            // combinations of the number of sequenced alleles for all taxa are dynamic

            target.MLNrOfSequencedAllelesPatterns = new ArrayList[nrOfMatrices * target.nrOfPatterns];
            target.MLNrOfSequencedAllelesWeights = new ArrayList[nrOfMatrices * target.nrOfPatterns];
        }

        // 0: for likelihood computation
        // 1: for debugging
        target.seqCovLikelihoods = new double[2][nrOfTaxa][nrOfMatrices * target.nrOfPatterns * modeledAllelesSize];

        ((PrivateAllelicSeqCovModel) target).allelicSeqCovPerPattern = new double[nrOfMatrices * target.nrOfPatterns];
        ((PrivateAllelicSeqCovModel) target).allelicSeqCovRawVarPerPattern = new double[nrOfMatrices * target.nrOfPatterns];
        ((PrivateAllelicSeqCovModel) target).storedAllelicSeqCovPerPattern = new double[nrOfMatrices * target.nrOfPatterns];
        ((PrivateAllelicSeqCovModel) target).storedAllelicSeqCovRawVarPerPattern = new double[nrOfMatrices * target.nrOfPatterns];

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {
            for (int locusIndex = 0; locusIndex < filteredLoci.length; locusIndex++) {

                // get pattern index
                int patternIndexOri = matrixIndex * this.nrOfPatterns + originalData.getPatternIndex(filteredLoci[locusIndex]);
                int patternIndex = matrixIndex * target.nrOfPatterns + targetData.getPatternIndex(locusIndex);

                ((PrivateAllelicSeqCovModel) target).allelicSeqCovPerPattern[patternIndex] = this.allelicSeqCovPerPattern[patternIndexOri];
                ((PrivateAllelicSeqCovModel) target).allelicSeqCovRawVarPerPattern[patternIndex] = this.allelicSeqCovRawVarPerPattern[patternIndexOri];

                ((PrivateAllelicSeqCovModel) target).storedAllelicSeqCovPerPattern[patternIndex] = this.storedAllelicSeqCovPerPattern[patternIndexOri];
                ((PrivateAllelicSeqCovModel) target).storedAllelicSeqCovRawVarPerPattern[patternIndex] = this.storedAllelicSeqCovRawVarPerPattern[patternIndexOri];

                if (matrixIndex == 0)
                    System.arraycopy(this.processedSeqCov, patternIndexOri * nrOfTaxa, target.processedSeqCov, patternIndex * nrOfTaxa, nrOfTaxa);

            }
        }
    } // duplicate

    @Override
    public void setSeqCovLikelihoodsNodeForUpdate(int nodeIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    } // setSeqCovLikelihoodsNodeForUpdate

    @Override
    public void storeLikelihoods(int taxonIndex) {
        // Nothing needs to be done.
    } // storeLikelihoods

    /**
     * Called during initialization of tree likelihoods for internal variables expanding / copying to other matrices.
     */
    @Override
    public void expandVariables(int taxonIndex) {
        for (int matrixIndex = 1; matrixIndex < nrOfMatrices; matrixIndex++) {
            System.arraycopy(
                    seqCovLikelihoods[0][taxonIndex],
                    0,
                    seqCovLikelihoods[1][taxonIndex],
                    matrixIndex * nrOfPatterns * modeledAllelesSize,
                    nrOfPatterns * modeledAllelesSize
            );
        }
    } // expandVariables

    /**
     * For the purpose of debugging.
     * Only the most recent accepted allelic sequencing coverage and raw variance should be stored.
     */
    @Override
    public void storeAllelicInfo() {
        System.arraycopy(allelicSeqCovPerPattern, 0, storedAllelicSeqCovPerPattern, 0, allelicSeqCovPerPattern.length);
        System.arraycopy(allelicSeqCovRawVarPerPattern, 0, storedAllelicSeqCovRawVarPerPattern, 0, allelicSeqCovRawVarPerPattern.length);
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
        final int index = matrixIndex * nrOfPatterns * modeledAllelesSize + patternIndex * modeledAllelesSize;

        for (int i = 0; i < modeledAllelesSize; i++) {
            if (inDebugMode)
                seqCovLikelihoods[1][taxonIndex][index + i] = computeSeqCovLikelihoodPerAllele(
                        seqCov,
                        modeledAlleles[i],
                        taxonIndex,
                        matrixIndex * nrOfPatterns + patternIndex
                );
            else
                seqCovLikelihoods[0][taxonIndex][index + i] = computeSeqCovLikelihoodPerAllele(
                        seqCov,
                        modeledAlleles[i],
                        taxonIndex,
                        matrixIndex * nrOfPatterns + patternIndex
                );
        }
    } // computeSeqCovLikelihood

    @Override
    public double getSeqCovLikelihood(
            final int taxonIndex,
            final int index
    ) {
        return inDebugMode ? seqCovLikelihoods[1][taxonIndex][index] : seqCovLikelihoods[0][taxonIndex][index];
    } // getSeqCovLikelihood

    /**
     * Compute allelic sequencing coverage and its raw variance.
     * Should be called at the post processing stage of computing tree likelihood.
     *
     * @param changedPatterns matrices and patterns which should be updated: (matrixIndex, patternIndex)
     */
    public void updateAllelicSeqCovAndRawVar(List<List<Integer>> changedPatterns) {
        Iterator<List<Integer>> itr = changedPatterns.iterator();

        while (itr.hasNext()) {
            // in a form of (matrixIndex, patternIndex)
            List<Integer> pair = itr.next();
            final int currentIndex = pair.get(0) * nrOfPatterns + pair.get(1);

            double sumSeqCov = 0.0;
            double sumSeqCovVar = 0.0;

            /*
             * 0: for allelic sequencing coverage
             * 1: for allelic sequencing coverage raw variance
             */
            int[] sumWeights = new int[2];

            for (int combIndex = 0; combIndex < MLNrOfSequencedAllelesPatterns[currentIndex].size(); combIndex++) {

                // the number of cells that have the same number of sequenced alleles for a pattern
                int[] h = new int[modeledAllelesSize];

                // the indices of taxa corresponding to each entry in h
                int[][] indices = new int[modeledAllelesSize][nrOfTaxa];

                // the number of unique occurrences of the number of sequenced alleles
                // first: possessed by at least 1 taxon
                // second: possessed by at least 2 taxa
                int[] N = new int[2];

                // intermediate values
                double[] q = new double[modeledAllelesSize];
                double[] w = new double[modeledAllelesSize];
                double[] z = new double[modeledAllelesSize];

                // set h and indices
                for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {

                    final int index = Arrays.binarySearch(modeledAlleles, MLNrOfSequencedAllelesPatterns[currentIndex].get(combIndex)[taxonIndex]);

                    indices[index][h[index]++] = taxonIndex;
                }

                // compute q, w, and z
                for (int allelesIndex = 0; allelesIndex < modeledAllelesSize; allelesIndex++) {

                    // compute q
                    if (h[allelesIndex] > 0) {
                        N[0]++;

                        double sumQ = 0.0;
                        int index;
                        for (int i = 0; i < h[allelesIndex]; i++) {
                            index = indices[allelesIndex][i];
                            sumQ += processedSeqCov[pair.get(1) * nrOfTaxa + index] / sizeFactors[index];
                        }

                        q[allelesIndex] = sumQ / h[allelesIndex];
                    }

                    // compute w and z
                    if (h[allelesIndex] > 1) {
                        N[1]++;

                        double sumW = 0.0;
                        double sumZ = 0.0;
                        int index;
                        for (int i = 0; i < h[allelesIndex]; i++) {
                            index = indices[allelesIndex][i];
                            sumW += Math.pow(processedSeqCov[pair.get(1) * nrOfTaxa + index] / sizeFactors[index] - q[allelesIndex], 2);
                            sumZ += (1 / sizeFactors[index]);
                        }

                        w[allelesIndex] = sumW / (h[allelesIndex] - 1);
                        z[allelesIndex] = w[allelesIndex] - q[allelesIndex] * sumZ / h[allelesIndex];

                        // sometimes z is negative, rule them out
                        if (z[allelesIndex] <= 0) {
                                /*
                                System.out.println("\nSize factors: " + Arrays.toString(sizeFactors));
                                System.out.println("Number of alleles: " + Arrays.toString(MLNrOfSequencedAllelesPatterns[currentIndex].get(combIndex)));
                                System.out.println("Number of alleles goes wrong: " + allelesIndex);
                                System.out.println("h: " + h[allelesIndex]);
                                System.out.println("q: " + q[allelesIndex]);
                                System.out.println("w: " + w[allelesIndex]);
                                System.out.println("z: " + z[allelesIndex] + "\n");
                                 */

                            z[allelesIndex] = 0.0;
                            N[1]--;
                        }
                    }
                }

                // compute allelic sequencing coverage and its raw variance
                double sumT = 0.0;
                double sumV = 0.0;
                for (int allelesIndex = 0; allelesIndex < modeledAllelesSize; allelesIndex++) {
                    sumT += q[allelesIndex] / (modeledAlleles[allelesIndex] + numericalStabilizer);
                    sumV += z[allelesIndex] / Math.pow(modeledAlleles[allelesIndex] + numericalStabilizer, 2);
                }

                sumSeqCov += MLNrOfSequencedAllelesWeights[currentIndex].get(combIndex) * sumT / N[0];
                sumWeights[0] += MLNrOfSequencedAllelesWeights[currentIndex].get(combIndex);
                if (N[1] > 0) {
                    sumSeqCovVar += MLNrOfSequencedAllelesWeights[currentIndex].get(combIndex) * sumV / N[1];
                    sumWeights[1] += MLNrOfSequencedAllelesWeights[currentIndex].get(combIndex);
                }
            }

            allelicSeqCovPerPattern[currentIndex] = sumSeqCov / sumWeights[0];
            if (sumSeqCovVar > 0)
                allelicSeqCovRawVarPerPattern[currentIndex] = sumSeqCovVar / sumWeights[1];

            // record matrices and patterns whose allelic sequencing coverage and raw variance are not changed
            if (allelicSeqCovPerPattern[currentIndex] == storedAllelicSeqCovPerPattern[currentIndex] &&
                    allelicSeqCovRawVarPerPattern[currentIndex] == storedAllelicSeqCovRawVarPerPattern[currentIndex])
                itr.remove();

        }
    } // updateAllelicSeqCovAndRawVar

    /**
     * log sampled allelic sequencing coverage and raw variance
     *
     * @param index which pattern is about to be logged?
     * @param out   apparently
     */
    public void logCovar(int index, PrintStream out) {
        for (int i = 0; i != nrOfMatrices; i++) {
            out.print(allelicSeqCovPerPattern[index + i * nrOfPatterns] + "," + allelicSeqCovRawVarPerPattern[index + i * nrOfPatterns]);

            if (i != nrOfMatrices - 1)
                out.print(";");
        }
    } // logCovar


    //*********************************************
    //*                  Methods                  *
    //*********************************************

    /**
     * Compute allelic sequencing coverage and its raw variance if the combination of the number of sequenced alleles
     * for all taxa is deterministic.
     */
    private void computeDeterministicAllelicSeqCovAndVar() {
        for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

            // intermediate values
            double q, w, z;
            q = w = z = 0.0;

            // compute q
            for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++)
                q += processedSeqCov[patternIndex * nrOfTaxa + taxonIndex] / sizeFactors[taxonIndex];

            q /= nrOfTaxa;

            // compute w and z
            if (nrOfTaxa > 1) {
                for (int taxonIndex = 0; taxonIndex < nrOfTaxa; taxonIndex++) {
                    w += Math.pow(processedSeqCov[patternIndex * nrOfTaxa + taxonIndex] / sizeFactors[taxonIndex] - q, 2);
                    z += (1 / sizeFactors[taxonIndex]);
                }

                w /= (nrOfTaxa - 1);
                z = w - q * z / nrOfTaxa;
            }

            // compute and save t and v
            allelicSeqCovPerPattern[patternIndex] = q / (modeledAlleles[0] + numericalStabilizer);

            if (z > 0)
                allelicSeqCovRawVarPerPattern[patternIndex] = z / Math.pow(modeledAlleles[0] + numericalStabilizer, 2);

        }
    } // computeDeterministicAllelicSeqCovAndVar

    /**
     * Compute sequencing coverage likelihood per allele under negative-binomial distribution.
     *
     * @param seqCov          sequencing coverage for current taxa under current pattern
     * @param seqAlleleNum    number of sequenced alleles
     * @param taxonIndex      index of taxa
     * @param allelicSeqIndex index to allelicSeqCovPerPattern and allelicSeqCovRawVarPerPattern ( = matrixIndex * nrOfPatterns + patternIndex)
     */
    private double computeSeqCovLikelihoodPerAllele(
            final int seqCov,
            final int seqAlleleNum,
            final int taxonIndex,
            final int allelicSeqIndex
    ) {
        final double mean = (seqAlleleNum + numericalStabilizer) * (inDebugMode ? storedAllelicSeqCovPerPattern[allelicSeqIndex] : allelicSeqCovPerPattern[allelicSeqIndex]) * sizeFactors[taxonIndex];
        final double var = mean + Math.pow(sizeFactors[taxonIndex], 2) * Math.pow(seqAlleleNum + numericalStabilizer, 2) * (inDebugMode ? storedAllelicSeqCovRawVarPerPattern[allelicSeqIndex] : allelicSeqCovRawVarPerPattern[allelicSeqIndex]);
        final double p = mean / var;
        final double r = Math.pow(mean, 2) / (var - mean);

        return useLogPartials ? seqCovDensity.logDensity(seqCov, r, p) : seqCovDensity.density(seqCov, r, p);
    } // computeSeqCovLikelihoodPerAllele

}
