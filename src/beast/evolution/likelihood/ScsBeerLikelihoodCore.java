package beast.evolution.likelihood;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.variantsinfo.GenericVariantsInfo;
import beast.math.util.MathFunctions;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

@Description("Implementation of sum-product algorithm and max-sum algorithm")
public class ScsBeerLikelihoodCore extends BeerLikelihoodCore {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    protected int nrOfLeafNodes;
    protected int nrOfInternalNodes;

    /**
     * store the maximum likelihood messages sending from bottom to top using
     * max-product algorithm (actually max-sum algorithm to avoid numerical problems)
     * refer to section 8.4.5 in Bishop's "Pattern Recognition and Machine Learning"
     */
    protected double[][][] MLPartials;

    /**
     * store the partial likelihood for constant site
     * 2 * #internal nodes * [#matrices * #patterns * #states]
     */
    protected double[][][] constPartials;

    /**
     * store the scaling factor for constant site, if using scaling
     * 2 * #internal nodes * #patterns
     */
    protected double[][][] constScalingFactors;

    /**
     * 2 * [nrOfNodes]
     * only for scaling
     * for each node
     * key: pattern
     * value: original scaling factor (before taking logarithm)
     * this is because Math.exp() is an approximated method which has numerical issues when dealing with very small
     * numbers (< 1.0E-100)
     */
    private Map<Integer, Double>[][] scalingFactorMap;

    /**
     * 2 * [nrOfInternalNodes]
     * only for scaling
     * for each internal node
     * key: pattern
     * value: original scaling factor (before taking logarithm)
     * this is because Math.exp() is an approximated method which has numerical issues when dealing with very small
     * numbers (< 1.0E-100)
     */
    private Map<Integer, Double>[][] constScalingFactorMap;

    final private double scalingThreshold = 1.0E-100;

    protected boolean useLogPartials;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public ScsBeerLikelihoodCore(int nrOfStates) {
        super(nrOfStates);
    }

    /**
     * initialize partial likelihood arrays
     *
     * @param nodeCount           the number of nodes in the tree
     * @param leafNodeCount       the number of external nodes in the tree
     * @param internalNodeCount   the number of internal nodes in the tree
     * @param patternCount        the number of patterns of the input data
     * @param matrixCount         the number of matrices (i.e., the number of categories)
     * @param stateCount          the number of states in the evolutionary model
     * @param integrateCategories whether sites are being integrated over all matrices
     * @param useLogPartials      whether to use log-partials or not
     */
    public void initialize(int nodeCount, int leafNodeCount, int internalNodeCount, int patternCount, int matrixCount,
                           int stateCount, boolean integrateCategories, boolean useLogPartials) {
        this.nrOfNodes = nodeCount;
        this.nrOfLeafNodes = leafNodeCount;
        this.nrOfInternalNodes = internalNodeCount;
        this.nrOfPatterns = patternCount;
        this.nrOfMatrices = matrixCount;
        this.nrOfStates = stateCount;
        this.integrateCategories = integrateCategories;
        this.useLogPartials = useLogPartials;

        if (integrateCategories)
            partialsSize = matrixCount * patternCount * nrOfStates;
        else
            partialsSize = patternCount * nrOfStates;

        matrixSize = nrOfStates * nrOfStates;
        matrices = new double[2][nrOfNodes - 1][matrixCount * matrixSize];
        currentMatrixIndex = new int[nrOfNodes - 1];
        storedMatrixIndex = new int[nrOfNodes - 1];

        // partial likelihood for internal nodes
        partials = new double[2][nrOfNodes][];
        MLPartials = new double[2][internalNodeCount][];
        currentPartialsIndex = new int[nrOfNodes];
        storedPartialsIndex = new int[nrOfNodes];
        for (int i = 0; i < nrOfNodes; i++) {
            partials[0][i] = null;
            partials[1][i] = null;
            if (i < internalNodeCount) {
                MLPartials[0][i] = null;
                MLPartials[1][i] = null;
            }
        }

        // partial likelihood for constant site (ascertainment bias correction)
        constPartials = new double[2][internalNodeCount][partialsSize];
    } // initialize

    /**
     * Allocates partials for an internal node
     */
    @Override
    public void createNodePartials(int nodeIndex) {
        this.partials[0][nodeIndex] = new double[partialsSize];
        this.partials[1][nodeIndex] = new double[partialsSize];

        if (nodeIndex >= nrOfLeafNodes) {
            this.MLPartials[0][nodeIndex - nrOfLeafNodes] = new double[partialsSize];
            this.MLPartials[1][nodeIndex - nrOfLeafNodes] = new double[partialsSize];
        }
    } // createInternalNodePartials

    /**
     * integrate likelihoods for each pattern across all site categories
     *
     * @param nodeIndex    root node index
     * @param proportions  proportions of site categories
     * @param rootGenotype genotype of the root node
     * @param outPartials  #patterns, integrated likelihoods for each pattern
     * @param constRoot    root likelihood for constant site (passed by reference)
     */
    public void integratePartials(int nodeIndex, double[] proportions, int rootGenotype, double[] outPartials, double[] constRoot) {
        calculateIntegratePartials(
                partials[currentPartialsIndex[nodeIndex]][nodeIndex],
                proportions,
                rootGenotype,
                outPartials,
                constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes],
                constRoot
        );
    } // integratePartials

    /**
     * Integrates partials (both variants and constant) across categories.
     *
     * @param inPartials    the array of partials to be integrated
     * @param proportions   the proportions of sites in each category
     * @param rootGenotype  genotype of the root node
     * @param outPartials   an array into which the partials will go
     * @param constPartials the arrays of partials of constant site
     * @param constRoot     root likelihood for constant site (passed by reference)
     */
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, int rootGenotype,
                                              double[] outPartials, double[] constPartials, double[] constRoot) {
        if (useLogPartials) {
            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                int inIndex = patternIndex * nrOfStates + rootGenotype;

                double[] out = new double[nrOfMatrices];
                double[] cst = new double[nrOfMatrices];

                for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {
                    out[matrixIndex] = Math.log(proportions[matrixIndex]) + inPartials[inIndex];
                    cst[matrixIndex] = Math.log(proportions[matrixIndex]) + constPartials[inIndex];

                    inIndex += nrOfPatterns * nrOfStates;
                }

                outPartials[patternIndex] = MathFunctions.logSumExp(out);
                constRoot[patternIndex] = MathFunctions.logSumExp(cst);
            }
        } else {
            int inIndex = rootGenotype;

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                outPartials[patternIndex] = inPartials[inIndex] * proportions[0];
                constRoot[patternIndex] = constPartials[inIndex] * proportions[0];

                inIndex += nrOfStates;
            }

            for (int matrixIndex = 1; matrixIndex < nrOfMatrices; matrixIndex++) {
                for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                    outPartials[patternIndex] += inPartials[inIndex] * proportions[matrixIndex];
                    constRoot[patternIndex] += constPartials[inIndex] * proportions[matrixIndex];

                    inIndex += nrOfStates;
                }
            }
        }
    } // calculateIntegratePartials

    /**
     * Integrates partials (only variants) across categories.
     *
     * @param inPartials   the array of partials to be integrated
     * @param proportions  the proportions of sites in each category
     * @param rootGenotype genotype of the root node
     * @param outPartials  an array into which the partials will go
     */
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, int rootGenotype,
                                              double[] outPartials) {
        int inIndex = rootGenotype;

        for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
            outPartials[patternIndex] = inPartials[inIndex] * proportions[0];

            inIndex += nrOfStates;
        }

        for (int matrixIndex = 1; matrixIndex < nrOfMatrices; matrixIndex++) {
            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
                outPartials[patternIndex] += inPartials[inIndex] * proportions[matrixIndex];

                inIndex += nrOfStates;
            }
        }
    } // calculateIntegratePartials

    /**
     * Compute log likelihood for each pattern.
     *
     * @param partials             #patterns, integrated root likelihoods
     * @param outLogLikelihoods    #patterns, output log likelihoods
     * @param constRoot            root likelihood for constant site (passed by reference)
     * @param logConstRoot         (scaled) root log-likelihood for constant site (passed by reference)
     */
    public void calculateLogLikelihoods(double[] partials, double[] outLogLikelihoods, double[] constRoot,
                                        double[] logConstRoot) {
        for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {
            outLogLikelihoods[patternIndex] = getLogScalingFactor(patternIndex);
            logConstRoot[patternIndex] = getConstLogScalingFactor(patternIndex);

            if (useLogPartials) {
                outLogLikelihoods[patternIndex] += partials[patternIndex];
                logConstRoot[patternIndex] += constRoot[patternIndex];
            } else {
                outLogLikelihoods[patternIndex] += Math.log(partials[patternIndex]);
                logConstRoot[patternIndex] += Math.log(constRoot[patternIndex]);
            }
        }
    } // calculateLogLikelihoods

    /**
     * to which category each pattern belongs?
     *
     * @param rootGenotype genotype of root node
     * @param root         root node index
     * @param out          output
     */
    public void getPatternCategories(int rootGenotype, int root, GenericVariantsInfo.Base out) {
        for (int i = 0; i < nrOfPatterns; i++) {
            List<Integer> matrices = new ArrayList<>();
            double max = -Double.MAX_VALUE;

            for (int j = 0; j < nrOfMatrices; j++) {
                int index = j * nrOfPatterns * nrOfStates + i * nrOfStates + rootGenotype;
                double tmp = partials[currentPartialsIndex[root]][root][index];

                if (max < tmp) {
                    max = tmp;
                    matrices.clear();
                    matrices.add(j);
                } else if (max == tmp) {
                    matrices.add(j);
                }
            }

            out.addMLCategory(i, matrices, true);
        }
    } // getPatternCategories

    public double getLogTransferProbability(
            Node node,
            int matrix,
            int parentState,
            int childState
    ) {
        if (node.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        } else {
            final int nodeIndex = node.getNr();

            return Math.log(matrices[currentMatrixIndex[nodeIndex]][nodeIndex][matrix * this.matrixSize + parentState * this.nrOfStates + childState]);
        }
    } // getLogTransferProbability

    public double[] getLogTransferProbabilityAll(
            Node node,
            int matrix,
            int parentState
    ) {
        if (node.isRoot()) {
            return null;
        } else {
            double[] results = new double[this.nrOfStates];

            for (int i = 0; i < this.nrOfStates; i++) {
                results[i] = getLogTransferProbability(
                        node,
                        matrix,
                        parentState,
                        i
                );
            }

            return results;
        }
    } // getLogTransferProbabilityAll

    public double getLogLikelihood(
            Node node,
            int matrix,
            int pattern,
            int state
    ) {
        int index = node.getNr();

        final double value = partials[currentPartialsIndex[index]][index][matrix * this.nrOfPatterns * this.nrOfStates + pattern * this.nrOfStates + state];

        return this.useLogPartials ? value : Math.log(value);
    } // getLikelihood

    public double[] getLogLikelihoodsAllStates(
            Node node,
            int matrix,
            int pattern
    ) {
        double[] results = new double[this.nrOfStates];

        for (int i = 0; i < this.nrOfStates; i++)
            results[i] = getLogLikelihood(node, matrix, pattern, i);

        return results;
    } // getLogLikelihoodsAllStates

    public double getLogLikelihoodConstantPattern(
            Node node,
            int matrix,
            int pattern,
            int state
    ) {
        int index = node.getNr();

        final double value = constPartials[currentPartialsIndex[index]][index - nrOfLeafNodes][matrix * this.nrOfPatterns * this.nrOfStates + pattern * this.nrOfStates + state];

        return this.useLogPartials ? value : Math.log(value);
    } // getLogLikelihoodConstantPattern

    /**
     * set up scaling
     *
     * @param scale value to compare
     */
    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][nrOfNodes][nrOfPatterns];
            scalingFactorMap = new HashMap[2][nrOfNodes];

            constScalingFactors = new double[2][nrOfInternalNodes][nrOfPatterns];
            constScalingFactorMap = new HashMap[2][nrOfInternalNodes];
        }
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws java.lang.Throwable {
        super.finalize();

        nrOfLeafNodes = 0;
        nrOfInternalNodes = 0;

        MLPartials = null;
        constPartials = null;
        constScalingFactors = null;
    } // finalize

    /**
     * Sets partials for a node
     */
    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {
        if (this.partials[0][nodeIndex] == null)
            createNodePartials(nodeIndex);

        if (partials.length < partialsSize) {
            int k = 0;
            for (int i = 0; i < nrOfMatrices; i++) {
                System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], k, partials.length);
                k += partials.length;
            }
        } else
            System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partials.length);

        if (useScaling)
            scalePartials(nodeIndex);
    } // setNodePartials

    /**
     * Duplicate leaf partials to the alternative position.
     * Should ONLY be called during initialization.
     *
     * @param nodeIndex apaprently
     */
    public void storeLeafPartials(int nodeIndex) {
        System.arraycopy(
                this.partials[currentPartialsIndex[nodeIndex]][nodeIndex],
                0,
                this.partials[1 - currentPartialsIndex[nodeIndex]][nodeIndex],
                0,
                partialsSize
        );
    } // storeLeafPartials

    /**
     * Sets partials for a node
     *
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> setNodePartials(int nodeIndex, double[] partials, final int[][] changedPatterns, final List<Integer> changedPatternsIndex) {

        if (this.partials[0][nodeIndex] == null) {
            createNodePartials(nodeIndex);
        }
        if (partials.length < partialsSize) {
            int k = 0;
            for (int i = 0; i < nrOfMatrices; i++) {
                System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], k, partials.length);
                k += partials.length;
            }
        } else {
            System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partials.length);
        }

        if (useScaling) {
            return scalePartials(nodeIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // setNodePartials

    @Override
    protected void scalePartials(int nodeIndex) {
        int u = 0;

        for (int i = 0; i < nrOfPatterns; i++) {

            double scaleFactor, constScaleFactor;
            scaleFactor = constScaleFactor = 0.0;
            int v = u;
            for (int k = 0; k < nrOfMatrices; k++) {
                for (int j = 0; j < nrOfStates; j++) {
                    if (partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] > scaleFactor) {
                        scaleFactor = partials[currentPartialsIndex[nodeIndex]][nodeIndex][v];
                    }
                    if (nodeIndex >= nrOfLeafNodes && constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][v] > constScaleFactor) {
                        constScaleFactor = constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][v];
                    }
                    v++;
                }
                v += (nrOfPatterns - 1) * nrOfStates;
            }

            if (scaleFactor < scalingThreshold) {

                if (scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex] == null) {
                    scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex] = new HashMap<>();
                }

                scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex].put(i, scaleFactor);

                v = u;
                for (int k = 0; k < nrOfMatrices; k++) {
                    for (int j = 0; j < nrOfStates; j++) {
                        if (useLogPartials) {
                            partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] -= scaleFactor;
                        } else {
                            partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] /= scaleFactor;
                        }

                        v++;
                    }
                    v += (nrOfPatterns - 1) * nrOfStates;
                }
                if (useLogPartials) {
                    scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = scaleFactor;
                } else {
                    scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = Math.log(scaleFactor);
                }

            } else {

                if (scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex] != null) {
                    scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex].remove(i);
                }

                scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = 0.0;
            }

            if (nodeIndex >= nrOfLeafNodes) {
                if (constScaleFactor < scalingThreshold) {

                    if (constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes] == null) {
                        constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes] = new HashMap<>();
                    }

                    constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes].put(i, scaleFactor);

                    v = u;
                    for (int k = 0; k < nrOfMatrices; k++) {
                        for (int j = 0; j < nrOfStates; j++) {
                            if (useLogPartials) {
                                constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][v] -= constScaleFactor;
                            } else {
                                constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][v] /= constScaleFactor;
                            }

                            v++;
                        }
                        v += (nrOfPatterns - 1) * nrOfStates;
                    }
                    if (useLogPartials) {
                        constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][i] = constScaleFactor;
                    } else {
                        constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][i] = Math.log(constScaleFactor);
                    }

                } else {

                    if (constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes] != null) {
                        constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes].remove(i);
                    }

                    constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][i] = 0.0;
                }
            }
            u += nrOfStates;

        }
    }
    /**
     * should only be called during post process when only some of the matrices and patterns are updated
     * <p>
     * changedPatterns could be updated in those cases that scaleFactor of a pattern is changed and not all matrices
     * for this pattern are in changedPatterns before
     *
     * @param nodeIndex            which node is being processed
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    protected List<int[]> scalePartials(final int nodeIndex, final int[][] changedPatterns,
                                        final List<Integer> changedPatternsIndex) {
        List<int[]> added = new ArrayList<>();

        // loop over all changed patterns
        for (int i = 0; i < changedPatternsIndex.size(); i++) {

            final int patternIndex = changedPatterns[changedPatternsIndex.get(i)][1];

            boolean scaledBefore, constScaledBefore = false, newItemsAdded = false;

            // prepare CHANGED matrices
            List<Integer> changedMatrices = new ArrayList<>();
            final int nextPIndex = i == changedPatternsIndex.size() - 1 ? changedPatterns.length : changedPatternsIndex.get(i + 1);
            for (int j = changedPatternsIndex.get(i); j < nextPIndex; j++) {
                changedMatrices.add(changedPatterns[j][0]);
            }

            // if current node was scaled at this pattern before, then recover those UNCHANGED matrices and genotypes
            // for normal
            scaledBefore = recoverUnchangedForPattern(patternIndex, changedMatrices,
                    scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][patternIndex],
                    scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex],
                    partials[currentPartialsIndex[nodeIndex]][nodeIndex]);
            // for constant
            if (nodeIndex >= nrOfLeafNodes) {
                constScaledBefore = recoverUnchangedForPattern(patternIndex, changedMatrices,
                        constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][patternIndex],
                        constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes],
                        constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes]);
            }

            // get scale factor
            double scaleFactor, constScaleFactor;
            scaleFactor = constScaleFactor = 0.0;
            int index = patternIndex * nrOfStates;
            for (int j = 0; j < nrOfMatrices; j++) {
                for (int k = 0; k < nrOfStates; k++) {
                    if (partials[currentPartialsIndex[nodeIndex]][nodeIndex][index] > scaleFactor) {
                        scaleFactor = partials[currentPartialsIndex[nodeIndex]][nodeIndex][index];
                    }
                    if (nodeIndex >= nrOfLeafNodes && constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][index] > constScaleFactor) {
                        constScaleFactor = constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][index];
                    }
                    index++;
                }
                index += (nrOfPatterns - 1) * nrOfStates;
            }

            // scaling for normal
            if (scaleFactor < scalingThreshold) {

                if (!scaledBefore || scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex].get(patternIndex) != scaleFactor) {
                    // if current pattern was not scaled before
                    // or
                    // if current pattern was scaled before and if current scale factor does not equal to the previous one
                    // add all UNCHANGED matrices of current pattern to the changedPatterns array

                    newItemsAdded = true;

                    for (int j = 0; j < nrOfMatrices; j++) {
                        if (!changedMatrices.contains(j)) {
                            added.add(IntStream.of(j, patternIndex).toArray());
                        }
                    }
                }

                // update key-value pair in scalingFactorMap
                if (scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex] == null) {
                    scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex] = new HashMap<>();
                }
                scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex].put(patternIndex, scaleFactor);

                // scale partials
                index = patternIndex * nrOfStates;
                for (int j = 0; j < nrOfMatrices; j++) {
                    for (int k = 0; k < nrOfStates; k++) {
                        if (useLogPartials) {
                            partials[currentPartialsIndex[nodeIndex]][nodeIndex][index] -= scaleFactor;
                        } else {
                            partials[currentPartialsIndex[nodeIndex]][nodeIndex][index] /= scaleFactor;
                        }

                        index++;
                    }
                    index += (nrOfPatterns - 1) * nrOfStates;
                }

                if (useLogPartials) {
                    scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][patternIndex] = scaleFactor;
                } else {
                    scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][patternIndex] = Math.log(scaleFactor);
                }

            } else {

                // if current pattern was scaled before, add all UNCHANGED matrices of current pattern to the
                // changedPatterns array
                if (scaledBefore) {
                    newItemsAdded = true;

                    for (int j = 0; j < nrOfMatrices; j++) {
                        if (!changedMatrices.contains(j)) {
                            added.add(IntStream.of(j, patternIndex).toArray());
                        }
                    }
                }

                // update key-value pair in scalingFactorMap
                if (scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex] != null) {
                    scalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex].remove(patternIndex);
                }

                scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][patternIndex] = 0.0;

            }

            // scaling for constant
            if (nodeIndex >= nrOfLeafNodes) {
                if (constScaleFactor > scalingThreshold) {

                    if (!newItemsAdded && (!constScaledBefore || constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes].get(patternIndex) != constScaleFactor)) {
                        // if UNCHANGED matrices are not added before
                        // and :
                        // if current pattern was not scaled before
                        // or
                        // if current pattern was scaled before and if current scale factor does not equal to the previous one
                        // add all UNCHANGED matrices of current pattern to the changedPatterns array

                        for (int j = 0; j < nrOfMatrices; j++) {
                            if (!changedMatrices.contains(j)) {
                                added.add(IntStream.of(j, patternIndex).toArray());
                            }
                        }
                    }

                    // update key-value pair in scalingFactorMap
                    if (constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes] == null) {
                        constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes] = new HashMap<>();
                    }
                    constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes].put(patternIndex, constScaleFactor);

                    // scale constant partials
                    index = patternIndex * nrOfStates;
                    for (int j = 0; j < nrOfMatrices; j++) {
                        for (int k = 0; k < nrOfStates; k++) {
                            if (useLogPartials) {
                                constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][index] -= constScaleFactor;
                            } else {
                                constPartials[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][index] /= constScaleFactor;
                            }

                            index++;
                        }
                        index += (nrOfPatterns - 1) * nrOfStates;
                    }

                    if (useLogPartials) {
                        constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][patternIndex] = constScaleFactor;
                    } else {
                        constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][patternIndex] = Math.log(constScaleFactor);
                    }

                } else {

                    // if UNCHANGED matrices are not added before
                    // and
                    // if current pattern was scaled before
                    // add all UNCHANGED matrices of current pattern to the changedPatterns array
                    if (!newItemsAdded && constScaledBefore) {
                        for (int j = 0; j < nrOfMatrices; j++) {
                            if (!changedMatrices.contains(j)) {
                                added.add(IntStream.of(j, patternIndex).toArray());
                            }
                        }
                    }

                    // update key-value pair in scalingFactorMap
                    if (constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes] != null) {
                        constScalingFactorMap[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes].remove(patternIndex);
                    }

                    constScalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex - nrOfLeafNodes][patternIndex] = 0.0;

                }
            }

        }

        if (added.size() > 0) {
            return added;
        } else {
            return null;
        }
    }

    /**
     * should only be called when scaling partials at post process stage
     * this method is ESSENTIAL to accurately compute likelihoods
     *
     * @param patternIndex    which pattern?
     * @param changedMatrices matrices listed in changed patterns already
     * @param scalingFactor   used to scale partials previously (log)
     * @param factorMap       storing original scale factor (before log)
     * @param partial         partial likelihoods
     * @return partials recovered or not?
     */
    private boolean recoverUnchangedForPattern(final int patternIndex, final List<Integer> changedMatrices,
                                               final double scalingFactor, final Map<Integer, Double> factorMap,
                                               double[] partial) {
        if (scalingFactor != 0.0) {
            assert factorMap != null;

            final double scaleFactor = factorMap.get(patternIndex);

            // recover those UNCHANGED matrices
            for (int i = 0; i < nrOfMatrices; i++) {
                if (!changedMatrices.contains(i)) {
                    final int pIndex = i * nrOfPatterns * nrOfStates + patternIndex * nrOfStates;
                    for (int j = 0; j < nrOfStates; j++) {
                        if (useLogPartials) {
                            partial[pIndex + j] += scaleFactor;
                        } else {
                            partial[pIndex + j] *= scaleFactor;
                        }
                    }
                }
            }

            return true;
        } else {
            return false;
        }
    } // recoverUnchangedForPattern

    /**
     * This function returns the scaling factor for constant site by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     *
     * @param patternIndex_ which pattern?
     * @return the log scaling factor
     */
    public double getConstLogScalingFactor(int patternIndex_) {
        double logScalingFactor = 0.0;
        if (useScaling) {
            for (int i = nrOfLeafNodes; i < nrOfNodes; i++) {
                logScalingFactor += constScalingFactors[currentPartialsIndex[i]][i - nrOfLeafNodes][patternIndex_];
            }
        }

        return logScalingFactor;
    } // getConstLogScalingFactor


    //****************************************************************
    //*            sum-product algorithm (normal partial)            *
    //****************************************************************

    /**
     * calculate partial likelihoods at a node
     * when the node is tree root, childIndex2 is set -1
     *
     * @param childIndex1   index of the first child
     * @param isLeaf1       is the first child is a leaf node?
     * @param childIndex2   index of the second child
     * @param isLeaf2       is the second child is a leaf node?
     * @param parentIndex   index of parent node
     * @param constGenotype genotype of constant site
     */
    public void calculatePartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            final int constGenotype
    ) {
        if (isLeaf1) {
            assert childIndex2 >= 0 : "two children required but only one provided";

            if (isLeaf2) {
                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype
                );
            } else {
                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype
                );
            }
        } else {
            if (isLeaf2) {
                assert childIndex2 >= 0 : "two children required but only one provided";

                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype
                );
            } else {
                if (childIndex2 >= 0) {
                    calculatePartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            partials[currentPartialsIndex[childIndex2]][childIndex2],
                            constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex2]][childIndex2],
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype
                    );
                } else if (childIndex2 == -1) {
                    calculatePartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            null,
                            null,
                            null,
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype
                    );
                } else {
                    throw new IllegalArgumentException("childIndex2 is only allowed to be a integer no smaller than -1 (" +
                            this.getClass().getName() + ")");
                }
            }
        }

        if (useScaling) {
            scalePartials(parentIndex);
        }
    } // calculatePartials

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param childPartialsIndex1 #matrices * #patterns * #states
     * @param childConstPartialsIndex1 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex1      #matrices * #states * #states
     * @param childPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param childConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2      #matrices * #states * #states, allowed to be null
     * @param parentPartialsIndex #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param constGenotype genotype of constant site
     */
    protected void calculatePartialPartialPruning(
            final double[] childPartialsIndex1,
            final double[] childConstPartialsIndex1,
            final double[] matricesIndex1,
            final double[] childPartialsIndex2,
            final double[] childConstPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            final int constGenotype
    ) {
        if ((childPartialsIndex2 == null && matricesIndex2 != null) ||
                (childPartialsIndex2 != null && matricesIndex2 == null)) {
            throw new IllegalArgumentException("childPartialsIndex2 and matricesIndex2 should be defined or be null " +
                    "synchronously (" + this.getClass().getName() + ")");
        }

        if (childPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("childPartialsIndex1 and matricesIndex1 should be defined instead of " +
                    "being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (childPartialsIndex2 != null && matricesIndex2 != null);

        double cst1, cst2;
        double sum1, sum2;
        double tmp1, tmp2;
        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    cst1 = cst2 = 0.0;
                    sum1 = sum2 = 0.0;
                    tmp1 = tmp2 = 0.0;

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        tmp1 = matricesIndex1[mIndex] * childPartialsIndex1[cIndex + cGenotypeIndex];

                        if (childConstPartialsIndex1 == null) {
                            // leaf

                            if (cGenotypeIndex == constGenotype) {
                                if (childPartialsIndex1[cIndex + cGenotypeIndex] <= 0.0) {
                                    cst1 = matricesIndex1[mIndex] * Double.MIN_VALUE;
                                } else {
                                    cst1 = tmp1;
                                }
                            }
                        } else {
                            // internal node

                            cst1 += matricesIndex1[mIndex] * childConstPartialsIndex1[cIndex + cGenotypeIndex];
                        }

                        sum1 += tmp1;

                        if (has2ndChild) {

                            tmp2 = matricesIndex2[mIndex] * childPartialsIndex2[cIndex + cGenotypeIndex];

                            if (childConstPartialsIndex2 == null) {
                                // leaf

                                if (cGenotypeIndex == constGenotype) {
                                    if (childPartialsIndex2[cIndex + cGenotypeIndex] <= 0.0) {
                                        cst2 = matricesIndex2[mIndex] * Double.MIN_VALUE;
                                    } else {
                                        cst2 = tmp2;
                                    }
                                }
                            } else {
                                // internal node

                                cst2 += matricesIndex2[mIndex] * childConstPartialsIndex2[cIndex + cGenotypeIndex];
                            }

                            sum2 += tmp2;

                        }

                        mIndex++;
                    }

                    if (has2ndChild) {
                        parentConstPartialsIndex[pIndex] = cst1 * cst2;
                        parentPartialsIndex[pIndex] = sum1 * sum2;
                    } else {
                        parentConstPartialsIndex[pIndex] = cst1;
                        parentPartialsIndex[pIndex] = sum1;
                    }

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculatePartialPartialPruning

    /**
     * calculate partial likelihoods at a node
     * when the node is tree root, childIndex2 is set -1
     *
     * @param childIndex1          index of the first child
     * @param isLeaf1 is the first child is a leaf node?
     * @param childIndex2          index of the second child
     * @param isLeaf2 is the second child is a leaf node?
     * @param parentIndex          index of parent node
     * @param constGenotype genotype of constant site
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> calculatePartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            final int constGenotype,
            final int[][] changedPatterns,
            final List<Integer> changedPatternsIndex
    ) {
        if (isLeaf1) {
            assert childIndex2 >= 0 : "two children required but only one provided";

            if (isLeaf2) {
                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype,
                        changedPatterns
                );
            } else {
                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype,
                        changedPatterns
                );
            }
        } else {
            if (isLeaf2) {
                assert childIndex2 >= 0 : "two children required but only one provided";

                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype,
                        changedPatterns
                );
            } else {
                if (childIndex2 >= 0) {
                    calculatePartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            partials[currentPartialsIndex[childIndex2]][childIndex2],
                            constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex2]][childIndex2],
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype,
                            changedPatterns
                    );
                } else if (childIndex2 == -1) {
                    calculatePartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            null,
                            null,
                            null,
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype,
                            changedPatterns
                    );
                } else {
                    throw new IllegalArgumentException("childIndex2 is only allowed to be a integer no smaller than -1 (" +
                            this.getClass().getName() + ")");
                }
            }
        }

        if (useScaling) {
            return scalePartials(parentIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // calculatePartials

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param childPartialsIndex1 #matrices * #patterns * #states
     * @param childConstPartialsIndex1 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex1      #matrices * #states * #states
     * @param childPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param childConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2      #matrices * #states * #states, allowed to be null
     * @param parentPartialsIndex #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param constGenotype genotype of constant site
     * @param changedPatterns     patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculatePartialPartialPruning(
            final double[] childPartialsIndex1,
            final double[] childConstPartialsIndex1,
            final double[] matricesIndex1,
            final double[] childPartialsIndex2,
            final double[] childConstPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            final int constGenotype,
            final int[][] changedPatterns
    ) {
        if ((childPartialsIndex2 == null && matricesIndex2 != null) ||
                (childPartialsIndex2 != null && matricesIndex2 == null)) {
            throw new IllegalArgumentException("childPartialsIndex2 and matricesIndex2 should be defined or be null " +
                    "synchronously (" + this.getClass().getName() + ")");
        }

        if (childPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("childPartialsIndex1 and matricesIndex1 should be defined instead of " +
                    "being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (childPartialsIndex2 != null && matricesIndex2 != null);

        double cst1, cst2;
        double sum1, sum2;
        double tmp1, tmp2;
        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                cst1 = cst2 = 0.0;
                sum1 = sum2 = 0.0;
                tmp1 = tmp2 = 0.0;

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    tmp1 = matricesIndex1[mIndex] * childPartialsIndex1[cIndex + cGenotypeIndex];

                    if (childConstPartialsIndex1 == null) {
                        if (cGenotypeIndex == constGenotype) {
                            if (childPartialsIndex1[cIndex + cGenotypeIndex] <= 0.0) {
                                cst1 = matricesIndex1[mIndex] * Double.MIN_VALUE;
                            } else {
                                cst1 = tmp1;
                            }
                        }
                    } else {
                        cst1 += matricesIndex1[mIndex] * childConstPartialsIndex1[cIndex + cGenotypeIndex];
                    }

                    sum1 += tmp1;

                    if (has2ndChild) {

                        tmp2 = matricesIndex2[mIndex] * childPartialsIndex2[cIndex + cGenotypeIndex];

                        if (childConstPartialsIndex2 == null) {
                            if (cGenotypeIndex == constGenotype) {
                                if (childPartialsIndex2[cIndex + cGenotypeIndex] <= 0.0) {
                                    cst2 = matricesIndex2[mIndex] * Double.MIN_VALUE;
                                } else {
                                    cst2 = tmp2;
                                }
                            }
                        } else {
                            cst2 += matricesIndex2[mIndex] * childConstPartialsIndex2[cIndex + cGenotypeIndex];
                        }

                        sum2 += tmp2;

                    }

                    mIndex++;
                }

                if (has2ndChild) {
                    parentConstPartialsIndex[pIndex] = cst1 * cst2;
                    parentPartialsIndex[pIndex] = sum1 * sum2;
                } else {
                    parentConstPartialsIndex[pIndex] = cst1;
                    parentPartialsIndex[pIndex] = sum1;
                }

                pIndex++;
            }
        }
    } // calculatePartialPartialPruning

    //*************************************************************
    //*            sum-product algorithm (log partial)            *
    //*************************************************************

    /**
     * calculate partial likelihoods at a node
     * when the node is tree root, childIndex2 is set -1
     *
     * @param childIndex1   index of the first child
     * @param isLeaf1       is the first child is a leaf node?
     * @param childIndex2   index of the second child
     * @param isLeaf2       is the second child is a leaf node?
     * @param parentIndex   index of parent node
     * @param constGenotype genotype of constant site
     */
    public void calculateLogPartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            final int constGenotype
    ) {
        if (isLeaf1) {
            assert childIndex2 >= 0 : "two children required but only one provided";

            if (isLeaf2) {
                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype
                );
            } else {
                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype
                );
            }
        } else {
            if (isLeaf2) {
                assert childIndex2 >= 0 : "two children required but only one provided";

                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype
                );
            } else {
                if (childIndex2 >= 0) {
                    calculateLogPartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            partials[currentPartialsIndex[childIndex2]][childIndex2],
                            constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex2]][childIndex2],
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype
                    );
                } else if (childIndex2 == -1) {
                    calculateLogPartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            null,
                            null,
                            null,
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype
                    );
                } else {
                    throw new IllegalArgumentException("childIndex2 is only allowed to be a integer no smaller than -1 (" +
                            this.getClass().getName() + ")");
                }
            }
        }

        if (useScaling) {
            scalePartials(parentIndex);
        }
    } // calculateLogPartials

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param childPartialsIndex1      #matrices * #patterns * #states
     * @param childConstPartialsIndex1 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex1           #matrices * #states * #states
     * @param childPartialsIndex2      #matrices * #patterns * #states, allowed to be null
     * @param childConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2           #matrices * #states * #states, allowed to be null
     * @param parentPartialsIndex      #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param constGenotype            genotype of constant site
     */
    protected void calculateLogPartialPartialPruning(
            final double[] childPartialsIndex1,
            final double[] childConstPartialsIndex1,
            final double[] matricesIndex1,
            final double[] childPartialsIndex2,
            final double[] childConstPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            final int constGenotype
    ) {
        if ((childPartialsIndex2 == null && matricesIndex2 != null) ||
                (childPartialsIndex2 != null && matricesIndex2 == null)) {
            throw new IllegalArgumentException("childPartialsIndex2 and matricesIndex2 should be defined or be null " +
                    "synchronously (" + this.getClass().getName() + ")");
        }

        if (childPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("childPartialsIndex1 and matricesIndex1 should be defined instead of " +
                    "being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (childPartialsIndex2 != null && matricesIndex2 != null);

        double cst1, cst2; // for leaf
        double[] cstArr1 = new double[nrOfStates]; // for internal node
        double[] cstArr2 = new double[nrOfStates]; // for internal node
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];
        double tmp1, tmp2;
        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    cst1 = cst2 = 0.0;
                    tmp1 = tmp2 = 0.0;

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        tmp1 = Math.log(matricesIndex1[mIndex]) + childPartialsIndex1[cIndex + cGenotypeIndex];

                        if (childConstPartialsIndex1 == null) {
                            // leaf

                            if (cGenotypeIndex == constGenotype)
                                cst1 = tmp1;
                        } else {
                            // internal node

                            cstArr1[cGenotypeIndex] = Math.log(matricesIndex1[mIndex]) + childConstPartialsIndex1[cIndex + cGenotypeIndex];
                        }

                        sp1[cGenotypeIndex] = tmp1;

                        if (has2ndChild) {

                            tmp2 = Math.log(matricesIndex2[mIndex]) + childPartialsIndex2[cIndex + cGenotypeIndex];

                            if (childConstPartialsIndex2 == null) {
                                // leaf

                                if (cGenotypeIndex == constGenotype)
                                    cst2 = tmp2;
                            } else {
                                // internal node

                                cstArr2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + childConstPartialsIndex2[cIndex + cGenotypeIndex];
                            }

                            sp2[cGenotypeIndex] = tmp2;

                        }

                        mIndex++;
                    }

                    if (has2ndChild) {

                        if (childConstPartialsIndex1 == null) {
                            if (childConstPartialsIndex2 == null) {
                                parentConstPartialsIndex[pIndex] = cst1 + cst2;
                            } else {
                                parentConstPartialsIndex[pIndex] = cst1 + MathFunctions.logSumExp(cstArr2);
                            }
                        } else {
                            if (childConstPartialsIndex2 == null) {
                                parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cstArr1) + cst2;
                            } else {
                                parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cstArr1) + MathFunctions.logSumExp(cstArr2);
                            }
                        }

                        parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                    } else {
                        if (childConstPartialsIndex1 == null) {
                            parentConstPartialsIndex[pIndex] = cst1;
                        } else {
                            parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cstArr1);
                        }

                        parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1);
                    }

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculateLogPartialPartialPruning

    /**
     * calculate partial likelihoods at a node
     * when the node is tree root, childIndex2 is set -1
     *
     * @param childIndex1          index of the first child
     * @param isLeaf1              is the first child is a leaf node?
     * @param childIndex2          index of the second child
     * @param isLeaf2              is the second child is a leaf node?
     * @param parentIndex          index of parent node
     * @param constGenotype        genotype of constant site
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> calculateLogPartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            final int constGenotype,
            final int[][] changedPatterns,
            final List<Integer> changedPatternsIndex
    ) {
        if (isLeaf1) {
            assert childIndex2 >= 0 : "two children required but only one provided";

            if (isLeaf2) {
                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype,
                        changedPatterns
                );
            } else {
                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        null,
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype,
                        changedPatterns
                );
            }
        } else {
            if (isLeaf2) {
                assert childIndex2 >= 0 : "two children required but only one provided";

                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        null,
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        constGenotype,
                        changedPatterns
                );
            } else {
                if (childIndex2 >= 0) {
                    calculateLogPartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            partials[currentPartialsIndex[childIndex2]][childIndex2],
                            constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex2]][childIndex2],
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype,
                            changedPatterns
                    );
                } else if (childIndex2 == -1) {
                    calculateLogPartialPartialPruning(
                            partials[currentPartialsIndex[childIndex1]][childIndex1],
                            constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                            matrices[currentMatrixIndex[childIndex1]][childIndex1],
                            null,
                            null,
                            null,
                            partials[currentPartialsIndex[parentIndex]][parentIndex],
                            constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                            constGenotype,
                            changedPatterns
                    );
                } else {
                    throw new IllegalArgumentException("childIndex2 is only allowed to be a integer no smaller than -1 (" +
                            this.getClass().getName() + ")");
                }
            }
        }

        if (useScaling) {
            return scalePartials(parentIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // calculateLogPartials

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param childPartialsIndex1      #matrices * #patterns * #states
     * @param childConstPartialsIndex1 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex1           #matrices * #states * #states
     * @param childPartialsIndex2      #matrices * #patterns * #states, allowed to be null
     * @param childConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2           #matrices * #states * #states, allowed to be null
     * @param parentPartialsIndex      #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param constGenotype            genotype of constant site
     * @param changedPatterns          patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculateLogPartialPartialPruning(
            final double[] childPartialsIndex1,
            final double[] childConstPartialsIndex1,
            final double[] matricesIndex1,
            final double[] childPartialsIndex2,
            final double[] childConstPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            final int constGenotype,
            final int[][] changedPatterns
    ) {
        if ((childPartialsIndex2 == null && matricesIndex2 != null) ||
                (childPartialsIndex2 != null && matricesIndex2 == null)) {
            throw new IllegalArgumentException("childPartialsIndex2 and matricesIndex2 should be defined or be null " +
                    "synchronously (" + this.getClass().getName() + ")");
        }

        if (childPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("childPartialsIndex1 and matricesIndex1 should be defined instead of " +
                    "being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (childPartialsIndex2 != null && matricesIndex2 != null);

        double cst1, cst2; // for leaf
        double[] cstArr1 = new double[nrOfStates]; // for internal node
        double[] cstArr2 = new double[nrOfStates]; // for internal node
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];
        double tmp1, tmp2;
        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                cst1 = cst2 = 0.0;
                tmp1 = tmp2 = 0.0;

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    tmp1 = Math.log(matricesIndex1[mIndex]) + childPartialsIndex1[cIndex + cGenotypeIndex];

                    if (childConstPartialsIndex1 == null) {
                        // leaf

                        if (cGenotypeIndex == constGenotype)
                            cst1 = tmp1;
                    } else {
                        // internal node

                        cstArr1[cGenotypeIndex] = Math.log(matricesIndex1[mIndex]) + childConstPartialsIndex1[cIndex + cGenotypeIndex];
                    }

                    sp1[cGenotypeIndex] = tmp1;

                    if (has2ndChild) {

                        tmp2 = Math.log(matricesIndex2[mIndex]) + childPartialsIndex2[cIndex + cGenotypeIndex];

                        if (childConstPartialsIndex2 == null) {
                            // leaf

                            if (cGenotypeIndex == constGenotype)
                                cst2 = tmp2;
                        } else {
                            // internal node

                            cstArr2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + childConstPartialsIndex2[cIndex + cGenotypeIndex];
                        }

                        sp2[cGenotypeIndex] = tmp2;

                    }

                    mIndex++;
                }

                if (has2ndChild) {

                    if (childConstPartialsIndex1 == null) {
                        if (childConstPartialsIndex2 == null) {
                            parentConstPartialsIndex[pIndex] = cst1 + cst2;
                        } else {
                            parentConstPartialsIndex[pIndex] = cst1 + MathFunctions.logSumExp(cstArr2);
                        }
                    } else {
                        if (childConstPartialsIndex2 == null) {
                            parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cstArr1) + cst2;
                        } else {
                            parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cstArr1) + MathFunctions.logSumExp(cstArr2);
                        }
                    }

                    parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);
                } else {

                    if (childConstPartialsIndex1 == null) {
                        parentConstPartialsIndex[pIndex] = cst1;
                    } else {
                        parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cstArr1);
                    }

                    parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1);
                }

                pIndex++;
            }
        }
    } // calculatePartialPartialPruning


    //****************************************************************
    //*      sum-product and max-sum algorithm (normal partial)      *
    //****************************************************************

    /**
     * calculate partial likelihoods at a node
     *
     * @param childIndex1       index of the first child
     * @param isLeaf1           is the first child is a leaf node?
     * @param childIndex2       index of the second child
     * @param isLeaf2           is the second child is a leaf node?
     * @param parentIndex       index of parent node
     * @param parentMLGenotypes record maximum likelihood genotype path at a parent node for its children
     *                          [#matrices * #patterns * #states] * 2
     * @param constGenotype genotype of constant site
     */
    public void calculatePartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype
    ) {
        if (isLeaf1) {
            // childIndex1 is a leaf node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLeafLeafPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype
                );
            } else {
                // childIndex2 is an internal node
                calculateLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype
                );
            }
        } else {
            // childIndex1 is an internal node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        true,
                        parentMLGenotypes,
                        constGenotype
                );
            } else {
                // childIndex2 is an internal node
                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes
                );
            }
        }

        if (useScaling) {
            scalePartials(parentIndex);
        }
    } // calculatePartials

    /**
     * calculate partial likelihoods at a node which has only one child
     * normally used for tree root and its only direct child
     * note that the child could only be internal
     *
     * @param childIndex        index of the child
     * @param parentIndex       index of the parent
     * @param parentMLGenotypes record maximum likelihood genotype path at a parent node for its children
     *                          [#matrices * #patterns * #states] * 2
     */
    public void calculatePartials(
            final int childIndex,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes
    ) {
        calculatePartialPartialPruning(
                partials[currentPartialsIndex[childIndex]][childIndex],
                constPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                matrices[currentMatrixIndex[childIndex]][childIndex],
                null,
                null,
                null,
                null,
                partials[currentPartialsIndex[parentIndex]][parentIndex],
                constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                false,
                parentMLGenotypes
        );

        if (useScaling) {
            scalePartials(parentIndex);
        }
    } // calculatePartials

    /**
     * Calculates partial likelihoods at a node when both children are leaf nodes.
     *
     * @param leafPartialsIndex1       #matrices * #patterns * #states
     * @param matricesIndex1           #matrices * #states * #states
     * @param leafPartialsIndex2       #matrices * #patterns * #states
     * @param matricesIndex2           #matrices * #states * #states
     * @param parentPartialsIndex      #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param parentMLPartialsIndex    #matrices * #patterns * #states
     * @param reversedChildrenOrder    whether the order of the children is reversed
     * @param parentMLGenotypes        record maximum likelihood genotype path at a parent node for its children
     *                                 [#matrices * #patterns * #states] * 2
     * @param constGenotype genotype of constant site
     */
    protected void calculateLeafLeafPruning(
            final double[] leafPartialsIndex1,
            final double[] matricesIndex1,
            final double[] leafPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // sum-product
        double sum1, sum2;

        // max-sum
        double max1, max2;

        // constant site
        double cst1, cst2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    // sum-product
                    sum1 = sum2 = 0.0;

                    // max-sum
                    max1 = max2 = 0.0;
                    if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                        parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild1Index].clear();
                    }
                    if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                        parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild2Index].clear();
                    }

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        tmp1 = matricesIndex1[mIndex] * leafPartialsIndex1[cIndex + cGenotypeIndex];
                        tmp2 = matricesIndex2[mIndex] * leafPartialsIndex2[cIndex + cGenotypeIndex];

                        // constant site
                        if (cGenotypeIndex == constGenotype) {
                            if (leafPartialsIndex1[cIndex + cGenotypeIndex] <= 0.0) {
                                cst1 = matricesIndex1[mIndex] * Double.MIN_VALUE;
                            } else {
                                cst1 = tmp1;
                            }

                            if (leafPartialsIndex2[cIndex + cGenotypeIndex] <= 0.0) {
                                cst2 = matricesIndex2[mIndex] * Double.MIN_VALUE;
                            } else {
                                cst2 = tmp2;
                            }

                            parentConstPartialsIndex[pIndex] = cst1 * cst2;
                        }

                        // sum-product
                        sum1 += tmp1;
                        sum2 += tmp2;

                        // max-sum
                        if (cGenotypeIndex == 0) {
                            max1 = Math.log(tmp1);
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            max2 = Math.log(tmp2);
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else {
                            if (max1 < Math.log(tmp1)) {
                                max1 = Math.log(tmp1);
                                parentMLGenotypes[pIndex][MLChild1Index].clear();
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            } else if (max1 == Math.log(tmp1)) {
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            }

                            if (max2 < Math.log(tmp2)) {
                                max2 = Math.log(tmp2);
                                parentMLGenotypes[pIndex][MLChild2Index].clear();
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            } else if (max2 == Math.log(tmp2)) {
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        }

                        mIndex++;
                    }

                    // sum-product
                    parentPartialsIndex[pIndex] = sum1 * sum2;

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1 + max2;

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculateLeafLeafPruning

    /**
     * Calculates partial likelihoods at a node when one child is a leaf node and the other is an internal node.
     *
     * @param leafPartialsIndex          #matrices * #patterns * #states
     * @param matricesIndex1             #matrices * #states * #states
     * @param internalPartialsIndex      #matrices * #patterns * #states
     * @param internalConstPartialsIndex #matrices * #patterns * #states
     * @param internalMLPartialsIndex    #matrices * #patterns * #states
     * @param matricesIndex2             #matrices * #states * #states
     * @param parentPartialsIndex        #matrices * #patterns * #states
     * @param parentConstPartialsIndex   #matrices * #patterns * #states
     * @param parentMLPartialsIndex      #matrices * #patterns * #states
     * @param reversedChildrenOrder      whether the order of the children is reversed
     * @param parentMLGenotypes          record maximum likelihood genotype path at a parent node for its children
     *                                   [#matrices * #patterns * #states] * 2
     * @param constGenotype              genotype of constant site
     */
    protected void calculateLeafPartialPruning(
            final double[] leafPartialsIndex,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex,
            final double[] internalConstPartialsIndex,
            final double[] internalMLPartialsIndex,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // constant site
        double cst1, cst2;

        // sum-product
        // sum1: for leaf node child
        // sum2: for internal node child
        double sum1, sum2;

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    // constant site
                    cst1 = cst2 = 0.0;

                    // sum-product
                    sum1 = sum2 = 0.0;

                    // max-sum
                    max1 = max2 = 0.0;
                    if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                        parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild1Index].clear();
                    }
                    if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                        parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild2Index].clear();
                    }

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        // for leaf node
                        tmp1 = matricesIndex1[mIndex] * leafPartialsIndex[cIndex + cGenotypeIndex];

                        // constant site
                        if (cGenotypeIndex == constGenotype) {
                            if (leafPartialsIndex[cIndex + cGenotypeIndex] <= 0.0) {
                                cst1 = matricesIndex1[mIndex] * Double.MIN_VALUE;
                            } else {
                                cst1 = tmp1;
                            }
                        }

                        // sum-product
                        sum1 += tmp1;

                        // max-sum
                        if (cGenotypeIndex == 0) {
                            max1 = Math.log(tmp1);
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else {
                            if (max1 < Math.log(tmp1)) {
                                max1 = Math.log(tmp1);
                                parentMLGenotypes[pIndex][MLChild1Index].clear();
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            } else if (max1 == Math.log(tmp1)) {
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            }
                        }


                        // for internal node

                        // constant site
                        cst2 += matricesIndex2[mIndex] * internalConstPartialsIndex[cIndex + cGenotypeIndex];

                        // sum-product
                        sum2 += matricesIndex2[mIndex] * internalPartialsIndex[cIndex + cGenotypeIndex];

                        // max-sum
                        tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex[cIndex + cGenotypeIndex];
                        if (cGenotypeIndex == 0) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else {
                            if (max2 < tmp2) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].clear();
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            } else if (max2 == tmp2) {
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        }

                        mIndex++;
                    }

                    // constant site
                    parentConstPartialsIndex[pIndex] = cst1 * cst2;

                    // sum-product
                    parentPartialsIndex[pIndex] = sum1 * sum2;

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1 + max2;

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculateLeafPartialPruning

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param internalPartialsIndex1      #matrices * #patterns * #states
     * @param internalConstPartialsIndex1 #matrices * #patterns * #states
     * @param internalMLPartialsIndex1    #matrices * #patterns * #states
     * @param matricesIndex1              #matrices  * #states * #states
     * @param internalPartialsIndex2      #matrices * #patterns * #states, allowed to be null
     * @param internalConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param internalMLPartialsIndex2    #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2              #matrices  * #states * #states, allowed to be null
     * @param parentPartialsIndex         #matrices * #patterns * #states
     * @param parentConstPartialsIndex    #matrices * #patterns * #states
     * @param parentMLPartialsIndex       #matrices * #patterns * #states
     * @param reversedChildrenOrder       whether the order of the children is reversed
     * @param parentMLGenotypes           record maximum likelihood genotype path at a parent node for its children
     *                                    [#matrices * #patterns * #states] * 2
     */
    protected void calculatePartialPartialPruning(
            final double[] internalPartialsIndex1,
            final double[] internalConstPartialsIndex1,
            final double[] internalMLPartialsIndex1,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex2,
            final double[] internalConstPartialsIndex2,
            final double[] internalMLPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        if (!(internalPartialsIndex2 == null && internalMLPartialsIndex2 == null && matricesIndex2 == null) &&
                !(internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null)) {
            throw new IllegalArgumentException("internalPartialsIndex2, internalMLPartialsIndex2 and matricesIndex2 " +
                    "should be defined or be null synchronously (" + this.getClass().getName() + ")");
        }

        if (internalPartialsIndex1 == null || internalMLPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("internalPartialsIndex1, internalMLPartialsIndex1 and " +
                    "matricesIndex1 should be defined instead of being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null);

        // constant site
        double cst1, cst2;

        // sum-product
        double sum1, sum2;

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    // constant site
                    cst1 = cst2 = 0.0;

                    // sum-product
                    sum1 = sum2 = 0.0;

                    // max-sum
                    max1 = max2 = 0.0;
                    if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                        parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild1Index].clear();
                    }
                    if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                        parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild2Index].clear();
                    }

                    // intermediate values
                    tmp2 = 0.0;

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        // constant site
                        cst1 += matricesIndex1[mIndex] * internalConstPartialsIndex1[cIndex + cGenotypeIndex];
                        if (has2ndChild) {
                            cst2 += matricesIndex2[mIndex] * internalConstPartialsIndex2[cIndex + cGenotypeIndex];
                        }

                        // sum-product
                        sum1 += matricesIndex1[mIndex] * internalPartialsIndex1[cIndex + cGenotypeIndex];
                        if (has2ndChild) {
                            sum2 += matricesIndex2[mIndex] * internalPartialsIndex2[cIndex + cGenotypeIndex];
                        }

                        // max-sum
                        tmp1 = Math.log(matricesIndex1[mIndex]) + internalMLPartialsIndex1[cIndex + cGenotypeIndex];
                        if (has2ndChild) {
                            tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex2[cIndex + cGenotypeIndex];
                        }
                        if (cGenotypeIndex == 0) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            if (has2ndChild) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        } else {
                            if (max1 < tmp1) {
                                max1 = tmp1;
                                parentMLGenotypes[pIndex][MLChild1Index].clear();
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            } else if (max1 == tmp1) {
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            }

                            if (has2ndChild) {
                                if (max2 < tmp2) {
                                    max2 = tmp2;
                                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                                    parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                                } else if (max2 == tmp2) {
                                    parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                                }
                            }
                        }

                        mIndex++;
                    }

                    if (has2ndChild) {
                        // constant site
                        parentConstPartialsIndex[pIndex] = cst1 * cst2;

                        // sum-product
                        parentPartialsIndex[pIndex] = sum1 * sum2;

                        // max-sum
                        parentMLPartialsIndex[pIndex] = max1 + max2;
                    } else {
                        // constant site
                        parentConstPartialsIndex[pIndex] = cst1;

                        // sum-product
                        parentPartialsIndex[pIndex] = sum1;

                        // max-sum
                        parentMLPartialsIndex[pIndex] = max1;
                    }

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculatePartialPartialPruning

    /**
     * calculate partial likelihoods at a node
     *
     * @param childIndex1          index of the first child
     * @param isLeaf1              is the first child is a leaf node?
     * @param childIndex2          index of the second child
     * @param isLeaf2              is the second child is a leaf node?
     * @param parentIndex          index of parent node
     * @param parentMLGenotypes    record maximum likelihood genotype path at a parent node for its children
     *                             [#matrices * #patterns * #states] * 2
     * @param constGenotype genotype of constant site
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> calculatePartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype,
            final int[][] changedPatterns,
            final List<Integer> changedPatternsIndex
    ) {
        if (isLeaf1) {
            // childIndex1 is a leaf node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLeafLeafPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype,
                        changedPatterns
                );
            } else {
                // childIndex2 is an internal node
                calculateLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype,
                        changedPatterns
                );
            }
        } else {
            // childIndex1 is an internal node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        true,
                        parentMLGenotypes,
                        constGenotype,
                        changedPatterns
                );
            } else {
                // childIndex2 is an internal node
                calculatePartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        changedPatterns
                );
            }
        }

        if (useScaling) {
            return scalePartials(parentIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // calculatePartials

    /**
     * calculate partial likelihoods at a node which has only one child
     * normally used for tree root and its only direct child
     * note that the child could only be internal
     *
     * @param childIndex           index of the child
     * @param parentIndex          index of the parent
     * @param parentMLGenotypes    record maximum likelihood genotype path at a parent node for its children
     *                             [#matrices * #patterns * #states] * 2
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> calculatePartials(
            final int childIndex,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes,
            final int[][] changedPatterns,
            final List<Integer> changedPatternsIndex
    ) {
        calculatePartialPartialPruning(
                partials[currentPartialsIndex[childIndex]][childIndex],
                constPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                matrices[currentMatrixIndex[childIndex]][childIndex],
                null,
                null,
                null,
                null,
                partials[currentPartialsIndex[parentIndex]][parentIndex],
                constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                false,
                parentMLGenotypes,
                changedPatterns
        );

        if (useScaling) {
            return scalePartials(parentIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // calculatePartials

    /**
     * Calculates partial likelihoods at a node when both children are leaf nodes.
     *
     * @param leafPartialsIndex1       #matrices * #patterns * #states
     * @param matricesIndex1           #matrices * #states * #states
     * @param leafPartialsIndex2       #matrices * #patterns * #states
     * @param matricesIndex2           #matrices * #states * #states
     * @param parentPartialsIndex      #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param parentMLPartialsIndex    #matrices * #patterns * #states
     * @param reversedChildrenOrder    whether the order of the children is reversed
     * @param parentMLGenotypes        record maximum likelihood genotype path at a parent node for its children
     *                                 [#matrices * #patterns * #states] * 2
     * @param constGenotype            genotype of constant site
     * @param changedPatterns          patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculateLeafLeafPruning(
            final double[] leafPartialsIndex1,
            final double[] matricesIndex1,
            final double[] leafPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype,
            final int[][] changedPatterns
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // sum-product
        double sum1, sum2;

        // max-sum
        double max1, max2;

        // constant site
        double cst1, cst2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                // sum-product
                sum1 = sum2 = 0.0;

                // max-sum
                max1 = max2 = 0.0;
                if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                    parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild1Index].clear();
                }
                if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                    parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                }

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    tmp1 = matricesIndex1[mIndex] * leafPartialsIndex1[cIndex + cGenotypeIndex];
                    tmp2 = matricesIndex2[mIndex] * leafPartialsIndex2[cIndex + cGenotypeIndex];

                    // constant site
                    if (cGenotypeIndex == constGenotype) {
                        if (leafPartialsIndex1[cIndex + cGenotypeIndex] <= 0.0) {
                            cst1 = matricesIndex1[mIndex] * Double.MIN_VALUE;
                        } else {
                            cst1 = tmp1;
                        }

                        if (leafPartialsIndex2[cIndex + cGenotypeIndex] <= 0.0) {
                            cst2 = matricesIndex2[mIndex] * Double.MIN_VALUE;
                        } else {
                            cst2 = tmp2;
                        }

                        parentConstPartialsIndex[pIndex] = cst1 * cst2;
                    }

                    // sum-product
                    sum1 += tmp1;
                    sum2 += tmp2;

                    // max-sum
                    if (cGenotypeIndex == 0) {
                        max1 = Math.log(tmp1);
                        parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        max2 = Math.log(tmp2);
                        parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                    } else {
                        if (max1 < Math.log(tmp1)) {
                            max1 = Math.log(tmp1);
                            parentMLGenotypes[pIndex][MLChild1Index].clear();
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else if (max1 == Math.log(tmp1)) {
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        }

                        if (max2 < Math.log(tmp2)) {
                            max2 = Math.log(tmp2);
                            parentMLGenotypes[pIndex][MLChild2Index].clear();
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else if (max2 == Math.log(tmp2)) {
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        }
                    }

                    mIndex++;
                }

                // sum-product
                parentPartialsIndex[pIndex] = sum1 * sum2;

                // max-sum
                parentMLPartialsIndex[pIndex] = max1 + max2;

                pIndex++;
            }
        }
    } // calculateLeafLeafPruning

    /**
     * Calculates partial likelihoods at a node when one child is a leaf node and the other is an internal node.
     *
     * @param leafPartialsIndex       #matrices * #patterns * #states
     * @param matricesIndex1          #matrices * #states * #states
     * @param internalPartialsIndex   #matrices * #patterns * #states
     * @param internalMLPartialsIndex #matrices * #patterns * #states
     * @param matricesIndex2          #matrices * #states * #states
     * @param parentPartialsIndex     #matrices * #patterns * #states
     * @param parentMLPartialsIndex   #matrices * #patterns * #states
     * @param reversedChildrenOrder   whether the order of the children is reversed
     * @param parentMLGenotypes       record maximum likelihood genotype path at a parent node for its children
     *                                [#matrices * #patterns * #states] * 2
     * @param constGenotype           genotype of constant site
     * @param changedPatterns         patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculateLeafPartialPruning(
            final double[] leafPartialsIndex,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex,
            final double[] internalConstPartialsIndex,
            final double[] internalMLPartialsIndex,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype,
            final int[][] changedPatterns
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // constant site
        double cst1, cst2;

        // sum-product
        // sum1: for leaf node child
        // sum2: for internal node child
        double sum1, sum2;

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                // constant site
                cst1 = cst2 = 0.0;

                // sum-product
                sum1 = sum2 = 0.0;

                // max-sum
                max1 = max2 = 0.0;
                if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                    parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild1Index].clear();
                }
                if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                    parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                }

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    // for leaf node
                    tmp1 = matricesIndex1[mIndex] * leafPartialsIndex[cIndex + cGenotypeIndex];

                    // constant site
                    if (cGenotypeIndex == constGenotype) {
                        if (leafPartialsIndex[cIndex + cGenotypeIndex] <= 0.0) {
                            cst1 = matricesIndex1[mIndex] * Double.MIN_VALUE;
                        } else {
                            cst1 = tmp1;
                        }
                    }

                    // sum-product
                    sum1 += tmp1;

                    // max-sum
                    if (cGenotypeIndex == 0) {
                        max1 = Math.log(tmp1);
                        parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                    } else {
                        if (max1 < Math.log(tmp1)) {
                            max1 = Math.log(tmp1);
                            parentMLGenotypes[pIndex][MLChild1Index].clear();
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else if (max1 == Math.log(tmp1)) {
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        }
                    }


                    // for internal node

                    // constant site
                    cst2 += matricesIndex2[mIndex] * internalConstPartialsIndex[cIndex + cGenotypeIndex];

                    // sum-product
                    sum2 += matricesIndex2[mIndex] * internalPartialsIndex[cIndex + cGenotypeIndex];

                    // max-sum
                    tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex[cIndex + cGenotypeIndex];
                    if (cGenotypeIndex == 0) {
                        max2 = tmp2;
                        parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                    } else {
                        if (max2 < tmp2) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].clear();
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else if (max2 == tmp2) {
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        }
                    }

                    mIndex++;
                }

                // constant site
                parentConstPartialsIndex[pIndex] = cst1 * cst2;

                // sum-product
                parentPartialsIndex[pIndex] = sum1 * sum2;

                // max-sum
                parentMLPartialsIndex[pIndex] = max1 + max2;

                pIndex++;
            }
        }
    } // calculateLeafPartialPruning

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param internalPartialsIndex1      #matrices * #patterns * #states
     * @param internalConstPartialsIndex1 #matrices * #patterns * #states
     * @param internalMLPartialsIndex1    #matrices * #patterns * #states
     * @param matricesIndex1              #matrices  * #states * #states
     * @param internalPartialsIndex2      #matrices * #patterns * #states, allowed to be null
     * @param internalConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param internalMLPartialsIndex2    #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2              #matrices  * #states * #states, allowed to be null
     * @param parentPartialsIndex         #matrices * #patterns * #states
     * @param parentConstPartialsIndex    #matrices * #patterns * #states
     * @param parentMLPartialsIndex       #matrices * #patterns * #states
     * @param reversedChildrenOrder       whether the order of the children is reversed
     * @param parentMLGenotypes           record maximum likelihood genotype path at a parent node for its children
     *                                    [#matrices * #patterns * #states] * 2
     * @param changedPatterns             patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculatePartialPartialPruning(
            final double[] internalPartialsIndex1,
            final double[] internalConstPartialsIndex1,
            final double[] internalMLPartialsIndex1,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex2,
            final double[] internalConstPartialsIndex2,
            final double[] internalMLPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int[][] changedPatterns
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        if (!(internalPartialsIndex2 == null && internalMLPartialsIndex2 == null && matricesIndex2 == null) &&
                !(internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null)) {
            throw new IllegalArgumentException("internalPartialsIndex2, internalMLPartialsIndex2 and matricesIndex2 " +
                    "should be defined or be null synchronously (" + this.getClass().getName() + ")");
        }

        if (internalPartialsIndex1 == null || internalMLPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("internalPartialsIndex1, internalMLPartialsIndex1 and " +
                    "matricesIndex1 should be defined instead of being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null);

        // constant site
        double cst1, cst2;

        // sum-product
        double sum1, sum2;

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                // constant site
                cst1 = cst2 = 0.0;

                // sum-product
                sum1 = sum2 = 0.0;

                // max-sum
                max1 = max2 = 0.0;
                if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                    parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild1Index].clear();
                }
                if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                    parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                }

                // intermediate values
                tmp2 = 0.0;

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    // constant site
                    cst1 += matricesIndex1[mIndex] * internalConstPartialsIndex1[cIndex + cGenotypeIndex];
                    if (has2ndChild) {
                        cst2 += matricesIndex2[mIndex] * internalConstPartialsIndex2[cIndex + cGenotypeIndex];
                    }

                    // sum-product
                    sum1 += matricesIndex1[mIndex] * internalPartialsIndex1[cIndex + cGenotypeIndex];
                    if (has2ndChild) {
                        sum2 += matricesIndex2[mIndex] * internalPartialsIndex2[cIndex + cGenotypeIndex];
                    }

                    // max-sum
                    tmp1 = Math.log(matricesIndex1[mIndex]) + internalMLPartialsIndex1[cIndex + cGenotypeIndex];
                    if (has2ndChild) {
                        tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex2[cIndex + cGenotypeIndex];
                    }
                    if (cGenotypeIndex == 0) {
                        max1 = tmp1;
                        parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        if (has2ndChild) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        }
                    } else {
                        if (max1 < tmp1) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].clear();
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else if (max1 == tmp1) {
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        }

                        if (has2ndChild) {
                            if (max2 < tmp2) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].clear();
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            } else if (max2 == tmp2) {
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        }
                    }

                    mIndex++;
                }

                if (has2ndChild) {
                    // constant site
                    parentConstPartialsIndex[pIndex] = cst1 * cst2;

                    // sum-product
                    parentPartialsIndex[pIndex] = sum1 * sum2;

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1 + max2;
                } else {
                    // constant site
                    parentConstPartialsIndex[pIndex] = cst1;

                    // sum-product
                    parentPartialsIndex[pIndex] = sum1;

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1;
                }

                pIndex++;
            }
        }
    } // calculatePartialPartialPruning


    //*************************************************************
    //*      sum-product and max-sum algorithm (log partial)      *
    //*************************************************************

    /**
     * calculate partial likelihoods at a node
     *
     * @param childIndex1       index of the first child
     * @param isLeaf1           is the first child is a leaf node?
     * @param childIndex2       index of the second child
     * @param isLeaf2           is the second child is a leaf node?
     * @param parentIndex       index of parent node
     * @param parentMLGenotypes record maximum likelihood genotype path at a parent node for its children
     *                          [#matrices * #patterns * #states] * 2
     * @param constGenotype     genotype of constant site
     */
    public void calculateLogPartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype
    ) {
        if (isLeaf1) {
            // childIndex1 is a leaf node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLogLeafLeafPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype
                );
            } else {
                // childIndex2 is an internal node
                calculateLogLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype
                );
            }
        } else {
            // childIndex1 is an internal node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLogLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        true,
                        parentMLGenotypes,
                        constGenotype
                );
            } else {
                // childIndex2 is an internal node
                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes
                );
            }
        }

        if (useScaling) {
            scalePartials(parentIndex);
        }
    } // calculateLogPartials

    /**
     * calculate partial likelihoods at a node which has only one child
     * normally used for tree root and its only direct child
     * note that the child could only be internal
     *
     * @param childIndex        index of the child
     * @param parentIndex       index of the parent
     * @param parentMLGenotypes record maximum likelihood genotype path at a parent node for its children
     *                          [#matrices * #patterns * #states] * 2
     */
    public void calculateLogPartials(
            final int childIndex,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes
    ) {
        calculateLogPartialPartialPruning(
                partials[currentPartialsIndex[childIndex]][childIndex],
                constPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                matrices[currentMatrixIndex[childIndex]][childIndex],
                null,
                null,
                null,
                null,
                partials[currentPartialsIndex[parentIndex]][parentIndex],
                constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                false,
                parentMLGenotypes
        );

        if (useScaling) {
            scalePartials(parentIndex);
        }
    } // calculateLogPartials

    /**
     * Calculates partial likelihoods at a node when both children are leaf nodes.
     *
     * @param leafPartialsIndex1       #matrices * #patterns * #states
     * @param matricesIndex1           #matrices * #states * #states
     * @param leafPartialsIndex2       #matrices * #patterns * #states
     * @param matricesIndex2           #matrices * #states * #states
     * @param parentPartialsIndex      #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param parentMLPartialsIndex    #matrices * #patterns * #states
     * @param reversedChildrenOrder    whether the order of the children is reversed
     * @param parentMLGenotypes        record maximum likelihood genotype path at a parent node for its children
     *                                 [#matrices * #patterns * #states] * 2
     * @param constGenotype            genotype of constant site
     */
    protected void calculateLogLeafLeafPruning(
            final double[] leafPartialsIndex1,
            final double[] matricesIndex1,
            final double[] leafPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // sum-product
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    // max-sum
                    max1 = max2 = 0.0;
                    if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                        parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild1Index].clear();
                    }
                    if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                        parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild2Index].clear();
                    }

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        tmp1 = Math.log(matricesIndex1[mIndex]) + leafPartialsIndex1[cIndex + cGenotypeIndex];
                        tmp2 = Math.log(matricesIndex2[mIndex]) + leafPartialsIndex2[cIndex + cGenotypeIndex];

                        // constant site
                        if (cGenotypeIndex == constGenotype) {
                            parentConstPartialsIndex[pIndex] = tmp1 + tmp2;
                        }

                        // sum-product
                        sp1[cGenotypeIndex] = tmp1;
                        sp2[cGenotypeIndex] = tmp2;

                        // max-sum
                        if (cGenotypeIndex == 0) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else {
                            if (max1 < tmp1) {
                                max1 = tmp1;
                                parentMLGenotypes[pIndex][MLChild1Index].clear();
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            } else if (max1 == tmp1) {
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            }

                            if (max2 < tmp2) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].clear();
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            } else if (max2 == tmp2) {
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        }

                        mIndex++;
                    }

                    // sum-product
                    parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1 + max2;

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculateLogLeafLeafPruning

    /**
     * Calculates partial likelihoods at a node when one child is a leaf node and the other is an internal node.
     *
     * @param leafPartialsIndex          #matrices * #patterns * #states
     * @param matricesIndex1             #matrices * #states * #states
     * @param internalPartialsIndex      #matrices * #patterns * #states
     * @param internalConstPartialsIndex #matrices * #patterns * #states
     * @param internalMLPartialsIndex    #matrices * #patterns * #states
     * @param matricesIndex2             #matrices * #states * #states
     * @param parentPartialsIndex        #matrices * #patterns * #states
     * @param parentConstPartialsIndex   #matrices * #patterns * #states
     * @param parentMLPartialsIndex      #matrices * #patterns * #states
     * @param reversedChildrenOrder      whether the order of the children is reversed
     * @param parentMLGenotypes          record maximum likelihood genotype path at a parent node for its children
     *                                   [#matrices * #patterns * #states] * 2
     * @param constGenotype              genotype of constant site
     */
    protected void calculateLogLeafPartialPruning(
            final double[] leafPartialsIndex,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex,
            final double[] internalConstPartialsIndex,
            final double[] internalMLPartialsIndex,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // constant site
        double cst1;
        double[] cst2 = new double[nrOfStates];

        // sum-product
        // sp1: for leaf node child
        // sp2: for internal node child
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    // constant site
                    cst1 = 0.0;

                    // max-sum
                    max1 = max2 = 0.0;
                    if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                        parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild1Index].clear();
                    }
                    if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                        parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild2Index].clear();
                    }

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        // for leaf node
                        tmp1 = Math.log(matricesIndex1[mIndex]) + leafPartialsIndex[cIndex + cGenotypeIndex];

                        // constant site
                        if (cGenotypeIndex == constGenotype) {
                            cst1 = tmp1;
                        }

                        // sum-product
                        sp1[cGenotypeIndex] = tmp1;

                        // max-sum
                        if (cGenotypeIndex == 0) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else {
                            if (max1 < tmp1) {
                                max1 = tmp1;
                                parentMLGenotypes[pIndex][MLChild1Index].clear();
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            } else if (max1 == tmp1) {
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            }
                        }


                        // for internal node

                        // constant site
                        cst2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalConstPartialsIndex[cIndex + cGenotypeIndex];

                        // sum-product
                        sp2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalPartialsIndex[cIndex + cGenotypeIndex];

                        // max-sum
                        tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex[cIndex + cGenotypeIndex];
                        if (cGenotypeIndex == 0) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else {
                            if (max2 < tmp2) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].clear();
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            } else if (max2 == tmp2) {
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        }

                        mIndex++;
                    }

                    // constant site
                    parentConstPartialsIndex[pIndex] = cst1 + MathFunctions.logSumExp(cst2);

                    // sum-product
                    parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1 + max2;

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculateLogLeafPartialPruning

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param internalPartialsIndex1      #matrices * #patterns * #states
     * @param internalConstPartialsIndex1 #matrices * #patterns * #states
     * @param internalMLPartialsIndex1    #matrices * #patterns * #states
     * @param matricesIndex1              #matrices  * #states * #states
     * @param internalPartialsIndex2      #matrices * #patterns * #states, allowed to be null
     * @param internalConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param internalMLPartialsIndex2    #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2              #matrices  * #states * #states, allowed to be null
     * @param parentPartialsIndex         #matrices * #patterns * #states
     * @param parentConstPartialsIndex    #matrices * #patterns * #states
     * @param parentMLPartialsIndex       #matrices * #patterns * #states
     * @param reversedChildrenOrder       whether the order of the children is reversed
     * @param parentMLGenotypes           record maximum likelihood genotype path at a parent node for its children
     *                                    [#matrices * #patterns * #states] * 2
     */
    protected void calculateLogPartialPartialPruning(
            final double[] internalPartialsIndex1,
            final double[] internalConstPartialsIndex1,
            final double[] internalMLPartialsIndex1,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex2,
            final double[] internalConstPartialsIndex2,
            final double[] internalMLPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        if (!(internalPartialsIndex2 == null && internalMLPartialsIndex2 == null && matricesIndex2 == null) &&
                !(internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null)) {
            throw new IllegalArgumentException("internalPartialsIndex2, internalMLPartialsIndex2 and matricesIndex2 " +
                    "should be defined or be null synchronously (" + this.getClass().getName() + ")");
        }

        if (internalPartialsIndex1 == null || internalMLPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("internalPartialsIndex1, internalMLPartialsIndex1 and " +
                    "matricesIndex1 should be defined instead of being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null);

        // constant site
        double[] cst1 = new double[nrOfStates];
        double[] cst2 = new double[nrOfStates];

        // sum-product
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;
        pIndex = cIndex = 0;

        for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

            for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                mIndex = matrixIndex * matrixSize;

                for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                    // max-sum
                    max1 = max2 = 0.0;
                    if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                        parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild1Index].clear();
                    }
                    if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                        parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                    } else {
                        parentMLGenotypes[pIndex][MLChild2Index].clear();
                    }

                    // intermediate values
                    tmp2 = 0.0;

                    for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                        // constant site
                        cst1[cGenotypeIndex] = Math.log(matricesIndex1[mIndex]) + internalConstPartialsIndex1[cIndex + cGenotypeIndex];
                        if (has2ndChild) {
                            cst2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalConstPartialsIndex2[cIndex + cGenotypeIndex];
                        }

                        // sum-product
                        sp1[cGenotypeIndex] = Math.log(matricesIndex1[mIndex]) + internalPartialsIndex1[cIndex + cGenotypeIndex];
                        if (has2ndChild) {
                            sp2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalPartialsIndex2[cIndex + cGenotypeIndex];
                        }

                        // max-sum
                        tmp1 = Math.log(matricesIndex1[mIndex]) + internalMLPartialsIndex1[cIndex + cGenotypeIndex];
                        if (has2ndChild) {
                            tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex2[cIndex + cGenotypeIndex];
                        }
                        if (cGenotypeIndex == 0) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            if (has2ndChild) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        } else {
                            if (max1 < tmp1) {
                                max1 = tmp1;
                                parentMLGenotypes[pIndex][MLChild1Index].clear();
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            } else if (max1 == tmp1) {
                                parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                            }

                            if (has2ndChild) {
                                if (max2 < tmp2) {
                                    max2 = tmp2;
                                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                                    parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                                } else if (max2 == tmp2) {
                                    parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                                }
                            }
                        }

                        mIndex++;
                    }

                    if (has2ndChild) {
                        // constant site
                        parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cst1) + MathFunctions.logSumExp(cst2);

                        // sum-product
                        parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                        // max-sum
                        parentMLPartialsIndex[pIndex] = max1 + max2;
                    } else {
                        // constant site
                        parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cst1);

                        // sum-product
                        parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1);

                        // max-sum
                        parentMLPartialsIndex[pIndex] = max1;
                    }

                    pIndex++;
                }

                cIndex += nrOfStates;
            }
        }
    } // calculateLogPartialPartialPruning


    /**
     * calculate partial likelihoods at a node
     *
     * @param childIndex1          index of the first child
     * @param isLeaf1              is the first child is a leaf node?
     * @param childIndex2          index of the second child
     * @param isLeaf2              is the second child is a leaf node?
     * @param parentIndex          index of parent node
     * @param parentMLGenotypes    record maximum likelihood genotype path at a parent node for its children
     *                             [#matrices * #patterns * #states] * 2
     * @param constGenotype        genotype of constant site
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> calculateLogPartials(
            final int childIndex1,
            final boolean isLeaf1,
            final int childIndex2,
            final boolean isLeaf2,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype,
            final int[][] changedPatterns,
            final List<Integer> changedPatternsIndex
    ) {
        if (isLeaf1) {
            // childIndex1 is a leaf node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLogLeafLeafPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype,
                        changedPatterns
                );
            } else {
                // childIndex2 is an internal node
                calculateLogLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        constGenotype,
                        changedPatterns
                );
            }
        } else {
            // childIndex1 is an internal node
            if (isLeaf2) {
                // childIndex2 is a leaf node
                calculateLogLeafPartialPruning(
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        true,
                        parentMLGenotypes,
                        constGenotype,
                        changedPatterns
                );
            } else {
                // childIndex2 is an internal node
                calculateLogPartialPartialPruning(
                        partials[currentPartialsIndex[childIndex1]][childIndex1],
                        constPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex1]][childIndex1 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex1]][childIndex1],
                        partials[currentPartialsIndex[childIndex2]][childIndex2],
                        constPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[childIndex2]][childIndex2 - nrOfLeafNodes],
                        matrices[currentMatrixIndex[childIndex2]][childIndex2],
                        partials[currentPartialsIndex[parentIndex]][parentIndex],
                        constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                        false,
                        parentMLGenotypes,
                        changedPatterns
                );
            }
        }

        if (useScaling) {
            return scalePartials(parentIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // calculateLogPartials

    /**
     * calculate partial likelihoods at a node which has only one child
     * normally used for tree root and its only direct child
     * note that the child could only be internal
     *
     * @param childIndex           index of the child
     * @param parentIndex          index of the parent
     * @param parentMLGenotypes    record maximum likelihood genotype path at a parent node for its children
     *                             [#matrices * #patterns * #states] * 2
     * @param changedPatterns      patterns have been changed (matrixIndex, patternIndex)
     * @param changedPatternsIndex first index of each pattern in changedPatterns
     */
    public List<int[]> calculateLogPartials(
            final int childIndex,
            final int parentIndex,
            List<Integer>[][] parentMLGenotypes,
            final int[][] changedPatterns,
            final List<Integer> changedPatternsIndex
    ) {
        calculateLogPartialPartialPruning(
                partials[currentPartialsIndex[childIndex]][childIndex],
                constPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[childIndex]][childIndex - nrOfLeafNodes],
                matrices[currentMatrixIndex[childIndex]][childIndex],
                null,
                null,
                null,
                null,
                partials[currentPartialsIndex[parentIndex]][parentIndex],
                constPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                MLPartials[currentPartialsIndex[parentIndex]][parentIndex - nrOfLeafNodes],
                false,
                parentMLGenotypes,
                changedPatterns
        );

        if (useScaling) {
            return scalePartials(parentIndex, changedPatterns, changedPatternsIndex);
        }

        return null;
    } // calculateLogPartials

    /**
     * Calculates partial likelihoods at a node when both children are leaf nodes.
     *
     * @param leafPartialsIndex1       #matrices * #patterns * #states
     * @param matricesIndex1           #matrices * #states * #states
     * @param leafPartialsIndex2       #matrices * #patterns * #states
     * @param matricesIndex2           #matrices * #states * #states
     * @param parentPartialsIndex      #matrices * #patterns * #states
     * @param parentConstPartialsIndex #matrices * #patterns * #states
     * @param parentMLPartialsIndex    #matrices * #patterns * #states
     * @param reversedChildrenOrder    whether the order of the children is reversed
     * @param parentMLGenotypes        record maximum likelihood genotype path at a parent node for its children
     *                                 [#matrices * #patterns * #states] * 2
     * @param constGenotype            genotype of constant site
     * @param changedPatterns          patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculateLogLeafLeafPruning(
            final double[] leafPartialsIndex1,
            final double[] matricesIndex1,
            final double[] leafPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int constGenotype,
            final int[][] changedPatterns
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // sum-product
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                // max-sum
                max1 = max2 = 0.0;
                if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                    parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild1Index].clear();
                }
                if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                    parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                }

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    tmp1 = Math.log(matricesIndex1[mIndex]) + leafPartialsIndex1[cIndex + cGenotypeIndex];
                    tmp2 = Math.log(matricesIndex2[mIndex]) + leafPartialsIndex2[cIndex + cGenotypeIndex];

                    // constant site
                    if (cGenotypeIndex == constGenotype) {
                        parentConstPartialsIndex[pIndex] = tmp1 + tmp2;
                    }

                    // sum-product
                    sp1[cGenotypeIndex] = tmp1;
                    sp2[cGenotypeIndex] = tmp2;

                    // max-sum
                    if (cGenotypeIndex == 0) {
                        max1 = tmp1;
                        parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        max2 = tmp2;
                        parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                    } else {
                        if (max1 < tmp1) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].clear();
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else if (max1 == tmp1) {
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        }

                        if (max2 < tmp2) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].clear();
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else if (max2 == tmp2) {
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        }
                    }

                    mIndex++;
                }

                // sum-product
                parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                // max-sum
                parentMLPartialsIndex[pIndex] = max1 + max2;

                pIndex++;
            }
        }
    } // calculateLogLeafLeafPruning

    /**
     * Calculates partial likelihoods at a node when one child is a leaf node and the other is an internal node.
     *
     * @param leafPartialsIndex       #matrices * #patterns * #states
     * @param matricesIndex1          #matrices * #states * #states
     * @param internalPartialsIndex   #matrices * #patterns * #states
     * @param internalMLPartialsIndex #matrices * #patterns * #states
     * @param matricesIndex2          #matrices * #states * #states
     * @param parentPartialsIndex     #matrices * #patterns * #states
     * @param parentMLPartialsIndex   #matrices * #patterns * #states
     * @param reversedChildrenOrder   whether the order of the children is reversed
     * @param parentMLGenotypes       record maximum likelihood genotype path at a parent node for its children
     *                                [#matrices * #patterns * #states] * 2
     * @param constGenotype           genotype of constant site
     * @param changedPatterns         patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculateLogLeafPartialPruning(
            final double[] leafPartialsIndex,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex,
            final double[] internalConstPartialsIndex,
            final double[] internalMLPartialsIndex,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            int constGenotype,
            final int[][] changedPatterns
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        // constant site
        double cst1;
        double[] cst2 = new double[nrOfStates];

        // sum-product
        // sp1: for leaf node child
        // sp2: for internal node child
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                // constant site
                cst1 = 0.0;

                // max-sum
                max1 = max2 = 0.0;
                if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                    parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild1Index].clear();
                }
                if (parentMLGenotypes[pIndex][MLChild2Index] == null) {
                    parentMLGenotypes[pIndex][MLChild2Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild2Index].clear();
                }

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    // for leaf node
                    tmp1 = Math.log(matricesIndex1[mIndex]) + leafPartialsIndex[cIndex + cGenotypeIndex];

                    // constant site
                    if (cGenotypeIndex == constGenotype) {
                        cst1 = tmp1;
                    }

                    // sum-product
                    sp1[cGenotypeIndex] = tmp1;

                    // max-sum
                    if (cGenotypeIndex == 0) {
                        max1 = tmp1;
                        parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                    } else {
                        if (max1 < tmp1) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].clear();
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else if (max1 == tmp1) {
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        }
                    }


                    // for internal node

                    // constant site
                    cst2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalConstPartialsIndex[cIndex + cGenotypeIndex];

                    // sum-product
                    sp2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalPartialsIndex[cIndex + cGenotypeIndex];

                    // max-sum
                    tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex[cIndex + cGenotypeIndex];
                    if (cGenotypeIndex == 0) {
                        max2 = tmp2;
                        parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                    } else {
                        if (max2 < tmp2) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].clear();
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        } else if (max2 == tmp2) {
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        }
                    }

                    mIndex++;
                }

                // constant site
                parentConstPartialsIndex[pIndex] = cst1 + MathFunctions.logSumExp(cst2);

                // sum-product
                parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                // max-sum
                parentMLPartialsIndex[pIndex] = max1 + max2;

                pIndex++;
            }
        }
    } // calculateLogLeafPartialPruning

    /**
     * Calculates partial likelihoods at a node when both children are internal nodes.
     *
     * @param internalPartialsIndex1      #matrices * #patterns * #states
     * @param internalConstPartialsIndex1 #matrices * #patterns * #states
     * @param internalMLPartialsIndex1    #matrices * #patterns * #states
     * @param matricesIndex1              #matrices  * #states * #states
     * @param internalPartialsIndex2      #matrices * #patterns * #states, allowed to be null
     * @param internalConstPartialsIndex2 #matrices * #patterns * #states, allowed to be null
     * @param internalMLPartialsIndex2    #matrices * #patterns * #states, allowed to be null
     * @param matricesIndex2              #matrices  * #states * #states, allowed to be null
     * @param parentPartialsIndex         #matrices * #patterns * #states
     * @param parentConstPartialsIndex    #matrices * #patterns * #states
     * @param parentMLPartialsIndex       #matrices * #patterns * #states
     * @param reversedChildrenOrder       whether the order of the children is reversed
     * @param parentMLGenotypes           record maximum likelihood genotype path at a parent node for its children
     *                                    [#matrices * #patterns * #states] * 2
     * @param changedPatterns             patterns have been changed (matrixIndex, patternIndex)
     */
    protected void calculateLogPartialPartialPruning(
            final double[] internalPartialsIndex1,
            final double[] internalConstPartialsIndex1,
            final double[] internalMLPartialsIndex1,
            final double[] matricesIndex1,
            final double[] internalPartialsIndex2,
            final double[] internalConstPartialsIndex2,
            final double[] internalMLPartialsIndex2,
            final double[] matricesIndex2,
            double[] parentPartialsIndex,
            final double[] parentConstPartialsIndex,
            double[] parentMLPartialsIndex,
            boolean reversedChildrenOrder,
            List<Integer>[][] parentMLGenotypes,
            final int[][] changedPatterns
    ) {
        final int MLChild1Index = reversedChildrenOrder ? 1 : 0;
        final int MLChild2Index = 1 - MLChild1Index;

        if (!(internalPartialsIndex2 == null && internalMLPartialsIndex2 == null && matricesIndex2 == null) &&
                !(internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null)) {
            throw new IllegalArgumentException("internalPartialsIndex2, internalMLPartialsIndex2 and matricesIndex2 " +
                    "should be defined or be null synchronously (" + this.getClass().getName() + ")");
        }

        if (internalPartialsIndex1 == null || internalMLPartialsIndex1 == null || matricesIndex1 == null) {
            throw new IllegalArgumentException("internalPartialsIndex1, internalMLPartialsIndex1 and " +
                    "matricesIndex1 should be defined instead of being null (" + this.getClass().getName() + ")");
        }

        final boolean has2ndChild = (internalPartialsIndex2 != null && internalMLPartialsIndex2 != null && matricesIndex2 != null);

        // constant site
        double[] cst1 = new double[nrOfStates];
        double[] cst2 = new double[nrOfStates];

        // sum-product
        double[] sp1 = new double[nrOfStates];
        double[] sp2 = new double[nrOfStates];

        // max-sum
        double max1, max2;

        // intermediate values
        double tmp1, tmp2;

        int pIndex, cIndex, mIndex;

        for (int[] pair : changedPatterns) {

            mIndex = pair[0] * matrixSize;
            pIndex = cIndex = pair[0] * nrOfPatterns * nrOfStates + pair[1] * nrOfStates;

            for (int pGenotypeIndex = 0; pGenotypeIndex < nrOfStates; pGenotypeIndex++) {

                // max-sum
                max1 = max2 = 0.0;
                if (parentMLGenotypes[pIndex][MLChild1Index] == null) {
                    parentMLGenotypes[pIndex][MLChild1Index] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][MLChild1Index].clear();
                }
                if (parentMLGenotypes[pIndex][1] == null) {
                    parentMLGenotypes[pIndex][1] = new ArrayList<>();
                } else {
                    parentMLGenotypes[pIndex][1].clear();
                }

                // intermediate values
                tmp2 = 0.0;

                for (int cGenotypeIndex = 0; cGenotypeIndex < nrOfStates; cGenotypeIndex++) {

                    // constant site
                    cst1[cGenotypeIndex] = Math.log(matricesIndex1[mIndex]) + internalConstPartialsIndex1[cIndex + cGenotypeIndex];
                    if (has2ndChild) {
                        cst2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalConstPartialsIndex2[cIndex + cGenotypeIndex];
                    }

                    // sum-product
                    sp1[cGenotypeIndex] = Math.log(matricesIndex1[mIndex]) + internalPartialsIndex1[cIndex + cGenotypeIndex];
                    if (has2ndChild) {
                        sp2[cGenotypeIndex] = Math.log(matricesIndex2[mIndex]) + internalPartialsIndex2[cIndex + cGenotypeIndex];
                    }

                    // max-sum
                    tmp1 = Math.log(matricesIndex1[mIndex]) + internalMLPartialsIndex1[cIndex + cGenotypeIndex];
                    if (has2ndChild) {
                        tmp2 = Math.log(matricesIndex2[mIndex]) + internalMLPartialsIndex2[cIndex + cGenotypeIndex];
                    }
                    if (cGenotypeIndex == 0) {
                        max1 = tmp1;
                        parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        if (has2ndChild) {
                            max2 = tmp2;
                            parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                        }
                    } else {
                        if (max1 < tmp1) {
                            max1 = tmp1;
                            parentMLGenotypes[pIndex][MLChild1Index].clear();
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        } else if (max1 == tmp1) {
                            parentMLGenotypes[pIndex][MLChild1Index].add(cGenotypeIndex);
                        }

                        if (has2ndChild) {
                            if (max2 < tmp2) {
                                max2 = tmp2;
                                parentMLGenotypes[pIndex][MLChild2Index].clear();
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            } else if (max2 == tmp2) {
                                parentMLGenotypes[pIndex][MLChild2Index].add(cGenotypeIndex);
                            }
                        }
                    }

                    mIndex++;
                }

                if (has2ndChild) {
                    // constant site
                    parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cst1) + MathFunctions.logSumExp(cst2);

                    // sum-product
                    parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1) + MathFunctions.logSumExp(sp2);

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1 + max2;
                } else {
                    // constant site
                    parentConstPartialsIndex[pIndex] = MathFunctions.logSumExp(cst1);

                    // sum-product
                    parentPartialsIndex[pIndex] = MathFunctions.logSumExp(sp1);

                    // max-sum
                    parentMLPartialsIndex[pIndex] = max1;
                }

                pIndex++;
            }
        }
    } // calculateLogPartialPartialPruning


    //**********************************************
    //*            Distribution methods            *
    //**********************************************

    @Override
    public void unstore() {
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes - 1);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);
    } // unstore

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes - 1);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);
    } // store


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    public int[] getCurrentMatrixIndex() {
        return currentMatrixIndex;
    } // getCurrentMatrixIndex

    public int[] getCurrentPartialsIndex() {
        return currentPartialsIndex;
    } // getCurrentPartialsIndex

}