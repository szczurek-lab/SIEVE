package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.ScsRandomLocalClockModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.rawreadcountsmodel.*;
import beast.evolution.sitemodel.ScsSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.ScsTree;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.variantsinfo.GenericVariantsInfo;
import beast.evolution.variantsinfo.VariantsInfoLog;
import beast.evolution.variantsinfo.VariantsInfoVCF;
import beast.math.util.MathFunctions;
import beast.util.ElementComparator;
import com.google.common.primitives.Ints;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see " +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class ScsTreeLikelihood extends ScsGenericTreeLikelihood {


    //**********************************************
    //*                   Inputs                   *
    //**********************************************

    public enum Scaling {none, always, _default}

    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " +
            Arrays.toString(ScsTreeLikelihood.Scaling.values()) +
            ". If not specified, the -beagle_scaling flag is used.",
            ScsTreeLikelihood.Scaling._default, ScsTreeLikelihood.Scaling.values());


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    /**
     * calculation engine *
     */
    protected LikelihoodCore likelihoodCore;

    /*
     * beagle likelihood core does not contain max-sum algorithm
     * extend it in the future
     */
    // protected BeagleScsTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected ScsSubstitutionModelBase substitutionModel;
    protected RawReadCountsModelInterface.Base rawReadCountsModel;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * whether to update partial likelihoods of leaves or not
     * should be true only when parameters in raw read counts model are being explored.
     */
    protected boolean updateLeaves;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;

    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;

    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

    protected int matrixSize;

    /**
     * Whether to tract the maximum likelihood genotypes or not.
     * Only tracing when in variant calling mode or the allelic sequencing coverage and its raw variance
     * in raw read counts model should be updated during not MCMC but post-processing.
     */
    protected boolean traceMLGenotypes;

    /**
     * [2] * [#nodes] * [#matrices * #patterns * #states] * [2 children]
     * Indicate the children's maximum likelihood genotype of a parent node.
     * <p>
     * For leaf node, `child 0` is used to store the index to the component
     * with the largest likelihood, therefore indicating the status of ADO.
     * <p>
     * For a pattern, a matrix, and a state.
     * Following a tree structure.
     */
    protected List<Integer>[][][][] MLGenotypesNodes;

    /**
     * the latest results of MLGenotypesNodes are retained
     */
    protected int[] currentMLGenotypesNodeIndex;
    protected int[] storedMLGenotypesNodeIndex;

    /**
     * used to update allelic sequencing coverage and its raw variance in the raw read counts model.
     * considering the weights for all possible maximum likelihood number of sequenced alleles.
     * [#matrices * #patterns]
     * all nodes (index from the number of external nodes) + the number of external nodes (for the number of sequenced alleles)
     */
    protected List<int[]>[] maxLikelihoodGenotypes;

    /**
     * store maximum likelihood genotypes for all nodes, all patterns, and all matrices
     * for a specific maximum likelihood path
     * <p>
     * [#matrices * #patterns * #nodes]
     */
    protected List<Integer>[] MLGenotypesCollection;

    /**
     * record cleanness status of each node during traverse
     */
    protected boolean[] isClean;

    /**
     * record patterns which require an update of allelicSeqCov & allelicSeqRawVar in raw read counts model
     * patterns inside require a complete update of partial likelihoods for all leaf nodes
     * <p>
     * converted from changedPatternsSet as reading from List is quicker than from Set
     * <p>
     * [#matrices, #patterns]
     */
    protected List<List<Integer>> changedPatternsList;
    protected Set<List<Integer>> changedPatternsSet;

    /**
     * containing the same elements as they are in changedPatternsList but more efficient
     */
    protected int[][] changedPatterns;

    /**
     * for scaling use
     * store the first index of each pattern in changedPatterns
     */
    protected List<Integer> changedPatternsIndex;

    /**
     * record changed patterns for a clean internal node
     * <p>
     * [#matrices, #patterns]
     */
    protected Set<List<Integer>>[] changedPatternsPerNode;

    /**
     * store root likelihood for constant site if using ascertainment bias correction
     * #patterns
     */
    protected double[] constRoot;

    /**
     * store (scaled) root log-likelihood for constant site if using ascertainment bias correction
     * #patterns
     */
    protected double[] logConstRoot;

    /**
     * whether to perform run time analysis
     */
    protected boolean runTimeAnalysis;

    /**
     * run time information for different part of likelihood computation
     * <p>
     * traverse time (0), including transition probability matrix computation time (1),
     * leaf partial assignment time (2), likelihood core running time (3), integrate across matrices (4),
     * and compute log likelihood of each pattern (5).
     * <p>
     * if the move is accepted:
     * maximum likelihood genotypes update time (6)
     * raw read counts model update time (7)
     * partly tree traverse time (8, 9 -> 2, 10 -> 3, 11 -> 4, 12 -> 5)
     */
    protected long[] runTime;

    // some frequently used variables
    protected int nrOfNodes;
    protected int nrOfInternalNodes;
    protected int nrOfExternalNodes;
    protected int nrOfMatrices;
    protected int nrOfPatterns;
    protected int nrOfStates;

    /**
     * Could be 2 (no ADO, one ADO) or 3 (no ADO, one ADO, two ADOs)
     */
    protected int nrOfAdoStates;

    /**
     * When branch rate model is not strict molecular clock model (Relaxed Molecular Clock),
     * turning on Felsenstein's bias correction will force the entire tree to be dirty to avoid numeric issues.
     */
    protected boolean isForcingTreeDirtyRMC;

    /**
     * Whether the maximum likelihood genotypes and ADO states are updated or not.
     */
    private boolean MLGenotypesAndAdosUpdated;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public ScsTreeLikelihood(final int nrOfNodes) {
        this.nrOfNodes = nrOfNodes;
    }

    public ScsTreeLikelihood() {
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // check whether inputs are valid
        inputsSanityCheck();

        // beagle likelihood computation
        /*
        beagle = null;
        beagle = new BeagleScsTreeLikelihood();
        try {
            beagle.initByName(
                    "scsData", scsDataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(), "scaling", scaling.get().toString());
            if (beagle.beagle != null) {
                //a Beagle instance was found, so we use it
                return;
            }
        } catch (Exception e) {
            // ignore
        }

        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;
        */

        // variables initialization
        nrOfNodes = treeInput.get().getNodeCount();
        m_siteModel = siteModelInput.get();
        ScsAlignment alignment = scsDataInput.get();
        ((ScsSiteModel) m_siteModel).setDataType(alignment.getDataType());
        substitutionModel = (ScsSubstitutionModelBase) m_siteModel.substModelInput.get();
        rawReadCountsModel = rawReadCountsModelInput.get();

        if (useLogPartialsInput.get() != null)
            useLogPartials = useLogPartialsInput.get();
        else
            useLogPartials = true;

        variablesSanityCheck();

        if (branchRateModelInput.get() != null)
            branchRateModel = branchRateModelInput.get();
        else
            branchRateModel = new StrictClockModel();

        m_branchLengths = new double[nrOfNodes];
        storedBranchLengths = new double[nrOfNodes];
        nrOfStates = substitutionModel.getStateCount();
        nrOfAdoStates = rawReadCountsModel.getModeledAllelesSize();
        nrOfPatterns = alignment.getPatternCount();
        nrOfMatrices = m_siteModel.getCategoryCount();
        nrOfInternalNodes = treeInput.get().getInternalNodeCount();
        nrOfExternalNodes = alignment.getTaxonCount();
        useAscBiasCorrection = alignment.isAscBiasCorrection();
        String className = getClass().getSimpleName();
        patternLogLikelihoods = new double[nrOfPatterns];
        m_fRootPartials = new double[nrOfPatterns];
        matrixSize = (nrOfStates + 1) * (nrOfStates + 1);
        probabilities = new double[(nrOfStates + 1) * (nrOfStates + 1)];
        Arrays.fill(probabilities, 1.0);

        isClean = new boolean[nrOfNodes];
        Arrays.fill(isClean, true);

        inVariantCallingMode = rawReadCountsModel.isInVariantCallingMode();
        traceMLGenotypes = rawReadCountsModel.updateSeqCovModel() || inVariantCallingMode || traceMLGenotypesInput.get();
        isForcingTreeDirtyRMC = alignment.isForcingTreeDirtyRMC();

        if (traceMLGenotypes) {
            MLGenotypesNodes = new ArrayList[2][nrOfNodes][nrOfMatrices * nrOfPatterns * nrOfStates][2];
            maxLikelihoodGenotypes = new ArrayList[nrOfMatrices * nrOfPatterns];
            MLGenotypesCollection = new ArrayList[nrOfMatrices * nrOfPatterns * nrOfNodes];
            currentMLGenotypesNodeIndex = new int[nrOfNodes];
            storedMLGenotypesNodeIndex = new int[nrOfNodes];
        }

        constRoot = new double[nrOfPatterns];
        logConstRoot = new double[nrOfPatterns];

        // if critical statistics in raw read counts model have not been computed, compute them
        if (!rawReadCountsModel.isDeeplyInitialized()) {
            try {
                rawReadCountsModel.deeplyInitialize(
                        alignment,
                        nrOfMatrices,
                        nrOfStates,
                        substitutionModel.getModeledAlleles(),
                        useLogPartials
                );
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        if (runTimeAnalysisInput.get() != null)
            runTimeAnalysis = runTimeAnalysisInput.get();
        else
            runTimeAnalysis = false;

        if (runTimeAnalysis)
            runTime = new long[13];

        resetVariables();

        likelihoodCore = new ScsBeerLikelihoodCore(nrOfStates);
        initCore();

        if (inVariantCallingMode) {
            if (variantsInfoInput.get() == null)
                throw new IllegalArgumentException("Error! 'variantsInfo' is missing in variant calling mode.");

            if (useOnlyBranchLengthInput.get() != null)
                useOnlyBranchLength = useOnlyBranchLengthInput.get();

            if (!useOnlyBranchLength && !(treeInput.get() instanceof ScsTree))
                throw new IllegalArgumentException("Error! The input tree should be of type 'beast.evolution.tree.ScsTree'");

            if (meanRateInput.get() != null)
                meanRate = meanRateInput.get().getValue();

            variantsInfo = variantsInfoInput.get();

            if (!(variantsInfo instanceof VariantsInfoVCF))
                throw new IllegalArgumentException("Error! 'variantsInfo' should be of type 'beast.evolution.variantsinfo.VariantsInfoVCF'");

            variantsInfo.initialise(nrOfPatterns, substitutionModel);
        } else if (traceMLGenotypesInput.get())
            variantsInfo = new VariantsInfoLog(alignment, treeInput.get(), substitutionModel);

        // print startup messages via Log.print
        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
    } // initAndValidate

    public void inputsSanityCheck() {
        if (scsDataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount())
            throw new IllegalArgumentException("Error! The number of nodes in the tree does not match the number of sequences");

        if (!(treeInput.get() instanceof ScsTree))
            throw new IllegalArgumentException("Error! 'treeInput' should be of type ScsTree (" + this.getClass().getName() + ")");

        if (!(siteModelInput.get() instanceof ScsSiteModel))
            throw new IllegalArgumentException("Error! 'siteModelInput' should be of type ScsSiteModel (" +
                    this.getClass().getName() + ")");

        if (!(siteModelInput.get().substModelInput.get() instanceof ScsSubstitutionModelBase))
            throw new IllegalArgumentException("Error! 'substModel' in 'siteModelInput' should be of type " +
                    "ScsSubstitutionModelBase (" + this.getClass().getName() + ")");
    } // inputsSanityCheck

    /**
     * make sure that the types of substitutionModel and rawReadCountsModel are matched.
     */
    protected void variablesSanityCheck() {
        if ((substitutionModel instanceof ScsFiniteMuModel) && (!(rawReadCountsModel instanceof RawReadCountsModelFiniteMuBB)))
            throw new IllegalArgumentException("Error! 'rawReadCountsModelInput' should be of type 'ScsErrorModelFiniteMuBB', rather than " +
                    "of the provided type " + rawReadCountsModel.getClass().getName() + ", because the 'substModel' in " +
                    "'siteModelInput' is of type 'ScsFiniteMuModel'. (" + this.getClass().getName() + ")");

        if ((substitutionModel instanceof ScsFiniteMuExtendedModel) && (!(rawReadCountsModel instanceof RawReadCountsModelFiniteMuDM)))
            throw new IllegalArgumentException("Error! 'rawReadCountsModelInput' should be of type 'RawReadCountsModelFiniteMuDM', rather than " +
                    "of the provided type " + rawReadCountsModel.getClass().getName() + ", because the 'substModel' in " +
                    "'siteModelInput' is of type 'ScsFiniteMuExtendedModel'. (" + this.getClass().getName() + ")");

        if ((substitutionModel instanceof ScsFiniteMuDelModel) && (!(rawReadCountsModel instanceof RawReadCountsModelFiniteMuDelBB)))
            throw new IllegalArgumentException("Error! 'rawReadCountsModelInput' should be of type 'ScsErrorModelFiniteMuDelBB', rather than " +
                    "of the provided type " + rawReadCountsModel.getClass().getName() + ", because the 'substModel' in " +
                    "'siteModelInput' is of type 'ScsFiniteMuDelModel'. (" + this.getClass().getName() + ")");

        if ((substitutionModel instanceof ScsFiniteMuDelInsModel) && (!(rawReadCountsModel instanceof RawReadCountsModelFiniteMuDelInsBB)))
            throw new IllegalArgumentException("Error! 'rawReadCountsModelInput' should be of type 'ScsErrorModelFiniteMuDelInsBB', rather than " +
                    "of the provided type " + rawReadCountsModel.getClass().getName() + ", because the 'substModel' in " +
                    "siteModelInput is of type 'ScsFiniteMuDelInsModel'. (" + this.getClass().getName() + ")");
    } // variablesSanityCheck

    protected void initCore() {
        ((ScsBeerLikelihoodCore) likelihoodCore).initialize(
                nrOfNodes,
                nrOfExternalNodes,
                nrOfInternalNodes,
                nrOfPatterns,
                nrOfMatrices,
                nrOfStates,
                true,
                useLogPartials
        );

        final long startTime = System.currentTimeMillis();
        for (Node i : treeInput.get().getExternalNodes()) {
            final int nodeIndex = i.getNr();

            likelihoodCore.setNodePartials(
                    nodeIndex,
                    (traceMLGenotypes ?
                            rawReadCountsModel.initializeLeafLikelihood(
                                    i,
                                    MLGenotypesNodes
                            ) :
                            rawReadCountsModel.initializeLeafLikelihood(i)
                    )
            );

            ((ScsBeerLikelihoodCore) likelihoodCore).storeLeafPartials(nodeIndex);
        }
        final long endTime = System.currentTimeMillis();
        if (runTimeAnalysis)
            System.out.println("Initially set all leaf partials: " + (endTime - startTime) + " milliseconds.");

        hasDirt = Tree.IS_FILTHY;
        updateLeaves = true;

        for (int i = 0; i < nrOfInternalNodes; i++)
            likelihoodCore.createNodePartials(i + nrOfExternalNodes);
    } // initCore

    /**
     * Update leaf likelihoods for changed patterns.
     * Should ONLY be called during post processing.
     *
     * @param node leaf node
     */
    protected void updateLeafLikelihoods(final Node node) {
        assert node.isLeaf();
        final int nodeIndex = node.getNr();

        double[] partials = new double[nrOfMatrices * nrOfPatterns * nrOfStates];
        likelihoodCore.getNodePartials(nodeIndex, partials);

        rawReadCountsModel.updatePartialLeafLikelihoods(
                node,
                changedPatterns,
                MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                partials
        );

        // copy to likelihood core
        List<int[]> added = ((ScsBeerLikelihoodCore) likelihoodCore).setNodePartials(
                nodeIndex,
                partials,
                changedPatterns,
                changedPatternsIndex
        );
        if (added != null && added.size() > 0)
            updateChangedPatterns(added);
    } // updateLeafLikelihoods

    /**
     * get the corresponding branch rate.
     * <p>
     * If in phylogenetic inference mode, this works the same as the original beast.
     * Otherwise, it is acquired from the given tree and clock rate.
     *
     * @param node which node is being checked?
     * @return branch rate
     */
    private double getBranchRate(Node node) {
        if (!this.inVariantCallingMode) return branchRateModel.getRateForBranch(node);

        // in variant calling mode
        if (useOnlyBranchLength) {
            // only use branch length, not considering any rates provided by the input tree.
            // usually happening when the rate and branch length cannot be differentiated,
            // and the branch length in the input tree is the product of real time and rate.

            if (branchRateModel instanceof ScsRandomLocalClockModel &&
                    !((ScsRandomLocalClockModel) branchRateModel).isScaling())
                return 1.0;

            return meanRate;
        } else {
            // use both the branch length and rate provided by the input tree.
            // usually happening when the data containing real time information or other constraints,
            // such as mean rate and MRCA prior, is available.

            if (branchRateModel instanceof ScsRandomLocalClockModel &&
                    !((ScsRandomLocalClockModel) branchRateModel).isScaling())
                return ((ScsTree) treeInput.get()).getBranchRate(node);

            return meanRate * ((ScsTree) treeInput.get()).getBranchRate(node);
        }
    } // getBranchRate

    /**
     * traverse the tree and update changed values
     *
     * @param node      a node on the tree, normally root node
     * @param isPartial whether to compute partial likelihoods and maximum likelihood genotypes only for a part of patterns and matrices
     * @return whether the tree (root node as the argument) or subtree (other nodes as the argument) likelihood needs to be updated (>0) or not (=0)
     */
    @SuppressWarnings("deprecation")
    public int traverse(final Node node, final boolean isPartial) {
        if (node == null)
            return Tree.IS_CLEAN;

        if (node.isRoot())
            MLGenotypesAndAdosUpdated = false;

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        // branch-wise rate multiplier
        final double branchRate = getBranchRate(node);

        // (adjusted) number of mutations (or substitutions, depending on the interpretation of evolutionary model) per site
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex]) && !isPartial) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < nrOfMatrices; i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;

                final long startTime1 = System.currentTimeMillis();
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                final long endTime1 = System.currentTimeMillis();
                if (runTimeAnalysis)
                    runTime[1] += (endTime1 - startTime1);

                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is a leaf, update partials if it is in post processing stage or a parameter of raw read counts model is being explored or computing tree likelihood at the first time
        if (node.isLeaf() && updateLeaves) {
            // if in post processing stage, update the index of current partial likelihoods
            if (!isPartial)
                setIndexForUpdate(nodeIndex);

            final long startTime2 = System.currentTimeMillis();
            if (isPartial)
                updateLeafLikelihoods(node);
            else
                likelihoodCore.setNodePartials(
                        nodeIndex,
                        (traceMLGenotypes ?
                                rawReadCountsModel.computeLeafLikelihood(
                                        node,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex]
                                ) :
                                rawReadCountsModel.computeLeafLikelihood(node)
                        )
                );
            final long endTime2 = System.currentTimeMillis();
            if (runTimeAnalysis) {
                if (isPartial)
                    runTime[9] += (endTime2 - startTime2);
                else
                    runTime[2] += (endTime2 - startTime2);
            }
        }

        // If the node is internal, update the partial likelihoods
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1, isPartial);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2, isPartial);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                // if in post processing stage, update contents of current partial likelihoods and maximum likelihood genotypes
                if (!isPartial)
                    setIndexForUpdate(nodeIndex);

                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY)
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);

                final int childNum1 = child1.getNr();
                final long startTime3 = System.currentTimeMillis();
                if (node.getChildCount() == 2) {
                    final int childNum2 = child2.getNr();

                    // calculate partial likelihoods by integrating over all categories

                    if (isPartial) {
                        // at post processing stage

                        if (rawReadCountsModel.updateSeqCovModel()) {
                            // combinations of maximum likelihood number of sequenced alleles for leaves are dynamic
                            // sum-product and max-sum

                            List<int[]> added;

                            if (useLogPartials)
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                                        substitutionModel.getConstGenotype(),
                                        changedPatterns,
                                        changedPatternsIndex
                                );
                            else
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                                        substitutionModel.getConstGenotype(),
                                        changedPatterns,
                                        changedPatternsIndex
                                );

                            if (added != null && added.size() > 0)
                                updateChangedPatterns(added);

                        } else {
                            // the combination of maximum likelihood number of sequenced alleles for leaves is deterministic
                            // sum-product

                            List<int[]> added;

                            if (useLogPartials)
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        substitutionModel.getConstGenotype(),
                                        changedPatterns,
                                        changedPatternsIndex
                                );
                            else
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        substitutionModel.getConstGenotype(),
                                        changedPatterns,
                                        changedPatternsIndex
                                );

                            if (added != null && added.size() > 0)
                                updateChangedPatterns(added);

                        }

                    } else {
                        // at tree likelihood computation stage

                        if (traceMLGenotypes) {
                            // combinations of maximum likelihood number of sequenced alleles for leaves are dynamic
                            // sum-product and max-sum
                            // constant site partials are computed at the same time

                            if (useLogPartials)
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                                        substitutionModel.getConstGenotype()
                                );
                            else
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                                        substitutionModel.getConstGenotype()
                                );

                        } else {
                            // the combination of maximum likelihood number of sequenced alleles for leaves is deterministic
                            // sum-product
                            // constant site partials are computed at the same time

                            if (useLogPartials)
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        substitutionModel.getConstGenotype()
                                );
                            else
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        child1.isLeaf(),
                                        childNum2,
                                        child2.isLeaf(),
                                        nodeIndex,
                                        substitutionModel.getConstGenotype()
                                );
                        }

                    }

                } else {
                    // calculate partial likelihoods by integrating over all categories

                    if (isPartial) {
                        // at post processing stage

                        if (rawReadCountsModel.updateSeqCovModel()) {
                            // combinations of maximum likelihood number of sequenced alleles for leaves are dynamic
                            // sum-product and max-sum

                            List<int[]> added;

                            if (useLogPartials)
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                                        changedPatterns, changedPatternsIndex
                                );
                            else
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex],
                                        changedPatterns, changedPatternsIndex
                                );

                            if (added != null && added.size() > 0)
                                updateChangedPatterns(added);

                        } else {
                            // the combination of maximum likelihood number of sequenced alleles for leaves is deterministic
                            // sum-product

                            List<int[]> added;

                            if (useLogPartials)
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        false,
                                        -1,
                                        false,
                                        nodeIndex,
                                        substitutionModel.getConstGenotype(),
                                        changedPatterns, changedPatternsIndex
                                );
                            else
                                added = ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        false,
                                        -1,
                                        false,
                                        nodeIndex,
                                        substitutionModel.getConstGenotype(),
                                        changedPatterns, changedPatternsIndex
                                );

                            if (added != null && added.size() > 0)
                                updateChangedPatterns(added);

                        }

                    } else {
                        // at tree likelihood computation stage

                        if (traceMLGenotypes) {
                            // combinations of maximum likelihood number of sequenced alleles for leaves are dynamic
                            // sum-product and max-sum
                            // constant site partials are computed at the same time

                            if (useLogPartials)
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex]
                                );
                            else
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        nodeIndex,
                                        MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex]
                                );
                        } else {
                            // the combination of maximum likelihood number of sequenced alleles for leaves is deterministic
                            // sum-product
                            // constant site partials are computed at the same time

                            if (useLogPartials)
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogPartials(
                                        childNum1,
                                        false,
                                        -1,
                                        false,
                                        nodeIndex,
                                        substitutionModel.getConstGenotype()
                                );
                            else
                                ((ScsBeerLikelihoodCore) likelihoodCore).calculatePartials(
                                        childNum1,
                                        false,
                                        -1,
                                        false,
                                        nodeIndex,
                                        substitutionModel.getConstGenotype()
                                );
                        }
                    }

                }
                final long endTime3 = System.currentTimeMillis();
                if (runTimeAnalysis) {
                    if (isPartial)
                        runTime[10] += (endTime3 - startTime3);
                    else
                        runTime[3] += (endTime3 - startTime3);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree
                    // The genotype of the root is fixed to a specific genotype

                    // integrate across all site categories
                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    final long startTime4 = System.currentTimeMillis();
                    ((ScsBeerLikelihoodCore) likelihoodCore).integratePartials(
                            node.getNr(),
                            proportions,
                            substitutionModel.getRootGenotype(),
                            m_fRootPartials,
                            constRoot
                    );
                    final long endTime4 = System.currentTimeMillis();
                    if (runTimeAnalysis) {
                        if (isPartial)
                            runTime[11] += (endTime4 - startTime4);
                        else
                            runTime[4] += (endTime4 - startTime4);
                    }

                    // find out to which category each pattern belongs
                    if (this.inVariantCallingMode || traceMLGenotypesInput.get())
                        ((ScsBeerLikelihoodCore) likelihoodCore).getPatternCategories(
                                substitutionModel.getRootGenotype(),
                                node.getNr(),
                                variantsInfo
                        );

                    // compute log likelihood for each pattern
                    final long startTime5 = System.currentTimeMillis();
                    ((ScsBeerLikelihoodCore) likelihoodCore).calculateLogLikelihoods(
                            m_fRootPartials,
                            patternLogLikelihoods,
                            constRoot,
                            logConstRoot
                    );
                    final long endTime5 = System.currentTimeMillis();
                    if (runTimeAnalysis) {
                        if (isPartial)
                            runTime[12] += (endTime5 - startTime5);
                        else
                            runTime[5] += (endTime5 - startTime5);
                    }
                }
            }
        }

        isClean[nodeIndex] = (update == Tree.IS_CLEAN);

        return update;
    } // traverseWithBRM

    /**
     * get maximum likelihood genotypes for a node
     * compute maximum likelihood number of sequenced alleles for leaf nodes for all patterns
     * <p>
     * if there are more than one group of maximum likelihood number of sequenced alleles found,
     * compute the patterns and take their weights into account when updating allelic sequencing
     * coverage and raw variance in raw read counts model.
     *
     * @param node node to be processed
     */
    protected void getMLGenotypes(final Node node) {
        final int nodeIndex = node.getNr();

        if (node.isRoot()) {
            MLGenotypesAndAdosUpdated = true;

            // for root node, initialize the maximum likelihood genotypes' data structure

            // reset some variables
            resetVariables();

            if (!isClean[nodeIndex]) {

                final int rootGenotype = substitutionModel.getRootGenotype();

                // nIndex = matrixIndex * nrOfPatterns * nrOfNodes + patternIndex * nrOfNodes
                // mpIndex = matrixIndex * nrOfPatterns + patternIndex
                // mlIndex = matrixIndex * nrOfPatterns * nrOfStates + patternIndex * nrOfStates + rootGenotype
                int nIndex, mpIndex, mlIndex;
                nIndex = mpIndex = 0;
                mlIndex = rootGenotype;

                // it is necessary to loop over every site for every matrix for the tree root
                for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

                    for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                        // update maxLikelihoodGenotypes for nodes having at least two children before proceeding
                        // this is ESSENTIAL as tree structure movements could change the parent-child relationships
                        if (node.getChildCount() > 1 &&
                                maxLikelihoodGenotypes[mpIndex] != null &&
                                maxLikelihoodGenotypes[mpIndex].size() > 1)
                            updateMaxLikelihoodGenotypes(
                                    maxLikelihoodGenotypes[mpIndex],
                                    node,
                                    rootGenotype
                            );

                        for (int childIndex = 0; childIndex < node.getChildCount(); childIndex++) {

                            final int childNodeIndex = node.getChild(childIndex).getNr();

                            // only for the first time run
                            if (MLGenotypesCollection[nIndex + nodeIndex] == null) {
                                changedPatternsPerNode[nodeIndex].add(Stream.of(matrixIndex, patternIndex).collect(Collectors.toList()));

                                MLGenotypesCollection[nIndex + nodeIndex] = new ArrayList<>();
                                MLGenotypesCollection[nIndex + nodeIndex].add(rootGenotype);

                                maxLikelihoodGenotypes[mpIndex] = new ArrayList<>();
                                maxLikelihoodGenotypes[mpIndex].add(new int[nrOfNodes + nrOfExternalNodes]);
                                Arrays.fill(maxLikelihoodGenotypes[mpIndex].get(0), -1);
                                maxLikelihoodGenotypes[mpIndex].get(0)[nodeIndex + nrOfExternalNodes] = rootGenotype;
                            }

                            if (MLGenotypesCollection[nIndex + childNodeIndex] == null) {
                                changedPatternsPerNode[childNodeIndex].add(Stream.of(matrixIndex, patternIndex).collect(Collectors.toList()));

                                MLGenotypesCollection[nIndex + childNodeIndex] = new ArrayList<>();
                            }

                            // get child's maximum likelihood genotype for this pattern and matrix
                            List<Integer> childGenotypes = MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex][mlIndex][childIndex];
                            if (childGenotypes.size() < 1)
                                throw new IllegalArgumentException("No maximum likelihood genotype for children " +
                                        "found! (" + this.getClass().getName() + ")");

                            // update if elements are not the same, otherwise nothing needs to be done
                            if (!childGenotypes.equals(MLGenotypesCollection[nIndex + childNodeIndex])) {
                                changedPatternsPerNode[childNodeIndex].add(Stream.of(matrixIndex, patternIndex).collect(Collectors.toList()));

                                List<Integer> deleted = new ArrayList<>();
                                List<Integer> added = new ArrayList<>();
                                getDeletedAndAddedGenotypes(
                                        MLGenotypesCollection[nIndex + childNodeIndex],
                                        deleted,
                                        childGenotypes,
                                        added
                                );

                                updateCollectionAndGenotypes(
                                        MLGenotypesCollection[nIndex + childNodeIndex],
                                        maxLikelihoodGenotypes[mpIndex],
                                        childNodeIndex + nrOfExternalNodes,
                                        deleted,
                                        added
                                );
                            }
                        }

                        nIndex += nrOfNodes;
                        mpIndex++;
                        mlIndex += nrOfStates;
                    }
                }

                // go downwards
                for (int childIndex = 0; childIndex < node.getChildCount(); childIndex++)
                    getMLGenotypes(node.getChild(childIndex));

                if (rawReadCountsModel.updateSeqCovModel()) {
                    // Updating allelic sequencing coverage and raw variance.

                    // gather changed patterns of every node together
                    for (int i = 0; i < nrOfNodes; i++)
                        changedPatternsSet.addAll(changedPatternsPerNode[i]);

                    changedPatternsList.addAll(changedPatternsSet);

                    /*
                     * store the maximum likelihood number of alleles for leaves
                     * only for the changed patterns and matrices
                     * do not use changedPatternsSet from now
                     */
                    Iterator<List<Integer>> itr = changedPatternsList.iterator();
                    while (itr.hasNext()) {
                        List<Integer> pair = itr.next();
                        mpIndex = pair.get(0) * nrOfPatterns + pair.get(1);

                        ListComparator comparator = new ListComparator();

                        /*
                         * convert maximum likelihood genotypes of all leaf nodes for this pattern and matrix
                         * to maximum likelihood number of sequenced alleles and store
                         */
                        List<int[]> patternMLNrOfSequencedAlleles = new ArrayList<>();

                        for (int combIndex = 0; combIndex < maxLikelihoodGenotypes[mpIndex].size(); combIndex++) {
                            patternMLNrOfSequencedAlleles.add(new int[nrOfExternalNodes]);

                            for (int taxonIndex = 0; taxonIndex < nrOfExternalNodes; taxonIndex++) {
                                final int[] pat = maxLikelihoodGenotypes[mpIndex].get(combIndex);

                                patternMLNrOfSequencedAlleles.get(patternMLNrOfSequencedAlleles.size() - 1)[taxonIndex] = substitutionModel.getNrOfAlleles(pat[taxonIndex + nrOfExternalNodes]) - pat[taxonIndex];
                            }
                        }

                        // sort the list of integer arrays
                        if (patternMLNrOfSequencedAlleles.size() > 1)
                            patternMLNrOfSequencedAlleles.sort(comparator);

                        // compute the weights of integer arrays
                        int patternsAlleles = 1;
                        int[] weightsAlleles = new int[patternMLNrOfSequencedAlleles.size()];
                        weightsAlleles[0] = 1;
                        List<int[]> finalPatternMLNrOfSequencedAlleles = new ArrayList<>();
                        finalPatternMLNrOfSequencedAlleles.add(patternMLNrOfSequencedAlleles.get(0));
                        for (int i = 1; i < patternMLNrOfSequencedAlleles.size(); i++) {
                            if (comparator.compare(patternMLNrOfSequencedAlleles.get(i - 1), patternMLNrOfSequencedAlleles.get(i)) != 0) {
                                patternsAlleles++;
                                finalPatternMLNrOfSequencedAlleles.add(patternMLNrOfSequencedAlleles.get(i));
                            }
                            weightsAlleles[patternsAlleles - 1]++;
                        }

                        // compare the difference from the previous calculation
                        boolean needToCopy = true;
                        List<int[]> previousPatternMLNrOfSequencedAlleles = rawReadCountsModel.getMLNrOfSequencedAllelesPatterns(mpIndex);
                        if (previousPatternMLNrOfSequencedAlleles != null && previousPatternMLNrOfSequencedAlleles.equals(finalPatternMLNrOfSequencedAlleles)) {
                            if (previousPatternMLNrOfSequencedAlleles.size() == 1) {
                                itr.remove();
                                needToCopy = false;
                            } else {
                                for (int i = 0; i < patternsAlleles; i++) {
                                    int previousWeight = rawReadCountsModel.getMLNrOfSequencedAllelesWeights(mpIndex, i);
                                    int finalWeight = weightsAlleles[finalPatternMLNrOfSequencedAlleles.indexOf(previousPatternMLNrOfSequencedAlleles.get(i))];

                                    // check whether the weights match or not
                                    if (previousWeight != finalWeight)
                                        break;

                                    if (i == patternsAlleles - 1) {
                                        itr.remove();
                                        needToCopy = false;
                                    }
                                }
                            }
                        }

                        // store information to rawReadCountsModel instance if necessary
                        if (needToCopy) {
                            if (previousPatternMLNrOfSequencedAlleles != null)
                                rawReadCountsModel.resetMLNrOfSequencedAllelesPatternsAndWeights(mpIndex);

                            for (int i = 0; i < patternsAlleles; i++) {
                                rawReadCountsModel.addMLNrOfExistingAllelesPatterns(mpIndex, finalPatternMLNrOfSequencedAlleles.get(i));
                                rawReadCountsModel.addMLNrOfSequencedAllelesWeights(mpIndex, weightsAlleles[i]);
                            }
                        }
                    }

                }

            }
        } else if (!node.isLeaf()) {
            // for internal but not root nodes

            // if the node is not clean, then compute the maximum likelihood genotypes path for every pattern and matrix
            if (!isClean[nodeIndex]) {

                // nIndex = matrixIndex * nrOfPatterns * nrOfNodes + patternIndex * nrOfNodes
                // mpIndex = matrixIndex * nrOfPatterns + patternIndex
                // mlIndex = matrixIndex * nrOfPatterns * nrOfStates + patternIndex * nrOfStates
                int nIndex, mpIndex, mlIndex;
                nIndex = mpIndex = mlIndex = 0;

                for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

                    for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                        final List<Integer> genotypes = MLGenotypesCollection[nIndex + nodeIndex];

                        // update maxLikelihoodGenotypes for nodes having at least two children before proceeding
                        // this is ESSENTIAL as tree structure movements could change the parent-child relationships
                        if (node.getChildCount() > 1 && maxLikelihoodGenotypes[mpIndex].size() > 1)
                            updateMaxLikelihoodGenotypes(
                                    maxLikelihoodGenotypes[mpIndex],
                                    node,
                                    genotypes
                            );

                        int pFlag = 0;
                        for (int genotypeIndex : genotypes) {

                            // loop over two children
                            for (int childIndex = 0; childIndex < node.getChildCount(); childIndex++) {

                                // get node index of child
                                int childNodeIndex = node.getChild(childIndex).getNr();

                                // only for the first time run
                                if (MLGenotypesCollection[nIndex + childNodeIndex] == null) {
                                    changedPatternsPerNode[childNodeIndex].add(Stream.of(matrixIndex, patternIndex).collect(Collectors.toList()));

                                    MLGenotypesCollection[nIndex + childNodeIndex] = new ArrayList<>();
                                }

                                List<Integer> previousChildGenotypes = getMLGenotypeEntries(
                                        maxLikelihoodGenotypes[mpIndex],
                                        nodeIndex + nrOfExternalNodes,
                                        genotypeIndex,
                                        childNodeIndex + nrOfExternalNodes
                                );

                                //assert previousChildGenotypes.size() > 0;

                                List<Integer> childGenotypes = MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex][mlIndex + genotypeIndex][childIndex];
                                if (childGenotypes.size() < 1)
                                    throw new IllegalArgumentException("No maximum likelihood genotype for children found! (" +
                                            this.getClass().getName() + ")");

                                if (!childGenotypes.equals(previousChildGenotypes)) {
                                    changedPatternsPerNode[childNodeIndex].add(Stream.of(matrixIndex, patternIndex).collect(Collectors.toList()));

                                    List<Integer> deleted = new ArrayList<>();
                                    List<Integer> added = new ArrayList<>();
                                    getDeletedAndAddedGenotypes(
                                            previousChildGenotypes,
                                            deleted,
                                            childGenotypes,
                                            added
                                    );

                                    updateCollectionAndGenotypes(
                                            previousChildGenotypes,
                                            maxLikelihoodGenotypes[mpIndex],
                                            nodeIndex + nrOfExternalNodes,
                                            genotypeIndex,
                                            childNodeIndex + nrOfExternalNodes,
                                            deleted,
                                            added
                                    );
                                }

                                // update MLGenotypesCollection
                                if (pFlag == 0)
                                    MLGenotypesCollection[nIndex + childNodeIndex].clear();

                                for (int i : previousChildGenotypes)
                                    if (!MLGenotypesCollection[nIndex + childNodeIndex].contains(i))
                                        MLGenotypesCollection[nIndex + childNodeIndex].add(i);
                            }

                            pFlag++;
                        }

                        nIndex += nrOfNodes;
                        mpIndex++;
                        mlIndex += nrOfStates;
                    }
                }
            }

            // if the node is clean, but some patterns and matrices should be recomputed
            if (isClean[nodeIndex] && changedPatternsPerNode[nodeIndex].size() > 0) {

                // nIndex = matrixIndex * nrOfPatterns * nrOfNodes + patternIndex * nrOfNodes
                // mpIndex = matrixIndex * nrOfPatterns + patternIndex
                // mlIndex = matrixIndex * nrOfPatterns * nrOfStates + patternIndex * nrOfStates
                int nIndex, mpIndex, mlIndex;

                for (List<Integer> pair : changedPatternsPerNode[nodeIndex]) {

                    nIndex = pair.get(0) * nrOfPatterns * nrOfNodes + pair.get(1) * nrOfNodes;
                    mpIndex = pair.get(0) * nrOfPatterns + pair.get(1);
                    mlIndex = pair.get(0) * nrOfPatterns * nrOfStates + pair.get(1) * nrOfStates;

                    final List<Integer> genotypes = MLGenotypesCollection[nIndex + nodeIndex];

                    // update maxLikelihoodGenotypes for nodes having at least two children before proceeding
                    // this is ESSENTIAL as tree structure movements could change the parent-child relationships
                    if (node.getChildCount() > 1 && maxLikelihoodGenotypes[mpIndex].size() > 1)
                        updateMaxLikelihoodGenotypes(
                                maxLikelihoodGenotypes[mpIndex],
                                node,
                                genotypes
                        );

                    int pFlag = 0;
                    for (int genotypeIndex : genotypes) {

                        // loop over two children
                        for (int childIndex = 0; childIndex < node.getChildCount(); childIndex++) {

                            // get node index of child
                            int childNodeIndex = node.getChild(childIndex).getNr();

                            // only for the first time run
                            if (MLGenotypesCollection[nIndex + childNodeIndex] == null) {
                                changedPatternsPerNode[childNodeIndex].add(new ArrayList<>(pair));

                                MLGenotypesCollection[nIndex + childNodeIndex] = new ArrayList<>();
                            }

                            List<Integer> previousChildGenotypes = getMLGenotypeEntries(
                                    maxLikelihoodGenotypes[mpIndex],
                                    nodeIndex + nrOfExternalNodes,
                                    genotypeIndex,
                                    childNodeIndex + nrOfExternalNodes
                            );

//                            assert previousChildGenotypes.size() > 0;

                            List<Integer> childGenotypes = MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex][mlIndex + genotypeIndex][childIndex];
                            if (childGenotypes.size() < 1)
                                throw new IllegalArgumentException("No maximum likelihood genotype for children found! (" +
                                        this.getClass().getName() + ")");

                            if (!childGenotypes.equals(previousChildGenotypes)) {
                                changedPatternsPerNode[childNodeIndex].add(new ArrayList<>(pair));

                                List<Integer> deleted = new ArrayList<>();
                                List<Integer> added = new ArrayList<>();
                                getDeletedAndAddedGenotypes(
                                        previousChildGenotypes,
                                        deleted,
                                        childGenotypes,
                                        added
                                );

                                updateCollectionAndGenotypes(
                                        previousChildGenotypes,
                                        maxLikelihoodGenotypes[mpIndex],
                                        nodeIndex + nrOfExternalNodes,
                                        genotypeIndex,
                                        childNodeIndex + nrOfExternalNodes,
                                        deleted,
                                        added
                                );
                            }

                            // update MLGenotypesCollection
                            if (pFlag == 0)
                                MLGenotypesCollection[nIndex + childNodeIndex].clear();

                            for (int i : previousChildGenotypes)
                                if (!MLGenotypesCollection[nIndex + childNodeIndex].contains(i))
                                    MLGenotypesCollection[nIndex + childNodeIndex].add(i);
                        }

                        pFlag++;
                    }
                }
            }

            // go downwards
            for (int childIndex = 0; childIndex < node.getChildCount(); childIndex++)
                getMLGenotypes(node.getChild(childIndex));

        } else if (node.isLeaf()) {
            // for leaf nodes find out if ADO is likely to happen or not / call ADO states

            // if the node is not clean, then compute the maximum likelihood number of sequenced alleles path for every pattern and matrix
            if (!isClean[nodeIndex]) {

                // nIndex = matrixIndex * nrOfPatterns * nrOfNodes + patternIndex * nrOfNodes
                // mpIndex = matrixIndex * nrOfPatterns + patternIndex
                // mlIndex = matrixIndex * nrOfPatterns * nrOfStates + patternIndex * nrOfStates
                int nIndex, mpIndex, mlIndex;
                nIndex = mpIndex = mlIndex = 0;

                for (int matrixIndex = 0; matrixIndex < nrOfMatrices; matrixIndex++) {

                    for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

                        final List<Integer> genotypes = MLGenotypesCollection[nIndex + nodeIndex];

                        for (int genotypeIndex : genotypes) {

                            List<Integer> previousSequencedAllelesIndices = getMLGenotypeEntries(
                                    maxLikelihoodGenotypes[mpIndex],
                                    nodeIndex + nrOfExternalNodes,
                                    genotypeIndex,
                                    nodeIndex
                            );

                            List<Integer> sequencedAllelesIndices = MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex][mlIndex + genotypeIndex][0];
                            assert sequencedAllelesIndices != null;
                            if (sequencedAllelesIndices.size() < 1)
                                throw new IllegalArgumentException("No maximum likelihood number of sequenced " +
                                        "alleles for a leaf node found! (" + this.getClass().getName() + ")");

                            if (!sequencedAllelesIndices.equals(previousSequencedAllelesIndices)) {
                                changedPatternsPerNode[nodeIndex].add(Stream.of(matrixIndex, patternIndex).collect(Collectors.toList()));

                                List<Integer> deleted = new ArrayList<>();
                                List<Integer> added = new ArrayList<>();
                                getDeletedAndAddedGenotypes(
                                        previousSequencedAllelesIndices,
                                        deleted,
                                        sequencedAllelesIndices,
                                        added
                                );

                                updateCollectionAndGenotypes(
                                        previousSequencedAllelesIndices,
                                        maxLikelihoodGenotypes[mpIndex],
                                        nodeIndex + nrOfExternalNodes,
                                        genotypeIndex,
                                        nodeIndex,
                                        deleted,
                                        added
                                );
                            }
                        }

                        nIndex += nrOfNodes;
                        mpIndex++;
                        mlIndex += nrOfStates;
                    }
                }
            }

            // if the node is clean, but some patterns and matrices should be recomputed
            if (isClean[nodeIndex] && changedPatternsPerNode[nodeIndex].size() > 0) {

                // nIndex = matrixIndex * nrOfPatterns * nrOfNodes + patternIndex * nrOfNodes
                // mpIndex = matrixIndex * nrOfPatterns + patternIndex
                // mlIndex = matrixIndex * nrOfPatterns * nrOfStates + patternIndex * nrOfStates
                int nIndex, mpIndex, mlIndex;

                for (List<Integer> pair : changedPatternsPerNode[nodeIndex]) {

                    nIndex = pair.get(0) * nrOfPatterns * nrOfNodes + pair.get(1) * nrOfNodes;
                    mpIndex = pair.get(0) * nrOfPatterns + pair.get(1);
                    mlIndex = pair.get(0) * nrOfPatterns * nrOfStates + pair.get(1) * nrOfStates;

                    final List<Integer> genotypes = MLGenotypesCollection[nIndex + nodeIndex];

                    for (int genotypeIndex : genotypes) {

                        List<Integer> previousSequencedAllelesIndices = getMLGenotypeEntries(
                                maxLikelihoodGenotypes[mpIndex],
                                nodeIndex + nrOfExternalNodes,
                                genotypeIndex,
                                nodeIndex
                        );

//                        assert previousSequencedAllelesIndices.size() > 0;

                        List<Integer> sequencedAllelesIndices = MLGenotypesNodes[currentMLGenotypesNodeIndex[nodeIndex]][nodeIndex][mlIndex + genotypeIndex][0];
                        if (sequencedAllelesIndices.size() < 1)
                            throw new IllegalArgumentException("No maximum likelihood genotype for children found! (" +
                                    this.getClass().getName() + ")");

                        if (!sequencedAllelesIndices.equals(previousSequencedAllelesIndices)) {
                            changedPatternsPerNode[nodeIndex].add(new ArrayList<>(pair));

                            List<Integer> deleted = new ArrayList<>();
                            List<Integer> added = new ArrayList<>();
                            getDeletedAndAddedGenotypes(
                                    previousSequencedAllelesIndices,
                                    deleted,
                                    sequencedAllelesIndices,
                                    added
                            );

                            updateCollectionAndGenotypes(
                                    previousSequencedAllelesIndices,
                                    maxLikelihoodGenotypes[mpIndex],
                                    nodeIndex + nrOfExternalNodes,
                                    genotypeIndex,
                                    nodeIndex,
                                    deleted,
                                    added
                            );
                        }
                    }
                }
            }

        }

    } // getMLGenotypes

    private void updateMaxLikelihoodGenotypes(
            List<int[]> MLGenotypes,
            final Node node,
            final List<Integer> genotypes
    ) {
        for (final int genotype : genotypes)
            updateMaxLikelihoodGenotypes(MLGenotypes, node, genotype);
    } // updateMaxLikelihoodGenotypes

    private void updateMaxLikelihoodGenotypes(
            List<int[]> MLGenotypes,
            final Node node,
            final int genotype
    ) {

        final int childCount = node.getChildCount();

        // store node numbers of every subtree of this node in pre-order
        List<Integer>[] subtreeLabels = new ArrayList[childCount];

        // store the unique root patterns excluding subtrees (subtrees marked with -1)
        List<int[]> uniqueRootPatterns = new ArrayList<>();

        // store the unique patterns of subtrees (in pre-order)
        List<int[]>[] uniqueSubtreePatterns = new ArrayList[childCount];

        // get indices in MLGenotypes which have "genotype" for "node"
        final List<Integer> pIndices = getMLGenotypeIndices(MLGenotypes, node.getNr() + nrOfExternalNodes, genotype);

        // the number of combinations when the node has some other genotypes
        final int restSize = MLGenotypes.size() - pIndices.size();

        // traverse subtrees and get node labels
        for (int i = 0; i < childCount; i++) {
            if (subtreeLabels[i] == null)
                subtreeLabels[i] = new ArrayList<>();

            traverseSubtree(node.getChild(i), subtreeLabels[i]);
        }

        // get unique root patterns
        for (int i : pIndices) {
            int[] patternCopy = new int[nrOfNodes + nrOfExternalNodes];
            System.arraycopy(MLGenotypes.get(i), 0, patternCopy, 0, nrOfNodes + nrOfExternalNodes);

            for (List<Integer> j : subtreeLabels) {
                for (int k : j) {
                    patternCopy[k + nrOfExternalNodes] = -1;

                    if (k < nrOfExternalNodes)
                        patternCopy[k] = -1;
                }
            }

            if (!ElementComparator.compIntArrays(uniqueRootPatterns, patternCopy))
                uniqueRootPatterns.add(patternCopy);
        }

        // get unique subtree patterns
        for (int i : pIndices) {
            for (int j = 0; j < childCount; j++) {
                if (uniqueSubtreePatterns[j] == null)
                    uniqueSubtreePatterns[j] = new ArrayList<>();

                int[] patternCopy = new int[nrOfNodes + nrOfExternalNodes];
                Arrays.fill(patternCopy, -1);
                for (int k : subtreeLabels[j]) {
                    patternCopy[k + nrOfExternalNodes] = MLGenotypes.get(i)[k + nrOfExternalNodes];

                    if (k < nrOfExternalNodes)
                        patternCopy[k] = MLGenotypes.get(i)[k];
                }

                if (!ElementComparator.compIntArrays(uniqueSubtreePatterns[j], patternCopy))
                    uniqueSubtreePatterns[j].add(patternCopy);
            }
        }

        // compute the total combinations of unique root patterns and unique subtree patterns
        int comb = uniqueRootPatterns.size();
        for (int i = 0; i < childCount; i++) {
            assert uniqueSubtreePatterns[i] != null;
            comb *= uniqueSubtreePatterns[i].size();
        }

        // combine each unique root pattern to each unique subtree pattern
        for (int[] rootP : uniqueRootPatterns) {

            // store the current processed unique subtree patterns
            int[] marker = new int[childCount];

            while (true) {

                boolean termination = false;

                // check and update marker
                for (int i = 0; i < childCount; i++) {
                    if (marker[i] == uniqueSubtreePatterns[i].size()) {
                        marker[i] = 0;

                        if (i < childCount - 1)
                            marker[i + 1]++;
                        else
                            termination = true;
                    }
                }

                if (termination)
                    break;

                // copy root pattern
                int[] rootPCopy = new int[nrOfNodes + nrOfExternalNodes];
                System.arraycopy(rootP, 0, rootPCopy, 0, nrOfNodes + nrOfExternalNodes);

                // copy subtree patterns
                for (int i = 0; i < childCount; i++)
                    for (int k = 0; k < nrOfNodes + nrOfExternalNodes; k++)
                        rootPCopy[k] = rootPCopy[k] >= 0 ? rootPCopy[k] : uniqueSubtreePatterns[i].get(marker[i])[k];

                // whether the current combination already exists or not
                if (!ElementComparator.compIntArrays(MLGenotypes, rootPCopy))
                    MLGenotypes.add(rootPCopy);

                // update marker
                marker[0]++;

            }
        }

        assert (MLGenotypes.size() - restSize) == comb;

    } // updateMaxLikelihoodGenotypes

    /**
     * traverse subtree rooted at a node, and store node numbers in pre-order
     *
     * @param node   subtree root
     * @param labels where to store node numbers
     */
    private void traverseSubtree(final Node node, List<Integer> labels) {
        if (labels == null)
            return;

        if (node != null) {
            labels.add(node.getNr());

            if (!node.isLeaf())
                for (int i = 0; i < node.getChildCount(); i++)
                    traverseSubtree(node.getChild(i), labels);
        }
    } // traverseSubtree

    /**
     * get the elements to be deleted and added
     *
     * @param oldList previous set to be replaced
     * @param deleted elements in oldSet to be deleted
     * @param newList newly gotten set
     * @param added   elements in newSet to be added
     */
    private void getDeletedAndAddedGenotypes(final List<Integer> oldList, List<Integer> deleted,
                                             final List<Integer> newList, List<Integer> added) {
        for (int i : oldList)
            if (!newList.contains(i) && !deleted.contains(i))
                deleted.add(i);

        for (int i : newList)
            if (!oldList.contains(i) && !added.contains(i))
                added.add(i);
    } // getDeletedAndAddedGenotypes

    /**
     * update maximum likelihood genotypes and collection based on elements to be deleted and added for a node
     *
     * @param collection MLGenotypesCollection
     * @param genotypes  maxLikelihoodGenotypes
     * @param nodeIndex  which node to update?
     * @param deleted    elements to be deleted
     * @param added      elements to be added
     */
    private void updateCollectionAndGenotypes(
            List<Integer> collection,
            List<int[]> genotypes,
            final int nodeIndex,
            final List<Integer> deleted,
            final List<Integer> added
    ) {

        Iterator<Integer> aItr = added.iterator();

        // update genotypes
        if (deleted.size() == added.size() && deleted.size() != 0) {
            // if deleted and added are of the same size

            for (int oldElement : deleted) {
                //replace

                final int newElement = aItr.next();
                List<Integer> indices = getMLGenotypeIndices(genotypes, nodeIndex, oldElement);

                for (int i : indices)
                    genotypes.get(i)[nodeIndex] = newElement;

                // update collection
                collection.removeIf(e -> e == oldElement);
                if (!collection.contains(newElement))
                    collection.add(newElement);
            }
        } else if (deleted.size() > added.size()) {
            // if more elements are deleted

            for (int oldElement : deleted) {
                if (aItr.hasNext()) {
                    //replace

                    final int newElement = aItr.next();
                    List<Integer> indices = getMLGenotypeIndices(genotypes, nodeIndex, oldElement);

                    for (int i : indices)
                        genotypes.get(i)[nodeIndex] = newElement;

                    // update collection
                    collection.removeIf(e -> e == oldElement);
                    if (!collection.contains(newElement))
                        collection.add(newElement);
                } else {
                    // delete

                    genotypes.removeIf(e -> e[nodeIndex] == oldElement);

                    // update collection
                    collection.removeIf(e -> e == oldElement);
                }
            }
        } else if (deleted.size() < added.size()) {
            // if more elements are added

            List<int[]> targets = new ArrayList<>();
            boolean needCopy = true;

            for (int oldElement : deleted) {
                //replace

                final int newElement = aItr.next();
                List<Integer> indices = getMLGenotypeIndices(genotypes, nodeIndex, oldElement);

                for (int i : indices) {
                    genotypes.get(i)[nodeIndex] = newElement;

                    if (needCopy) {
                        targets.add(new int[genotypes.get(i).length]);
                        System.arraycopy(genotypes.get(i), 0, targets.get(targets.size() - 1), 0, genotypes.get(i).length);
                    }
                }

                needCopy = false;

                // update collection
                collection.removeIf(e -> e == oldElement);
                if (!collection.contains(newElement))
                    collection.add(newElement);
            }

            // add
            while (aItr.hasNext()) {
                final int newElement = aItr.next();

                // update collection
                if (!collection.contains(newElement))
                    collection.add(newElement);

                if (needCopy) {
                    needCopy = false;

                    List<Integer> entries = getMLGenotypeIndices(genotypes, nodeIndex, -2);

                    List<Integer> indices;
                    if (entries.contains(-1)) {
                        // replace

                        indices = getMLGenotypeIndices(genotypes, nodeIndex, -1);

                        for (int i : indices) {
                            genotypes.get(i)[nodeIndex] = newElement;

                            targets.add(new int[genotypes.get(i).length]);
                            System.arraycopy(genotypes.get(i), 0, targets.get(targets.size() - 1), 0, genotypes.get(i).length);
                        }

                        continue;
                    } else {
                        // copy one of them

                        indices = getMLGenotypeIndices(genotypes, nodeIndex, entries.iterator().next());

                        for (int i : indices) {
                            targets.add(new int[genotypes.get(i).length]);
                            System.arraycopy(genotypes.get(i), 0, targets.get(targets.size() - 1), 0, genotypes.get(i).length);
                        }
                    }
                }

                for (int[] i : targets) {
                    genotypes.add(new int[i.length]);
                    i[nodeIndex] = newElement;
                    System.arraycopy(i, 0, genotypes.get(genotypes.size() - 1), 0, i.length);
                }
            }
        }
    } // updateCollectionAndGenotypes

    /**
     * update maximum likelihood genotypes and collection based on elements to be deleted and added for a node
     *
     * @param collection       MLGenotypesCollection
     * @param genotypes        maxLikelihoodGenotypes
     * @param nodeIndex1       which node to refer to? (usually the parent node)
     * @param matchEntry       which entry you are searching for?
     * @param nodeIndex2       which node to update? (usually a child node)
     * @param deleted          elements to be deleted
     * @param added            elements to be added
     */
    private void updateCollectionAndGenotypes(
            List<Integer> collection,
            List<int[]> genotypes,
            final int nodeIndex1,
            final int matchEntry,
            final int nodeIndex2,
            final List<Integer> deleted,
            final List<Integer> added
    ) {

        Iterator<Integer> aItr = added.iterator();

        // update genotypes
        if (deleted.size() == added.size() && deleted.size() != 0) {
            // if deleted and added are of the same size

            for (int oldElement : deleted) {
                //replace

                final int newElement = aItr.next();

                List<Integer> indices = getMLGenotypeIndices(
                        genotypes,
                        nodeIndex1,
                        matchEntry,
                        nodeIndex2,
                        oldElement
                );

                for (int i : indices)
                    genotypes.get(i)[nodeIndex2] = newElement;

                // update collection
                collection.removeIf(e -> e == oldElement);
                if (!collection.contains(newElement))
                    collection.add(newElement);
            }
        } else if (deleted.size() > added.size()) {
            // if more elements are deleted

            for (int oldElement : deleted) {
                if (aItr.hasNext()) {
                    //replace

                    final int newElement = aItr.next();
                    List<Integer> indices = getMLGenotypeIndices(genotypes, nodeIndex1, matchEntry, nodeIndex2, oldElement);

                    for (int i : indices)
                        genotypes.get(i)[nodeIndex2] = newElement;

                    // update collection
                    collection.removeIf(e -> e == oldElement);
                    if (!collection.contains(newElement))
                        collection.add(newElement);
                } else {
                    // delete

                    genotypes.removeIf(e -> (e[nodeIndex2] == oldElement && e[nodeIndex1] == matchEntry));

                    // update collection
                    collection.removeIf(e -> e == oldElement);
                }
            }
        } else if (deleted.size() < added.size()) {
            // if more elements are added
            List<int[]> targets = new ArrayList<>();
            boolean needCopy = true;

            for (int oldElement : deleted) {
                //replace

                final int newElement = aItr.next();
                List<Integer> indices = getMLGenotypeIndices(genotypes, nodeIndex1, matchEntry, nodeIndex2, oldElement);

                for (int i : indices) {
                    genotypes.get(i)[nodeIndex2] = newElement;

                    if (needCopy) {
                        targets.add(new int[genotypes.get(i).length]);
                        System.arraycopy(genotypes.get(i), 0, targets.get(targets.size() - 1), 0, genotypes.get(i).length);
                    }
                }

                needCopy = false;

                // update collection
                collection.removeIf(e -> e == oldElement);
                if (!collection.contains(newElement))
                    collection.add(newElement);
            }

            // add
            while (aItr.hasNext()) {
                final int newElement = aItr.next();

                // update collection
                if (!collection.contains(newElement))
                    collection.add(newElement);

                if (needCopy) {
                    needCopy = false;

                    List<Integer> entries = getMLGenotypeIndices(genotypes, nodeIndex1, matchEntry, nodeIndex2, -2);
                    assert entries.size() > 0;

                    List<Integer> indices;
                    if (entries.contains(-1)) {
                        // replace

                        indices = getMLGenotypeIndices(genotypes, nodeIndex1, matchEntry, nodeIndex2, -1);

                        for (int i : indices) {
                            genotypes.get(i)[nodeIndex2] = newElement;

                            targets.add(new int[genotypes.get(i).length]);
                            System.arraycopy(genotypes.get(i), 0, targets.get(targets.size() - 1), 0, genotypes.get(i).length);
                        }

                        continue;
                    } else {
                        // copy one of them

                        indices = getMLGenotypeIndices(genotypes, nodeIndex1, matchEntry, nodeIndex2, entries.iterator().next());

                        for (int i : indices) {
                            targets.add(new int[genotypes.get(i).length]);
                            System.arraycopy(genotypes.get(i), 0, targets.get(targets.size() - 1), 0, genotypes.get(i).length);
                        }
                    }
                }

                for (int[] i : targets) {
                    genotypes.add(new int[i.length]);
                    i[nodeIndex2] = newElement;
                    System.arraycopy(i, 0, genotypes.get(genotypes.size() - 1), 0, i.length);
                }
            }
        }
    } // updateCollectionAndGenotypes

    /**
     * find indices in a list of integer arrays given the node index and match entry
     * if matchEntry == -2, then a set of all possible match entries for nodeIndex is returned
     *
     * @param arr a list of integer arrays
     * @param nodeIndex    index of interested node
     * @param matchEntry   which value you want to match?
     * @return a list of found indices
     */
    private List<Integer> getMLGenotypeIndices(
            final List<int[]> arr,
            final int nodeIndex,
            final int matchEntry
    ) {
        List<Integer> results = new ArrayList<>();

        for (int i = 0; i < arr.size(); i++) {
            if (matchEntry == -2 && !results.contains(arr.get(i)[nodeIndex]))
                results.add(arr.get(i)[nodeIndex]);
            else if (matchEntry >= -1) {
                if (arr.get(i)[nodeIndex] == matchEntry)
                    results.add(i);
            }
        }

        return results;
    } // getMLGenotypeIndex

    /**
     * find indices in a list of integer arrays given the node index and match entry
     * if matchEntry2 == -2, then a set of all possible match entries for nodeIndex2 is returned
     *
     * @param arr          a list of integer arrays
     * @param nodeIndex1   index to interested node 1
     * @param matchEntry1  which value you want to match for node 1?
     * @param nodeIndex2   index to interested node 2
     * @param matchEntry2  which value you want to match for node 2?
     * @return a list of found indices
     */
    private List<Integer> getMLGenotypeIndices(
            final List<int[]> arr,
            final int nodeIndex1,
            final int matchEntry1,
            final int nodeIndex2,
            final int matchEntry2
    ) {
        List<Integer> results = new ArrayList<>();

        for (int i = 0; i < arr.size(); i++) {
            if (matchEntry2 == -2) {
                if (arr.get(i)[nodeIndex1] == matchEntry1 && !results.contains(arr.get(i)[nodeIndex2]))
                    results.add(arr.get(i)[nodeIndex2]);
            } else if (matchEntry2 >= -1) {
                if (arr.get(i)[nodeIndex1] == matchEntry1 && arr.get(i)[nodeIndex2] == matchEntry2)
                    results.add(i);
            }
        }

        return results;
    } // getMLGenotypeIndex

    /**
     * get all genotypes of nodeIndex2 given that the genotype of nodeIndex1 (usually the parent of nodeIndex2) is matchEntry
     *
     * @param arr a list of integer arrays
     * @param nodeIndex1   node to match
     * @param matchEntry   genotype of node 1
     * @param nodeIndex2   node to search
     * @return a list of found genotypes
     */
    private List<Integer> getMLGenotypeEntries(
            final List<int[]> arr,
            final int nodeIndex1,
            final int matchEntry,
            final int nodeIndex2
    ) {
        List<Integer> results = new ArrayList<>();

        for (int[] i : arr)
            if (i[nodeIndex1] == matchEntry && i[nodeIndex2] >= 0 && !results.contains(i[nodeIndex2]))
                results.add(i[nodeIndex2]);

        return results;
    } // getMLGenotypeEntries

    /**
     * for variant calling
     * get maximum likelihood genotypes for all cells and all sites
     */
    public void collectMLGenotypesAndLogLikelihoods() {
        final List<Node> treeTips = treeInput.get().getExternalNodes();
        final Node[] treeNodes = treeInput.get().getNodesAsArray();

        int duplicateLines = 0;

        for (int patternIndex = 0; patternIndex < nrOfPatterns; patternIndex++) {

            // initialise variantsInfo member variables at current pattern
            this.variantsInfo.initialiseAPattern(patternIndex);

            // maximum likelihood categories for current pattern
            List<Integer> mlCategory = this.variantsInfo.getMLCategory(patternIndex);

            // candidate maximum likelihood genotypes and ADO states
            List<int[]> candidateGenotypes = new ArrayList<>();

            // collect maximum likelihood genotypes under maximum likelihood categories
            for (int matrixIndex : mlCategory)
                candidateGenotypes.addAll(this.maxLikelihoodGenotypes[matrixIndex * this.nrOfPatterns + patternIndex]);

            if (candidateGenotypes.size() == 1) {
                // hooray! only one combination of maximum likelihood genotypes found

                final int[] maxLikelihoodGenotypes = Arrays.copyOfRange(
                        candidateGenotypes.get(0),
                        nrOfExternalNodes,
                        candidateGenotypes.get(0).length
                );
                final int[] maxLikelihoodAdos = Arrays.copyOfRange(
                        candidateGenotypes.get(0),
                        0,
                        nrOfExternalNodes
                );

                final int[] mlGenotypesTips = getMLStatesForTips(treeTips, maxLikelihoodGenotypes);
                final int[] mlAdo = getMLStatesForTips(treeTips, maxLikelihoodAdos);
                double[] genotypeLogLikelihoodsTips = new double[this.nrOfExternalNodes];
                double[] genotypeLogLikelihoodsNodes = new double[this.nrOfNodes];

                getMLGenotypesLogLikelihood(
                        patternIndex,
                        mlCategory.get(0),
                        maxLikelihoodGenotypes,
                        treeNodes,
                        genotypeLogLikelihoodsTips,
                        genotypeLogLikelihoodsNodes
                );

                this.variantsInfo.addGenotypeLogLikelihoodConstantPattern(
                        patternIndex,
                        getLogLikelihoodConstantPattern(
                                patternIndex,
                                mlCategory
                        ),
                        true
                );
                this.variantsInfo.addMLGenotypesNodes(
                        patternIndex,
                        maxLikelihoodGenotypes,
                        true
                );
                this.variantsInfo.addMLGenotypesTips(patternIndex, mlGenotypesTips, true);
                this.variantsInfo.addMLAdo(patternIndex, mlAdo, true);
                this.variantsInfo.addGenotypeLogLikelihoodsNodes(patternIndex, genotypeLogLikelihoodsNodes, true);
                this.variantsInfo.addGenotypeLogLikelihoodsTips(patternIndex, genotypeLogLikelihoodsTips, true);

                addLogLikelihoodsForAllGenotypes(treeTips, patternIndex, mlCategory.get(0), true);

                addLogTransferProbabilitiesForAllGenotypes(treeTips, mlCategory.get(0), maxLikelihoodGenotypes, true);

                addAdoLogLikelihoods(treeTips, patternIndex, mlCategory.get(0), maxLikelihoodGenotypes, true);

            } else if (candidateGenotypes.size() > 1) {
                // oops! multiple combinations of maximum likelihood genotypes found (which is rare)
                // choose those with the largest sum of log likelihood of all nodes

                // combinations of genotypes which have the maximum sum of log likelihood
                candidateGenotypes.clear();

                List<int[]> maxLikelihoodGenotypesAll = new ArrayList<>();
                List<int[]> mlGenotypesTipsAll = new ArrayList<>();
                List<int[]> mlAdoAll = new ArrayList<>();
                List<double[]> genotypeLogLikelihoodsTipsAll = new ArrayList<>();
                List<double[]> genotypeLogLikelihoodsNodesAll = new ArrayList<>();
                List<Integer> categories = new ArrayList<>();

                double max = -Double.MAX_VALUE;

                for (int category : mlCategory) {

                    List<int[]> genotypesAll = this.maxLikelihoodGenotypes[category * this.nrOfPatterns + patternIndex];

                    for (int[] genotypes : genotypesAll) {

                        final int[] maxLikelihoodGenotypes = Arrays.copyOfRange(
                                genotypes,
                                nrOfExternalNodes,
                                genotypes.length
                        );
                        final int[] maxLikelihoodAdos = Arrays.copyOfRange(
                                genotypes,
                                0,
                                nrOfExternalNodes
                        );

                        final int[] mlGenotypesTips = getMLStatesForTips(treeTips, maxLikelihoodGenotypes);
                        final int[] mlAdo = getMLStatesForTips(treeTips, maxLikelihoodAdos);
                        double[] genotypeLogLikelihoodsTips = new double[this.nrOfExternalNodes];
                        double[] genotypeLogLikelihoodsNodes = new double[this.nrOfNodes];

                        final double logValue = getMLGenotypesLogLikelihood(
                                patternIndex,
                                category,
                                maxLikelihoodGenotypes,
                                treeNodes,
                                genotypeLogLikelihoodsTips,
                                genotypeLogLikelihoodsNodes
                        );

                        if (max < logValue) {

                            max = logValue;

                            candidateGenotypes.clear();
                            candidateGenotypes.add(genotypes);

                            maxLikelihoodGenotypesAll.clear();
                            maxLikelihoodGenotypesAll.add(maxLikelihoodGenotypes);

                            mlGenotypesTipsAll.clear();
                            mlGenotypesTipsAll.add(mlGenotypesTips);

                            mlAdoAll.clear();
                            mlAdoAll.add(mlAdo);

                            genotypeLogLikelihoodsTipsAll.clear();
                            genotypeLogLikelihoodsTipsAll.add(genotypeLogLikelihoodsTips);

                            genotypeLogLikelihoodsNodesAll.clear();
                            genotypeLogLikelihoodsNodesAll.add(genotypeLogLikelihoodsNodes);

                            categories.clear();
                            categories.add(category);

                        } else if (max == logValue) {

                            candidateGenotypes.add(genotypes);

                            maxLikelihoodGenotypesAll.add(maxLikelihoodGenotypes);

                            mlGenotypesTipsAll.add(mlGenotypesTips);

                            mlAdoAll.add(mlAdo);

                            genotypeLogLikelihoodsTipsAll.add(genotypeLogLikelihoodsTips);

                            genotypeLogLikelihoodsNodesAll.add(genotypeLogLikelihoodsNodes);

                            categories.add(category);

                        }

                    }

                }

                if (candidateGenotypes.size() == 1) {
                    // hooray! only one combination of maximum likelihood genotypes with the largest log likelihood found

                    this.variantsInfo.addMLCategory(patternIndex, categories.get(0), true);

                    this.variantsInfo.addGenotypeLogLikelihoodConstantPattern(
                            patternIndex,
                            getLogLikelihoodConstantPattern(
                                    patternIndex,
                                    categories
                            ),
                            true
                    );
                    this.variantsInfo.addMLGenotypesNodes(
                            patternIndex,
                            maxLikelihoodGenotypesAll.get(0),
                            true
                    );
                    this.variantsInfo.addMLGenotypesTips(patternIndex, mlGenotypesTipsAll.get(0), true);
                    this.variantsInfo.addMLAdo(patternIndex, mlAdoAll.get(0), true);
                    this.variantsInfo.addGenotypeLogLikelihoodsNodes(patternIndex, genotypeLogLikelihoodsNodesAll.get(0), true);
                    this.variantsInfo.addGenotypeLogLikelihoodsTips(patternIndex, genotypeLogLikelihoodsTipsAll.get(0), true);

                    addLogLikelihoodsForAllGenotypes(treeTips, patternIndex, categories.get(0), true);

                    addLogTransferProbabilitiesForAllGenotypes(treeTips, categories.get(0), maxLikelihoodGenotypesAll.get(0), true);

                    addAdoLogLikelihoods(treeTips, patternIndex, categories.get(0), maxLikelihoodGenotypesAll.get(0), true);

                } else if (candidateGenotypes.size() > 1) {
                    // oops! multiple combinations of maximum likelihood genotypes with the largest log likelihood found (which is much rarer)

                    // note that 'categories' and 'candidateGenotypes' share the same size
                    for (int i = 0; i < categories.size(); i++) {

                        addLogLikelihoodsForAllGenotypes(
                                treeTips,
                                patternIndex,
                                categories.get(i),
                                i == 0
                        );

                        addLogTransferProbabilitiesForAllGenotypes(
                                treeTips,
                                categories.get(i),
                                maxLikelihoodGenotypesAll.get(i),
                                i == 0
                        );


                        addAdoLogLikelihoods(
                                treeTips,
                                patternIndex,
                                categories.get(i),
                                maxLikelihoodGenotypesAll.get(i),
                                i == 0
                        );

                    }

                    this.variantsInfo.addMLCategory(patternIndex, categories, true);

                    this.variantsInfo.addGenotypeLogLikelihoodConstantPattern(
                            patternIndex,
                            getLogLikelihoodConstantPattern(
                                    patternIndex,
                                    categories
                            ),
                            true
                    );
                    this.variantsInfo.addMLGenotypesNodes(
                            patternIndex,
                            maxLikelihoodGenotypesAll,
                            true
                    );
                    this.variantsInfo.addMLGenotypesTips(patternIndex, mlGenotypesTipsAll, true);
                    this.variantsInfo.addMLAdo(patternIndex, mlAdoAll, true);
                    this.variantsInfo.addGenotypeLogLikelihoodsNodes(patternIndex, genotypeLogLikelihoodsNodesAll, true);
                    this.variantsInfo.addGenotypeLogLikelihoodsTips(patternIndex, genotypeLogLikelihoodsTipsAll, true);

                    duplicateLines++;
                    System.err.print("More than one variant calling results for " + duplicateLines + " patterns. " +
                            "There will be lines in the final variant calling document with duplicate locus names " +
                            "but different variant calling results for some cells. Which one to be used is up to " +
                            "users.\r");


                } else
                    throw new NullPointerException("Error: no genotypes found.");

            } else
                throw new NullPointerException("Error: no genotypes found.");

        }
    } // collectMLGenotypesAndLogLikelihoods

    /**
     * get log likelihood of pattern if being constant (matrix-wise).
     *
     * @param patternIndex  apparently
     * @param matrixIndices apparently
     * @return a list of log likelihoods
     */
    public List<Double> getLogLikelihoodConstantPattern(
            int patternIndex,
            final List<Integer> matrixIndices
    ) {
        List<Double> results = new ArrayList<>();

        for (int matrixIndex : matrixIndices) {
            results.add(
                    ((ScsBeerLikelihoodCore) this.likelihoodCore).getLogLikelihoodConstantPattern(
                            this.treeInput.get().getRoot(),
                            matrixIndex,
                            patternIndex,
                            this.substitutionModel.getRootGenotype()
                    )
            );
        }

        return results;
    } // getLogLikelihoodConstantPattern

    /**
     * get log likelihoods for all genotypes
     *
     * @param tips         apparently
     * @param patternIndex apparently
     * @param matrixIndex  apparently
     * @param overwrite    apparently
     */
    private void addLogLikelihoodsForAllGenotypes(
            final List<Node> tips,
            int patternIndex,
            int matrixIndex,
            boolean overwrite
    ) {
        for (Node i : tips) {
            this.variantsInfo.addGenotypeLogLikelihoodsAll(
                    patternIndex,
                    i.getNr(),
                    ((ScsBeerLikelihoodCore) this.likelihoodCore).getLogLikelihoodsAllStates(
                            i,
                            matrixIndex,
                            patternIndex
                    ),
                    overwrite
            );
        }
    } // addLogLikelihoodsForAllGenotypes

    private void addLogTransferProbabilitiesForAllGenotypes(
            final List<Node> tips,
            int matrixIndex,
            final int[] mlGenotypesNode,
            boolean overwrite
    ) {
        for (Node i : tips) {
            final int parentGenotype = mlGenotypesNode[i.getParent().getNr()];

            this.variantsInfo.addGenotypeLogTransferProbabilitiesAll(
                    i.getNr(),
                    parentGenotype,
                    this.nrOfStates,
                    ((ScsBeerLikelihoodCore) this.likelihoodCore).getLogTransferProbabilityAll(
                            i,
                            matrixIndex,
                            parentGenotype
                    ),
                    overwrite
            );
        }
    } // addLogTransferProbabilitiesForAllGenotypes

    private void addAdoLogLikelihoods(
            final List<Node> tips,
            int patternIndex,
            int matrixIndex,
            final int[] mlGenotypesNode,
            boolean overwrite
    ) {
        for (Node i : tips) {
            this.variantsInfo.addAdoLogLikelihoods(
                    patternIndex,
                    i.getNr(),
                    nrOfAdoStates,
                    rawReadCountsModel.getAdoLogLikelihoods(
                            matrixIndex,
                            patternIndex,
                            i.getNr(),
                            mlGenotypesNode[i.getNr()]
                    ),
                    overwrite
            );
        }
    } // addAdoLogLikelihoods

    /**
     * get the log likelihood for a group of specific maximum likelihood genotypes / ADO states
     *
     * @param patternIndex                which pattern?
     * @param matrixIndex                 which matrix?
     * @param mlGenotypes                 maximum likelihood genotypes
     * @param nodes                       tree nodes (internal and external)
     * @param genotypeLogLikelihoodsTips  genotype log likelihoods of tips
     * @param genotypeLogLikelihoodsNodes genotype log likelihoods of nodes
     * @return log likelihood
     */
    public double getMLGenotypesLogLikelihood(
            int patternIndex,
            int matrixIndex,
            int[] mlGenotypes,
            final Node[] nodes,
            double[] genotypeLogLikelihoodsTips,
            double[] genotypeLogLikelihoodsNodes
    ) {
        for (Node i : nodes) {
            final int index = i.getNr();

            genotypeLogLikelihoodsNodes[index] = ((ScsBeerLikelihoodCore) likelihoodCore).getLogLikelihood(
                    i,
                    matrixIndex,
                    patternIndex,
                    mlGenotypes[index]
            );

            if (i.isLeaf()) {
                if (index < genotypeLogLikelihoodsTips.length)
                    genotypeLogLikelihoodsTips[index] = genotypeLogLikelihoodsNodes[index];
                else
                    throw new IllegalArgumentException("Leaf index out of bound. Indices of tree tips should be rearranged.");
            }
        }

        return MathFunctions.logSumExp(genotypeLogLikelihoodsNodes);
    } // getMLGenotypesLogLikelihood

    /**
     * get maximum likelihood states (either genotypes or ADOs) for tips
     *
     * @param tips apparently
     * @param all  genotypes for all nodes
     * @return genotypes for tips
     */
    public int[] getMLStatesForTips(
            final List<Node> tips,
            final int[] all
    ) {
        int[] result = new int[tips.size()];

        for (Node i : tips) {
            final int index = i.getNr();

            if (index >= tips.size())
                throw new IllegalArgumentException("Indices of tips are not within [0," + (tips.size() - 1) + "]. " +
                        "Labels of tips are not updated.");

            result[index] = all[index];
        }

        return result;
    } // getMLStatesForTips

    @Override
    public void callVariants() {
        if (!MLGenotypesAndAdosUpdated)
            getMLGenotypes(treeInput.get().getRoot());

        collectMLGenotypesAndLogLikelihoods();
    } // callVariants

    /**
     * perform post process
     */
    private void postProcess() {
        // 1. compute the maximum likelihood genotype and number of sequenced alleles for each node on the tree
        final long startTime1 = System.currentTimeMillis();
        getMLGenotypes(treeInput.get().getRoot());
        final long endTime1 = System.currentTimeMillis();
        if (runTimeAnalysis) {
            runTime[6] = endTime1 - startTime1;

            System.out.println("1. Compute maximum likelihood genotypes time: " + runTime[6] + " milliseconds.");
            System.out.println("     (Pre)  Changed matrices and patterns: " + changedPatternsList.size() + " out of " + nrOfMatrices * nrOfPatterns);
        }

        // 2.0 update allelic sequencing coverage and raw variance and remove some unchanged patterns
        final long startTime2 = System.currentTimeMillis();
        rawReadCountsModel.updateAllelicSeqCovAndRawVar(changedPatternsList);
        final long endTime2 = System.currentTimeMillis();
        if (runTimeAnalysis) {
            runTime[7] = endTime2 - startTime2;

            System.out.println("2.0 Update raw read counts model time: " + runTime[7] + " milliseconds.");
            System.out.println("     (Post) Changed matrices and patterns: " + changedPatternsList.size() + " out of " + nrOfMatrices * nrOfPatterns);
        }

        // proceed only if there are matrices and patterns changed
        if (changedPatternsList.size() > 0) {

            // 2.1 store changed patterns into an integer array to gain efficiency
            final long startTime3 = System.currentTimeMillis();
            changedPatterns = new int[changedPatternsList.size()][2];
            for (int i = 0; i < changedPatternsList.size(); i++)
                changedPatterns[i] = Ints.toArray(changedPatternsList.get(i));
            final long endTime3 = System.currentTimeMillis();
            if (runTimeAnalysis)
                System.out.println("2.1 Conversion from list to array time: " + (endTime3 - startTime3) + " milliseconds.");

            // 2.2 sort changedPatterns according to the value of pattern
            final long startTime4 = System.currentTimeMillis();
            updateChangedPatterns(null);
            final long endTime4 = System.currentTimeMillis();
            if (runTimeAnalysis)
                System.out.println("2.2 Process changed patterns time: " + (endTime4 - startTime4) + " milliseconds.");

            // 3. traverse the entire tree for changed patterns
            hasDirt = Tree.IS_DIRTY;
            updateLeaves = true;
            final long startTime5 = System.currentTimeMillis();
            traverse(treeInput.get().getRoot(), true);
            final long endTime5 = System.currentTimeMillis();
            if (runTimeAnalysis) {
                runTime[8] = endTime5 - startTime5;

                System.out.println("3. Partially traverse time: " + runTime[8] + " milliseconds.");
                System.out.println("   Leaf partial assignment time: " + runTime[9] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[9] / runTime[8])));
                System.out.println("   Likelihood core run time: " + runTime[10] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[10] / runTime[8])));
                System.out.println("   Integrate across matrices time: " + runTime[11] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[11] / runTime[8])));
                System.out.println("   Log likelihood of each pattern time: " + runTime[12] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[12] / runTime[8])));
                final long tmp = runTime[8] - runTime[9] - runTime[10] - runTime[11] - runTime[12];
                System.out.println("   Time for other things: " + tmp + " milliseconds, proportion: " + String.format("%.5f", ((double) tmp / runTime[8])) + "\n");
            }

        } else {
            if (runTimeAnalysis)
                System.out.println("No changed patterns detected.");
        }

    } // postProcess

    /**
     * should only be called during post process when scaling is turned on
     * <p>
     * if addedChangedPatterns is not null or empty, expand current changedPatterns.
     *
     * @param addedChangedPatterns potentially added changed patterns; allowed to be null
     */
    private void updateChangedPatterns(final List<int[]> addedChangedPatterns) {
        if (addedChangedPatterns != null && addedChangedPatterns.size() > 0) {
            final int newLen = changedPatterns.length + addedChangedPatterns.size();
            int[][] newChangedPatterns = new int[newLen][2];
            System.arraycopy(changedPatterns, 0, newChangedPatterns, 0, changedPatterns.length);

            for (int i = changedPatterns.length; i < newLen; i++)
                System.arraycopy(addedChangedPatterns.get(i - changedPatterns.length), 0, newChangedPatterns[i], 0, 2);

            changedPatterns = newChangedPatterns;

            if (runTimeAnalysis)
                System.out.println("     Changed matrices and patterns updated: " + newLen + " out of " + nrOfMatrices * nrOfPatterns);
        }

        PatternComparator comparator = new PatternComparator();
        Arrays.sort(changedPatterns, comparator);
        changedPatternsIndex.clear();
        changedPatternsIndex.add(0);
        for (int i = 1; i < changedPatterns.length; i++)
            if (comparator.compare(changedPatterns[i - 1], changedPatterns[i]) != 0)
                changedPatternsIndex.add(i);
    } // updateChangedPatterns

    double m_fScale = 1.01;

    /**
     * compute the log likelihood of the current state
     * including background information if it is available
     *
     * @return log tree likelihood
     */
    @Override
    public double calculateLogP() {
        /*
        if (beagle != null) {
            logP = beagle.calculateLogP();
            return logP;
        }
         */

        final TreeInterface tree = treeInput.get();

        try {

            if (runTimeAnalysis) {
                runTime = null;
                runTime = new long[13];
            }

            final long startTime = System.currentTimeMillis();
            if (traverse(tree.getRoot(), false) != Tree.IS_CLEAN) {
                final long endTime = System.currentTimeMillis();
                if (runTimeAnalysis) {
                    runTime[0] = (endTime - startTime);

                    // print some run time information
                    System.out.println("Tree likelihood run time: " + runTime[0] + " milliseconds.");

                    System.out.println("  Transition probability matrix update time: " + runTime[1] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[1] / runTime[0])));
                    System.out.println("  Leaf partial assignment time: " + runTime[2] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[2] / runTime[0])));
                    System.out.println("  Likelihood core run time: " + runTime[3] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[3] / runTime[0])));
                    System.out.println("  Integrate across matrices time: " + runTime[4] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[4] / runTime[0])));
                    System.out.println("  Log likelihood of each pattern time: " + runTime[5] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[5] / runTime[0])));

                    final long tmp = runTime[0] - runTime[1] - runTime[2] - runTime[3] - runTime[4] - runTime[5];
                    System.out.println("  Time for other things: " + tmp + " milliseconds, proportion: " + String.format("%.5f", ((double) tmp / runTime[0])) + "\n");
                }
                calcLogP(false);
            }
        } catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }

        if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 &&
                !scaling.get().equals(ScsTreeLikelihood.Scaling.none)) {
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            unstore();
            hasDirt = Tree.IS_FILTHY;
            updateLeaves = true;
            traverse(tree.getRoot(), false);
            calcLogP(false);
            return logP;
        }

        return logP;
    } // calculateLogP

    /**
     * @param returnConstSum whether to return the sum of the likelihoods of constant sites
     */
    private void calcLogP(boolean returnConstSum) {
        logP = 0.0;

        List<Double> logConstRootAll = new ArrayList<>();
        for (int i = 0; i < scsDataInput.get().getPatternCount(); i++) {
            final int weight = scsDataInput.get().getPatternWeight(i);

            logP += patternLogLikelihoods[i] * weight;
            logConstRootAll.addAll(Collections.nCopies(weight, logConstRoot[i]));
        }

        if (useAscBiasCorrection) {
            biasCorr = scsDataInput.get().getAscBiasCorrection(
                    logConstRootAll.stream().mapToDouble(Double::doubleValue).toArray(),
                    returnConstSum
            );

            if (!scsDataInput.get().isReturnConstSum()) {
                logP += biasCorr;
                biasCorr = 0.0;

                if (runTimeAnalysis) {
                    System.out.println("Ascertainment bias correction: ");
                    System.out.println("  Uncorrected tree likelihood: " + (logP - biasCorr));
                    System.out.println("  Corrected tree likelihood: " + logP);
                    System.out.println("  Correction term: " + biasCorr);
                }
            }
        }

    } // calcLogP

    public double calculateLogPByThread() {
        /*
        if (beagle != null) {
            logP = beagle.calculateLogP();
            return logP;
        }
         */

        final TreeInterface tree = treeInput.get();

        try {

            if (runTimeAnalysis) {
                runTime = null;
                runTime = new long[13];
            }

            final long startTime = System.currentTimeMillis();
            if (traverse(tree.getRoot(), false) != Tree.IS_CLEAN) {
                final long endTime = System.currentTimeMillis();
                if (runTimeAnalysis) {
                    runTime[0] = (endTime - startTime);

                    // print some run time information
                    System.out.println("Tree likelihood run time: " + runTime[0] + " milliseconds.");

                    System.out.println("  Transition probability matrix update time: " + runTime[1] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[1] / runTime[0])));
                    System.out.println("  Leaf partial assignment time: " + runTime[2] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[2] / runTime[0])));
                    System.out.println("  Likelihood core run time: " + runTime[3] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[3] / runTime[0])));
                    System.out.println("  Integrate across matrices time: " + runTime[4] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[4] / runTime[0])));
                    System.out.println("  Log likelihood of each pattern time: " + runTime[5] + " milliseconds, proportion: " + String.format("%.5f", ((double) runTime[5] / runTime[0])));

                    final long tmp = runTime[0] - runTime[1] - runTime[2] - runTime[3] - runTime[4] - runTime[5];
                    System.out.println("  Time for other things: " + tmp + " milliseconds, proportion: " + String.format("%.5f", ((double) tmp / runTime[0])) + "\n");
                }
                calcLogP(true);
            }
        } catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }

        if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 &&
                !scaling.get().equals(ScsTreeLikelihood.Scaling.none)) {
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            unstore();
            hasDirt = Tree.IS_FILTHY;
            updateLeaves = true;
            traverse(tree.getRoot(), false);
            calcLogP(true);
            return logP;
        }

        return logP;
    } // calculateLogP

    /**
     * reset some variables for getMLGenotypes
     */
    private void resetVariables() {
        changedPatternsList = null;
        changedPatternsList = new ArrayList<>();

        changedPatternsSet = null;
        changedPatternsSet = new HashSet<>();

        changedPatterns = null;

        changedPatternsIndex = null;
        changedPatternsIndex = new ArrayList<>();

        changedPatternsPerNode = null;
        changedPatternsPerNode = new HashSet[nrOfNodes];
        for (int i = 0; i < nrOfNodes; i++)
            changedPatternsPerNode[i] = new HashSet<>();
    } // resetVariables


    //************************************************
    //*                Nested classes                *
    //************************************************

    static class ListComparator implements Comparator<int[]> {

        @Override
        public int compare(int[] o1, int[] o2) {
            if (o1.length != o2.length) {
                throw new IllegalArgumentException("Integer arrays to be sorted should be of the same length!");
            }
            for (int i = 0; i < o1.length; i++) {
                if (o1[i] > o2[i]) {
                    return 1;
                }
                if (o1[i] < o2[i]) {
                    return -1;
                }
            }

            return 0;
        } // compare

    } // class ListComparator

    // sort (matrix, pattern) according to pattern
    static class PatternComparator implements Comparator<int[]> {

        @Override
        public int compare(int[] o1, int[] o2) {
            if (o1.length != 2 || o2.length != 2) {
                throw new IllegalArgumentException("Expect arrays in a form of (matrix, pattern)!");
            }

            if (o1[1] > o2[1]) {
                return 1;
            } else if (o1[1] < o2[1]) {
                return -1;
            } else {
                return 0;
            }
        } // compare

    } // class ListComparator


    //***********************************************
    //*          Calculation nodes methods          *
    //***********************************************

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        /*
        // todo: beagle likelihood
        if (beagle != null) {
            return beagle.requiresRecalculation();
        }
         */

        hasDirt = Tree.IS_CLEAN;
        updateLeaves = false;

        // raw read counts model should be checked first
        if (rawReadCountsModel.isDirtyCalculation() || (rawReadCountsModel.isInDebugMode() && rawReadCountsModel.updateSeqCovModel())) {
            hasDirt = Tree.IS_DIRTY;
            updateLeaves = true;
            return true;
        }

        if (scsDataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }

        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }

        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            if (isForcingTreeDirtyRMC) {
                hasDirt = Tree.IS_DIRTY;
            }

            return true;
        }

        return treeInput.get().somethingIsDirty();
    } // requiresRecalculation

    @Override
    public void store() {
        /*
        if (beagle != null) {
            beagle.store();
            super.store();
            return;
        }
         */

        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);

        if (traceMLGenotypes)
            System.arraycopy(currentMLGenotypesNodeIndex, 0, storedMLGenotypesNodeIndex, 0, currentMLGenotypesNodeIndex.length);
    } // store

    public void unstore() {
        if (likelihoodCore != null) {
            likelihoodCore.unstore();
        }

        System.arraycopy(storedBranchLengths, 0, m_branchLengths, 0, storedBranchLengths.length);

        if (traceMLGenotypes)
            System.arraycopy(storedMLGenotypesNodeIndex, 0, currentMLGenotypesNodeIndex, 0, storedMLGenotypesNodeIndex.length);
    }

    @Override
    public void restore() {
        /*
        if (beagle != null) {
            beagle.restore();
            super.restore();
            return;
        }
         */

        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp1 = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp1;

        if (traceMLGenotypes) {
            int[] tmp2 = currentMLGenotypesNodeIndex;
            currentMLGenotypesNodeIndex = storedMLGenotypesNodeIndex;
            storedMLGenotypesNodeIndex = tmp2;
        }
    } // restore

    /**
     * the proposal is accepted
     * <p>
     * some other work should be processed at this stage in order to boost efficiency, including:
     * 1. computing maximum likelihood genotypes for all matrices and all patterns; meanwhile, finding out which
     * patterns for which matrices have been changed.
     * 2. for those changed patterns and matrices, updating allelic sequencing coverage and raw variance in raw read counts model.
     * 3. for those changed patterns and matrices, traversing the entire tree and recomputing each node's partial
     * likelihoods and maximum likelihood genotypes.
     */
    @Override
    public void accept() {
        // if not in debug mode, perform post processing
        if (!this.rawReadCountsModel.isInDebugMode()) {

            if (runTimeAnalysis) {
                System.out.print("Move is accepted. ");

                if (this.rawReadCountsModel.updateSeqCovModel())
                    System.out.println("Performing post process...");
                else
                    System.out.println("Nothing else needs to be processed.");
            }

            if (this.rawReadCountsModel.updateSeqCovModel()) {

                // only store the most recent accepted allelic sequencing coverage and raw variance
                this.rawReadCountsModel.storeSeqCovInfo();

                postProcess();
            }
        }

        super.accept();
    } // accept


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    /**
     * Get the variantsInfo object. Only used by multithreading when maximum likelihood genotypes and ADO states
     * are being traced and logged during MCMC.
     *
     * @return reference to `variantsInfo`
     */
    public GenericVariantsInfo.Base getVariantsInfo() {
        return variantsInfo;
    } // getVariantsInfo

    /**
     * Is tracing maximum likelihood genotypes or not.
     *
     * @return yes or no
     */
    public boolean isTraceMLGenotypes() {
        return traceMLGenotypesInput.get();
    } // isTraceMLGenotypes

    private void setMLGenotypesNodeForUpdate(int nodeIndex) {
        currentMLGenotypesNodeIndex[nodeIndex] = 1 - currentMLGenotypesNodeIndex[nodeIndex];
    } // setMLGenotypesNodeForUpdate

    /**
     * Set the active index marker for node partials and, if applicable, maximum likelihood genotypes for update.
     *
     * @param nodeIndex apparently
     */
    private void setIndexForUpdate(int nodeIndex) {
        if (traceMLGenotypes)
            setMLGenotypesNodeForUpdate(nodeIndex);

        likelihoodCore.setNodePartialsForUpdate(nodeIndex);
    } // setForUpdate

    /* return copy of pattern log likelihoods for each of the patterns in the alignment */
    public double[] getPatternLogLikelihoods() {
        /*
        if (beagle != null) {
            return beagle.getPatternLogLikelihoods();
        }
         */
        return patternLogLikelihoods.clone();
    } // getPatternLogLikelihoods

    /**
     * for the purpose of debugging (to accurately compute tree likelihood)
     *
     * @param inDebugMode in debug mode or not
     */
    @Override
    public void updateInDebugMode(final boolean inDebugMode) {
        rawReadCountsModel.setInDebugMode(inDebugMode);
    } // updateInDebugMode

    public int[] getCurrentMLGenotypesNodeIndex() {
        return currentMLGenotypesNodeIndex;
    } // getCurrentMLGenotypesNodeIndex

    /**
     * Get the bias correction term, usually the sum of the constant sites.
     *
     * @return bias correction term
     */
    public double getBiasCorr() {
        return biasCorr;
    }

    public String getSortedCellNamesFromVarInfo(String separator) {
        return variantsInfo.getSortedTipNames(separator);
    } // getSortedCellNamesFromVarInfo

    /**
     * initialize header of logger
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    @Override
    public void initSiteHeader(int start, PrintStream out) {
        scsDataInput.get().logSiteHeader(start, out);
    } // initSiteHeader


    //******************************************
    //*           Logger methods for           *
    //*        genotypes and ADO states        *
    //******************************************

    /**
     * close
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    public void closeGenotypeAndAdo(int start, PrintStream out) {
        scsDataInput.get().logSiteMap(start, out);
    } // closeGenotypeAndAdo


    //*******************************************
    //*           Logger methods for            *
    //*    allelic coverage and raw variance    *
    //*******************************************

    /**
     * log sampled allelic sequencing coverage and raw variance for each site
     *
     * @param out apparently
     */
    @Override
    public void logCovar(PrintStream out) {
        for (int i = 0; i != scsDataInput.get().getLociNr(); i++) {
            rawReadCountsModelInput.get().logCovar(
                    scsDataInput.get().getPatternIndex(i),
                    out
            );

            if (i < scsDataInput.get().getLociNr() - 1)
                out.print("\t");
        }
    } // logCovar

    /**
     * close
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    @Override
    public void closeCovar(int start, PrintStream out) {
        scsDataInput.get().logSiteMap(start, out);
    } // closeCovar

}