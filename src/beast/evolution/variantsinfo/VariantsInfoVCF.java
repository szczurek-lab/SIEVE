package beast.evolution.variantsinfo;

import beast.evolution.alignment.VariantSiteInfo;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import beast.evolution.tree.Node;
import beast.evolution.tree.ScsTree;
import beast.evolution.variantsinfo.vcfentry.CandidateAltNuc;
import beast.evolution.variantsinfo.vcfentry.CellLocusGenotype;
import beast.evolution.variantsinfo.vcfentry.VCFEntry;
import beast.evolution.variantsinfo.vcfentry.VCFInfo;
import com.google.common.primitives.Ints;
import org.jetbrains.annotations.NotNull;

import java.io.PrintStream;
import java.util.*;

import static beast.math.util.MathFunctions.convertLogE2RoundedPhredScaled;

public class VariantsInfoVCF extends GenericVariantsInfoVCF {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    /**
     * for variant calling
     * which category each pattern belongs to?
     * [#patterns]
     */
    protected List<Integer>[] MLCategory;

    /**
     * for variant calling and all tips
     * which genotype of a specific cell at a specific pattern is inferred to be?
     * [#patterns]
     */
    protected List<int[]>[] MLGenotypesTips;

    /**
     * for variant calling and all nodes
     * which genotype of a specific node at a specific pattern is inferred to be?
     * [#patterns]
     */
    protected List<int[]>[] MLGenotypesNodes;

    /**
     * for variant calling and all tips
     * The maximum likelihood ADO state for each tip
     * [#patterns]
     */
    protected List<int[]>[] MLAdo;

    /**
     * for variant calling
     * the log likelihood of a variant pattern being constant (matrix-wise)
     * [#patterns]
     */
    protected List<Double>[] genotypeLogLikelihoodConstantPattern;

    /**
     * for variant calling and all tips
     * stores the maximum log likelihoods
     * [#patterns]
     */
    protected List<double[]>[] genotypeLogLikelihoodsTips;

    /**
     * for variant calling and all nodes
     * stores the maximum log likelihoods
     * [#patterns]
     */
    protected List<double[]>[] genotypeLogLikelihoodsNodes;

    /**
     * for variant calling and all tips
     * stores the log likelihoods for all genotypes
     * [#patterns] * [#tips]
     */
    protected List<double[]>[][] genotypeLogLikelihoodsAll;

    /**
     * for variant calling and all tips
     * stores the log transfer probability from the parent with the maximum likelihood genotype (can be from different
     * matrices / categories) to the tip itself with all possible genotypes
     * [#tips] * [#genotypes (parent)]
     */
    protected List<double[]>[][] genotypeLogTransferProbabilityAll;

    /**
     * for variant calling and all tips
     * stores the log likelihoods for ADO states
     * [#patterns] * [#tips]
     */
    protected List<double[]>[][] adoLogLikelihoods;

    /**
     * a list of variant sites
     */
    protected List<VCFEntry> vcfEntries;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public VariantsInfoVCF() {
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initialise(int numOfPatterns, ScsSubstitutionModelBase substModel) {
        this.substModel = substModel;

        this.MLCategory = new ArrayList[numOfPatterns];
        this.MLGenotypesTips = new ArrayList[numOfPatterns];
        this.MLGenotypesNodes = new ArrayList[numOfPatterns];
        this.MLAdo = new ArrayList[numOfPatterns];
        this.genotypeLogLikelihoodConstantPattern = new ArrayList[numOfPatterns];
        this.genotypeLogLikelihoodsTips = new ArrayList[numOfPatterns];
        this.genotypeLogLikelihoodsNodes = new ArrayList[numOfPatterns];

        this.vcfEntries = new ArrayList<>();

        this.tree = ((ScsTree) treeInput.get()).copy();
        ((ScsTree) this.tree).initAndValidate();
        this.numOfTips = this.tree.getLeafNodeCount();
        this.numOfNodes = this.tree.getNodeCount();

        this.genotypes = new ArrayList[this.numOfNodes];

        this.genotypeLogLikelihoodsAll = new ArrayList[numOfPatterns][this.numOfTips];
        this.genotypeLogTransferProbabilityAll = new ArrayList[this.numOfTips][];
        this.adoLogLikelihoods = new ArrayList[numOfPatterns][this.numOfTips];

        this.sortedTipNames = new String[this.numOfTips];
        String[] original = new String[this.numOfTips];
        this.originalTipNamesIndices = new int[this.numOfTips];
        Arrays.fill(this.originalTipNamesIndices, -1);

        Set<String> existedTipNames = new HashSet<>();
        for (Node i : this.tree.getExternalNodes()) {
            this.sortedTipNames[i.getNr()] = i.getID();
            existedTipNames.add(i.getID());
        }

        if (existedTipNames.size() != this.numOfTips)
            throw new RuntimeException("Duplicated cell names found.");

        System.arraycopy(this.sortedTipNames, 0, original, 0, this.numOfTips);
        Arrays.sort(this.sortedTipNames);
        for (int i = 0; i < this.numOfTips; i++) {
            for (int j = 0; j < this.numOfTips; j++) {
                if (this.sortedTipNames[i].equals(original[j])) {
                    this.originalTipNamesIndices[i] = j;
                    break;
                }
            }
        }

        VCFEntry.addCellNames(
                Arrays.asList(this.sortedTipNames),
                VCFEntry.getCellNames().size() == 0
        );
    } // initialise

    @Override
    public void deeplyInitialise() {
        throw new IllegalArgumentException("Unsupported function.");
    } // deeplyInitialise

    @Override
    public void addMLCategory(int patternIndex, List<Integer> values, boolean overwrite) {
        if (this.MLCategory[patternIndex] == null)
            this.MLCategory[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLCategory[patternIndex].clear();

        this.MLCategory[patternIndex].addAll(values);
    } // addMLCategory

    @Override
    public void addMLCategory(int patternIndex, int value, boolean overwrite) {
        if (this.MLCategory[patternIndex] == null)
            this.MLCategory[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLCategory[patternIndex].clear();

        this.MLCategory[patternIndex].add(value);
    } // addMLCategory

    @Override
    public void addVariantsInfo(int index, Base value) {
        throw new IllegalArgumentException("Unsupported function.");
    } // addVariantsInfo

    @Override
    public List<Integer> getMLCategory(int patternIndex) {
        return this.MLCategory[patternIndex];
    } // getMLCategory

    @Override
    public void initialiseAPattern(int patternIndex) {
        if (this.MLGenotypesTips[patternIndex] == null)
            this.MLGenotypesTips[patternIndex] = new ArrayList<>();

        if (this.MLGenotypesNodes[patternIndex] == null)
            this.MLGenotypesNodes[patternIndex] = new ArrayList<>();

        if (this.MLAdo[patternIndex] == null)
            this.MLAdo[patternIndex] = new ArrayList<>();

        if (this.genotypeLogLikelihoodConstantPattern[patternIndex] == null)
            this.genotypeLogLikelihoodConstantPattern[patternIndex] = new ArrayList<>();

        if (this.genotypeLogLikelihoodsTips[patternIndex] == null)
            this.genotypeLogLikelihoodsTips[patternIndex] = new ArrayList<>();

        if (this.genotypeLogLikelihoodsNodes[patternIndex] == null)
            this.genotypeLogLikelihoodsNodes[patternIndex] = new ArrayList<>();
    } // initialiseAPattern

    /**
     * Add maximum likelihood genotypes of tips.
     *
     * @param patternIndex which pattern?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addMLGenotypesTips(int patternIndex, final int[] values, boolean overwrite) {
        if (this.MLGenotypesTips[patternIndex] == null)
            this.MLGenotypesTips[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLGenotypesTips[patternIndex].clear();

        if (values.length != this.numOfTips)
            throw new IllegalArgumentException("The expected number of tips is " + this.numOfTips + ", but " +
                    values.length + " tips are provided.");

        this.MLGenotypesTips[patternIndex].add(values);
    } // addMLGenotypesTips

    @Override
    public void addMLGenotypesTips(int patternIndex, final List<int[]> values, boolean overwrite) {
        if (this.MLGenotypesTips[patternIndex] == null)
            this.MLGenotypesTips[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLGenotypesTips[patternIndex].clear();

        if (values.get(0).length != this.numOfTips)
            throw new IllegalArgumentException("The expected number of tips is " + this.numOfTips + ", but " +
                    values.get(0).length + " tips are provided.");

        this.MLGenotypesTips[patternIndex].addAll(values);
    } // addMLGenotypesTips

    /**
     * Add maximum likelihood genotypes of all nodes.
     *
     * @param patternIndex which pattern?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addMLGenotypesNodes(int patternIndex, final int[] values, boolean overwrite) {
        if (this.MLGenotypesNodes[patternIndex] == null)
            this.MLGenotypesNodes[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLGenotypesNodes[patternIndex].clear();

        if (values.length != this.numOfNodes)
            throw new IllegalArgumentException("The expected number of nodes is " + this.numOfNodes + ", but " +
                    values.length + " nodes are provided.");

        this.MLGenotypesNodes[patternIndex].add(values);
    } // addMLGenotypesNodes

    @Override
    public void addMLGenotypesNodes(int patternIndex, final List<int[]> values, boolean overwrite) {
        if (this.MLGenotypesNodes[patternIndex] == null)
            this.MLGenotypesNodes[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLGenotypesNodes[patternIndex].clear();

        if (values.get(0).length != this.numOfNodes)
            throw new IllegalArgumentException("The expected number of nodes is " + this.numOfNodes + ", but " +
                    values.get(0).length + " nodes are provided.");

        this.MLGenotypesNodes[patternIndex].addAll(values);
    } // addMLGenotypesNodes

    /**
     * Add maximum likelihood ADO state of tips.
     *
     * @param patternIndex which pattern?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addMLAdo(int patternIndex, final int[] values, boolean overwrite) {
        if (this.MLAdo[patternIndex] == null)
            this.MLAdo[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLAdo[patternIndex].clear();

        if (values.length != this.numOfTips)
            throw new IllegalArgumentException("The expected number of tips is " + this.numOfTips + ", but " +
                    values.length + " tips are provided.");

        this.MLAdo[patternIndex].add(values);
    } // addMLAdo

    @Override
    public void addMLAdo(int patternIndex, final List<int[]> values, boolean overwrite) {
        if (this.MLAdo[patternIndex] == null)
            this.MLAdo[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.MLAdo[patternIndex].clear();

        if (values.get(0).length != this.numOfTips)
            throw new IllegalArgumentException("The expected number of tips is " + this.numOfTips + ", but " +
                    values.get(0).length + " tips are provided.");

        this.MLAdo[patternIndex].addAll(values);
    } // addMLAdo

    /**
     * Add log likelihood of a pattern being constant (matrix wise).
     *
     * @param patternIndex which pattern?
     * @param value        apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addGenotypeLogLikelihoodConstantPattern(int patternIndex, double value, boolean overwrite) {
        if (this.genotypeLogLikelihoodConstantPattern[patternIndex] == null)
            this.genotypeLogLikelihoodConstantPattern[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodConstantPattern[patternIndex].clear();

        this.genotypeLogLikelihoodConstantPattern[patternIndex].add(value);
    } // addGenotypeLogLikelihoodConstantPattern

    @Override
    public void addGenotypeLogLikelihoodConstantPattern(int patternIndex, List<Double> values, boolean overwrite) {
        if (this.genotypeLogLikelihoodConstantPattern[patternIndex] == null)
            this.genotypeLogLikelihoodConstantPattern[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodConstantPattern[patternIndex].clear();

        this.genotypeLogLikelihoodConstantPattern[patternIndex].addAll(values);
    } // addGenotypeLogLikelihoodConstantPattern

    /**
     * Add log likelihoods of tips.
     *
     * @param patternIndex which pattern?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addGenotypeLogLikelihoodsTips(int patternIndex, final double[] values, boolean overwrite) {
        if (this.genotypeLogLikelihoodsTips[patternIndex] == null)
            this.genotypeLogLikelihoodsTips[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodsTips[patternIndex].clear();

        if (values.length != this.numOfTips)
            throw new IllegalArgumentException("The expected number of tips is " + this.numOfTips + ", but " +
                    values.length + " tips are provided.");

        this.genotypeLogLikelihoodsTips[patternIndex].add(values);
    } // addGenotypeLogLikelihoodsTips

    @Override
    public void addGenotypeLogLikelihoodsTips(int patternIndex, final List<double[]> values, boolean overwrite) {
        if (this.genotypeLogLikelihoodsTips[patternIndex] == null)
            this.genotypeLogLikelihoodsTips[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodsTips[patternIndex].clear();

        if (values.get(0).length != this.numOfTips)
            throw new IllegalArgumentException("The expected number of tips is " + this.numOfTips + ", but " +
                    values.get(0).length + " tips are provided.");

        this.genotypeLogLikelihoodsTips[patternIndex].addAll(values);
    } // addGenotypeLogLikelihoodsTips

    /**
     * Add log likelihoods of all nodes.
     *
     * @param patternIndex which pattern?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addGenotypeLogLikelihoodsNodes(int patternIndex, final double[] values, boolean overwrite) {
        if (this.genotypeLogLikelihoodsNodes[patternIndex] == null)
            this.genotypeLogLikelihoodsNodes[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodsNodes[patternIndex].clear();

        if (values.length != this.numOfNodes)
            throw new IllegalArgumentException("The expected number of nodes is " + this.numOfNodes + ", but " +
                    values.length + " nodes are provided.");

        this.genotypeLogLikelihoodsNodes[patternIndex].add(values);
    } // addGenotypeLogLikelihoodsTips

    @Override
    public void addGenotypeLogLikelihoodsNodes(int patternIndex, final List<double[]> values, boolean overwrite) {
        if (this.genotypeLogLikelihoodsNodes[patternIndex] == null)
            this.genotypeLogLikelihoodsNodes[patternIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodsNodes[patternIndex].clear();

        if (values.get(0).length != this.numOfNodes)
            throw new IllegalArgumentException("The expected number of nodes is " + this.numOfNodes + ", but " +
                    values.get(0).length + " nodes are provided.");

        this.genotypeLogLikelihoodsNodes[patternIndex].addAll(values);
    } // addGenotypeLogLikelihoodsTips

    /**
     * Add log likelihoods of all genotypes for a specific tip at a specific pattern
     *
     * @param patternIndex which pattern?
     * @param tipIndex     which tip?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addGenotypeLogLikelihoodsAll(int patternIndex, int tipIndex, final double[] values, boolean overwrite) {
        if (this.genotypeLogLikelihoodsAll[patternIndex][tipIndex] == null)
            this.genotypeLogLikelihoodsAll[patternIndex][tipIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodsAll[patternIndex][tipIndex].clear();

        this.genotypeLogLikelihoodsAll[patternIndex][tipIndex].add(values);
    } // addGenotypeLogLikelihoodsAll

    @Override
    public void addGenotypeLogLikelihoodsAll(int patternIndex, int tipIndex, final List<double[]> values, boolean overwrite) {
        if (this.genotypeLogLikelihoodsAll[patternIndex][tipIndex] == null)
            this.genotypeLogLikelihoodsAll[patternIndex][tipIndex] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogLikelihoodsAll[patternIndex][tipIndex].clear();

        this.genotypeLogLikelihoodsAll[patternIndex][tipIndex].addAll(values);
    } // addGenotypeLogLikelihoodsAll

    /**
     * Add log transfer probabilities of a parent node from its maximum likelihood genotype to all genotypes of a child tip.
     *
     * @param tipIndex       which tip?
     * @param parentGenotype maximum likelihood genotype of parent node
     * @param numOfGenotypes the number of genotypes
     * @param values         apparently
     * @param overwrite      overwrite current values or not
     */
    @Override
    public void addGenotypeLogTransferProbabilitiesAll(int tipIndex, int parentGenotype, int numOfGenotypes, double[] values, boolean overwrite) {
        if (this.genotypeLogTransferProbabilityAll[tipIndex] == null)
            this.genotypeLogTransferProbabilityAll[tipIndex] = new ArrayList[numOfGenotypes];

        if (this.genotypeLogTransferProbabilityAll[tipIndex][parentGenotype] == null)
            this.genotypeLogTransferProbabilityAll[tipIndex][parentGenotype] = new ArrayList<>();
        else if (overwrite)
            this.genotypeLogTransferProbabilityAll[tipIndex][parentGenotype].clear();

        if (values.length != numOfGenotypes)
            throw new IllegalArgumentException("Error: the length of transfer probability vector should be " +
                    numOfGenotypes + ", not " + values.length);

        this.genotypeLogTransferProbabilityAll[tipIndex][parentGenotype].add(values);
    } // addGenotypeLogTransferProbabilitiesAll

    /**
     * Add log likelihoods of ADOs.
     *
     * @param patternIndex which pattern?
     * @param tipIndex     which tip?
     * @param numOfAdos    the number of allowed ADOs
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    @Override
    public void addAdoLogLikelihoods(int patternIndex, int tipIndex, int numOfAdos, double[] values, boolean overwrite) {
        if (this.adoLogLikelihoods[patternIndex][tipIndex] == null)
            this.adoLogLikelihoods[patternIndex][tipIndex] = new ArrayList<>();
        else if (overwrite)
            this.adoLogLikelihoods[patternIndex][tipIndex].clear();

        if (values.length != numOfAdos)
            throw new IllegalArgumentException("Error: the length of ADO likelihoods vector should be " +
                    numOfAdos + ", not " + values.length);

        this.adoLogLikelihoods[patternIndex][tipIndex].add(values);
    } // addAdoLogLikelihoods

    /**
     * write header information
     *
     * @param vcfOut         vcf log stream
     * @param cellNamesOut   cell names log stream
     * @param probsOut       posterior probability log stream
     * @param treeOut        annotated tree log stream
     * @param allelicInfoOut allelic information log stream
     * @param sizeFactorOut  size factor log stream
     */
    @Override
    public void init(
            @NotNull PrintStream vcfOut,
            PrintStream cellNamesOut,
            PrintStream probsOut,
            PrintStream treeOut,
            PrintStream allelicInfoOut,
            PrintStream sizeFactorOut
    ) {
        initHeader(vcfOut, cellNamesOut, probsOut, treeOut, allelicInfoOut);

        if (allelicInfoOut != null)
            logAllelicInfoHelper(allelicInfoOut, 0);

        if (sizeFactorOut != null)
            logSizeFactors(sizeFactorOut);
    } // init

    /**
     * Helper function to log allelic information.
     *
     * @param allelicInfoOut  allelic information log stream
     * @param startLocusIndex starting locus index
     */
    public void logAllelicInfoHelper(
            @NotNull PrintStream allelicInfoOut,
            final int startLocusIndex
    ) {
        if (startLocusIndex > 0 && !logSiteWiseAllelicInfo)
            return;

        int[] cellNameIndicesInData = new int[numOfTips];
        for (int i = 0; i < numOfTips; i++)
            cellNameIndicesInData[i] = scsData.getTaxonIndex(sortedTipNames[i]);

        if (logSiteWiseAllelicInfo) {
            for (int locusIndex = 0; locusIndex < scsData.getLociNr(); locusIndex++) {
                final int patternIndex = scsData.getPatternIndex(locusIndex);
                final int mlCategory = MLCategory[patternIndex].get(0);
                final int matrixIndex = mlCategory >= allelicSeqCovArray.length ? 0 : mlCategory;

                for (int cellIndex = 0; cellIndex < numOfTips; cellIndex++) {
                    allelicInfoOut.print(allelicSeqCovArray[matrixIndex][startLocusIndex + locusIndex] * sizeFactors[cellNameIndicesInData[cellIndex]]);
                    allelicInfoOut.print(",");
                    allelicInfoOut.print(allelicSeqCovRawVarArray[matrixIndex][startLocusIndex + locusIndex] * Math.pow(sizeFactors[cellNameIndicesInData[cellIndex]], 2));

                    if (cellIndex < numOfTips - 1)
                        allelicInfoOut.print("\t");
                }

                if (locusIndex < scsData.getLociNr() - 1)
                    allelicInfoOut.println();
            }
        } else {
            for (int i = 0; i < numOfTips; i++) {
                allelicInfoOut.print(allelicSeqCovArray[0][0] * sizeFactors[cellNameIndicesInData[i]]);
                allelicInfoOut.print("\t");
                allelicInfoOut.print(allelicSeqCovRawVarArray[0][0] * Math.pow(sizeFactors[cellNameIndicesInData[i]], 2));

                if (i < numOfTips - 1)
                    allelicInfoOut.println();
            }
        }
    } // logAllelicInfoHelper

    /**
     * Helper function to log size factors.
     *
     * @param sizeFactorOut size factor log stream
     */
    public void logSizeFactors(@NotNull PrintStream sizeFactorOut) {
        int[] cellNameIndicesInData = new int[numOfTips];
        for (int i = 0; i < numOfTips; i++)
            cellNameIndicesInData[i] = scsData.getTaxonIndex(sortedTipNames[i]);

        for (int i = 0; i < numOfTips; i++) {
            sizeFactorOut.print(sizeFactors[cellNameIndicesInData[i]]);

            if (i < numOfTips - 1)
                sizeFactorOut.println();
        }
    } // sizeFactorOut

    /**
     * log each entry to PrintStream
     *
     * @param vcfOut       vcf log stream
     * @param adoOut       ado state log stream
     * @param genotypesOut inferred genotypes at tips log stream
     * @param ternaryOut   ternary log stream
     * @param lociInfoOut  loci information log stream
     * @param probsOut     posterior probability log stream
     */
    @Override
    public void log(
            @NotNull PrintStream vcfOut,
            PrintStream adoOut,
            PrintStream genotypesOut,
            PrintStream ternaryOut,
            PrintStream lociInfoOut,
            PrintStream probsOut
    ) {
        // gather VCF entries
        for (int locusIndex = 0; locusIndex < scsData.getLociNr(); locusIndex++) {
            final VariantSiteInfo locusInfo = scsData.getLociInfo(locusIndex);
            final int patternIndex = scsData.getPatternIndex(locusIndex);

            final int loopSize = adjustAmbiguousGenotypes(this.MLGenotypesTips[patternIndex].size());

            for (int itemIndex = 0; itemIndex < loopSize; itemIndex++) {

                int variantTips = 0;

                for (int i : this.MLGenotypesTips[patternIndex].get(itemIndex)) {
                    if (this.substModel.isVariant(i))
                        variantTips++;
                }

                // fail to pass the consensus filter
                if (variantTips < this.cellThreshold) continue;

                this.vcfEntries.add(
                        getVCFInfoAcrossCells(
                                locusInfo,
                                locusIndex,
                                patternIndex,
                                itemIndex
                        )
                );

                // Save genotypes of all nodes as meta data
                for (Node node : this.tree.getNodesAsArray()) {
                    final int nodeNr = node.getNr();
                    if (this.genotypes[nodeNr] == null)
                        this.genotypes[nodeNr] = new ArrayList<>();

                    this.genotypes[nodeNr].add(this.MLGenotypesNodes[patternIndex].get(itemIndex)[nodeNr]);
                }
            }

        }

        // log VCF entries to VCF file and other specified files
        for (VCFEntry i : this.vcfEntries) {
            vcfOut.println(i.toString());

            if (adoOut != null)
                adoOut.println(i.getAdoStateAsString("\t"));

            if (genotypesOut != null)
                genotypesOut.println(i.getCellLocusOriGenotypeAsString("\t"));

            if (ternaryOut != null)
                ternaryOut.println(i.getCellLocusTernaryGenotypeAsString("\t"));

            if (lociInfoOut != null)
                lociInfoOut.println(i.getLociInfoAsString("\t", ","));

            if (probsOut != null)
                probsOut.println(i.getCellLocusGenotypeLikelihoodAsString("\t", ","));
        }

        // add meta data to the tree
        for (Node node : this.tree.getNodesAsArray())
            node.setMetaData(META_DATA_GENOTYPES, this.genotypes[node.getNr()].toArray(new Integer[0]));

    } // log

    /**
     * Gather VCF information for a locus across all cells.
     *
     * @param locusInfo        locus information
     * @param locusIndex       which locus?
     * @param patternIndex     which pattern?
     * @param combIndex        which combination of variant calling results?
     * @return an VCFEntry instance
     */
    private VCFEntry getVCFInfoAcrossCells(
            VariantSiteInfo locusInfo,
            int locusIndex,
            int patternIndex,
            int combIndex
    ) {
        final int[] MLGenotypesTips = this.MLGenotypesTips[patternIndex].get(combIndex);
        final int[] MLGenotypesNodes = this.MLGenotypesNodes[patternIndex].get(combIndex);

        List<CandidateAltNuc> candidateAltNucs = new ArrayList<>();
        int totalAlleleNums = 0;
        int totalReadDepth = 0;

        boolean isLocusInfoRenewed = false;

        CellLocusGenotype[] cellLocusGenotypes = new CellLocusGenotype[this.numOfTips];

        // candidateAltNucs initially contains alt nucs from locusInfo
        for (char altNuc : locusInfo.getAltNucs()) {
            candidateAltNucs.add(
                    new CandidateAltNuc(altNuc)
            );
        }

        // loop through all cells
        for (int i = 0; i < this.numOfTips; i++) {
            // label of the cell on beast.tree
            final int tipIndex = this.originalTipNamesIndices[i];

            // taxon index of this cell in scsData
            final int taxonIndex = getTaxonIndex(this.sortedTipNames[i], this.scsData);

            // get genotype of the parent node
            final int parentGenotype = MLGenotypesNodes[this.tree.getNode(tipIndex).getParent().getNr()];

            // get number of genotypes
            final int numOfStates = this.substModel.getStateCount();

            // get original genotype
            final String genotypeOri = this.substModel.getAlphabetGenotype(MLGenotypesTips[tipIndex]);

            // get adapted genotype
            final String genotypeAdapted = this.substModel.getGenotypeForVCF(
                    MLGenotypesTips[tipIndex],
                    candidateAltNucs,
                    this.scsData.getAltNucs(taxonIndex, locusIndex),
                    '.'
            );

            // get locus-wise significant alt nucs
            char[] altNucs = locusInfo.getAltNucs();
            if (candidateAltNucs.size() != altNucs.length) {
                // locusInfo should be updated
                isLocusInfoRenewed = true;

                altNucs = new char[candidateAltNucs.size()];

                for (int j = 0; j < altNucs.length; j++) {
                    altNucs[j] = candidateAltNucs.get(j).getNuc();
                }

                locusInfo = new VariantSiteInfo(
                        locusInfo.getChromosome(),
                        locusInfo.getPosition(),
                        locusInfo.getRefNuc(),
                        altNucs
                );
            }

            // get cell read depth
            final int readDepth = this.scsData.getSequencingCoverageOfLocus(taxonIndex, locusIndex);

            // get cell allelic depth (ref nuc and all alt nucs in candidateAltNucs)
            final int[] alleleDepth = this.scsData.getAlleleDepth(taxonIndex, locusIndex, altNucs);

            // get genotype log likelihoods
            double[] logLikelihood = new double[numOfStates];
            for (int j = 0; j < numOfStates; j++) {
                logLikelihood[j] = this.genotypeLogLikelihoodsAll[patternIndex][tipIndex].get(combIndex)[j] +
                        this.genotypeLogTransferProbabilityAll[tipIndex][parentGenotype].get(combIndex)[j];
            }

            // update totalAlleleNums
            totalAlleleNums += this.substModel.getNrOfAlleles(MLGenotypesTips[tipIndex]);

            // update totalReadDepth
            totalReadDepth += readDepth;

            // save cellLocusGenotype information
            cellLocusGenotypes[i] = new CellLocusGenotype(
                    genotypeOri,
                    genotypeAdapted,
                    this.substModel.getTernaryCode(MLGenotypesTips[tipIndex]),
                    alleleDepth,
                    readDepth,
                    logLikelihood,
                    this.MLAdo[patternIndex].get(combIndex)[tipIndex],
                    this.adoLogLikelihoods[patternIndex][tipIndex].get(combIndex)
            );
        }

        // updating alleleDepth for each member in cellLocusGenotypes is necessary
        if (isLocusInfoRenewed) {
            for (int i = 0; i < cellLocusGenotypes.length; i++) {
                if (cellLocusGenotypes[i].getAlleleDepthLength() - 1 < locusInfo.getAltNucsLength())
                    cellLocusGenotypes[i].setAlleleDepth(
                            this.scsData.getAlleleDepth(
                                    getTaxonIndex(this.sortedTipNames[i], this.scsData),
                                    locusIndex,
                                    locusInfo.getAltNucs()
                            )
                    );
            }
        }

        // gather information for VCFInfo
        List<Integer> alleleCount = new ArrayList<>();
        for (CandidateAltNuc i : candidateAltNucs)
            alleleCount.add(i.getCount());

        // new VCFInfo instance
        VCFInfo vcfInfo = new VCFInfo(
                Ints.toArray(alleleCount),
                totalAlleleNums,
                totalReadDepth
        );

        // new VCFEntry instance
        return new VCFEntry(
                locusInfo,
                MISSING_VALUE_IN_VCF,
                convertLogE2RoundedPhredScaled(this.genotypeLogLikelihoodConstantPattern[patternIndex].get(combIndex)),
                "PASS",
                vcfInfo,
                cellLocusGenotypes
        );
    } // getVCFInfoAcrossCells

    /**
     * log each entry to PrintStream
     *
     * @param vcfOut vcf log stream
     */
    @Override
    public void close(@NotNull PrintStream vcfOut) {
        // nothing needs to be closed
    } // close

}
