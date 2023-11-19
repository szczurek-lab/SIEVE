package beast.evolution.variantsinfo;

import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import beast.evolution.tree.Node;
import beast.evolution.tree.ScsTree;
import org.jetbrains.annotations.NotNull;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class ThreadedVariantsInfoVCF extends GenericVariantsInfoVCF {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    VariantsInfoVCF[] variantsInfoVCFs;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public ThreadedVariantsInfoVCF() {
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initialise(int numOfThreads, ScsSubstitutionModelBase substModel) {
        assert numOfThreads > 0;
        this.variantsInfoVCFs = new VariantsInfoVCF[numOfThreads];
        this.tree = ((ScsTree) treeInput.get()).copy();
        ((ScsTree) this.tree).initAndValidate();

        this.genotypes = new ArrayList[this.tree.getNodeCount()];
    } // initialise

    @Override
    public void deeplyInitialise() {
        this.substModel = this.variantsInfoVCFs[0].substModel;

        this.numOfTips = this.variantsInfoVCFs[0].numOfTips;
        this.numOfNodes = this.variantsInfoVCFs[0].numOfNodes;

        this.sortedTipNames = new String[this.numOfTips];
        this.originalTipNamesIndices = new int[this.numOfTips];

        System.arraycopy(this.variantsInfoVCFs[0].sortedTipNames, 0, this.sortedTipNames, 0, this.numOfTips);
        System.arraycopy(this.variantsInfoVCFs[0].originalTipNamesIndices, 0, this.originalTipNamesIndices, 0, this.numOfTips);
    } // deeplyInitialise

    @Override
    public void addMLCategory(int patternIndex, List<Integer> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
    } // addMLCategory

    @Override
    public void addMLCategory(int patternIndex, int value, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
    } // addMLCategory

    @Override
    public void addVariantsInfo(int index, Base value) {
        assert value instanceof VariantsInfoVCF;
        this.variantsInfoVCFs[index] = (VariantsInfoVCF) value;
    } // addVariantsInfo

    @Override
    public List<Integer> getMLCategory(int patternIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    } // getMLCategory

    @Override
    public void initialiseAPattern(int patternIndex) {
        for (VariantsInfoVCF i : variantsInfoVCFs)
            i.initialiseAPattern(patternIndex);
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
        throw new IllegalArgumentException("Unsupported function.");
    } // addMLGenotypesTips

    @Override
    public void addMLGenotypesTips(int patternIndex, final List<int[]> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
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
        throw new IllegalArgumentException("Unsupported function.");
    } // addMLGenotypesNodes

    @Override
    public void addMLGenotypesNodes(int patternIndex, final List<int[]> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
    } // addMLGenotypesNodes

    /**
     * Add maximum likelihood ADO state of tips.
     *
     * @param patternIndex which pattern?
     * @param values       apparently
     * @param overwrite    overwrite current values or not
     */
    public void addMLAdo(int patternIndex, final int[] values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
    } // addMLAdo

    public void addMLAdo(int patternIndex, final List<int[]> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
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
        throw new IllegalArgumentException("Unsupported function.");
    } // addGenotypeLogLikelihoodConstantPattern

    @Override
    public void addGenotypeLogLikelihoodConstantPattern(int patternIndex, List<Double> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
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
        throw new IllegalArgumentException("Unsupported function.");
    } // addGenotypeLogLikelihoodsTips

    @Override
    public void addGenotypeLogLikelihoodsTips(int patternIndex, final List<double[]> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
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
        throw new IllegalArgumentException("Unsupported function.");
    } // addGenotypeLogLikelihoodsNodes

    @Override
    public void addGenotypeLogLikelihoodsNodes(int patternIndex, final List<double[]> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
    } // addGenotypeLogLikelihoodsNodes

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
        throw new IllegalArgumentException("Unsupported function.");
    } // addGenotypeLogLikelihoodsAll

    @Override
    public void addGenotypeLogLikelihoodsAll(int patternIndex, int tipIndex, final List<double[]> values, boolean overwrite) {
        throw new IllegalArgumentException("Unsupported function.");
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
        throw new IllegalArgumentException("Unsupported function.");
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
        throw new IllegalArgumentException("Unsupported function.");
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

        if (allelicInfoOut != null) {
            int start = 0;

            for (VariantsInfoVCF i : this.variantsInfoVCFs) {
                i.logAllelicInfoHelper(allelicInfoOut, start);
                start += scsDataInput.get().getLociNr();
            }
        }

        if (sizeFactorOut != null)
            this.variantsInfoVCFs[0].logSizeFactors(sizeFactorOut);
    } // init

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
        for (int i = 0; i < this.variantsInfoVCFs.length; i++) {
            this.variantsInfoVCFs[i].log(vcfOut, adoOut, genotypesOut, ternaryOut, lociInfoOut, probsOut);

            // copy tree meta data
            for (Node node : this.tree.getNodesAsArray()) {
                final int nodeNr = node.getNr();
                if (this.genotypes[nodeNr] == null)
                    this.genotypes[nodeNr] = new ArrayList<>();

                this.genotypes[nodeNr].addAll(this.variantsInfoVCFs[i].genotypes[nodeNr]);
            }

            if (i < this.variantsInfoVCFs.length - 1) {
                vcfOut.println();

                if (adoOut != null)
                    adoOut.println();

                if (genotypesOut != null)
                    genotypesOut.println();

                if (ternaryOut != null)
                    ternaryOut.println();

                if (lociInfoOut != null)
                    lociInfoOut.println();
            }
        }

        // add meta data to the tree
        for (Node node : this.tree.getNodesAsArray())
            node.setMetaData(META_DATA_GENOTYPES, this.genotypes[node.getNr()].toArray(new Integer[0]));
    } // log

    /**
     * log each entry to PrintStream
     *
     * @param vcfOut vcf log stream
     */
    @Override
    public void close(@NotNull PrintStream vcfOut) {
        for (VariantsInfoVCF i : this.variantsInfoVCFs) i.close(vcfOut);
    } // close

}
