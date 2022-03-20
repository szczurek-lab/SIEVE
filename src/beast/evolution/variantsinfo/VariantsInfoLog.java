package beast.evolution.variantsinfo;

import beast.evolution.alignment.ScsAlignment;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import org.jetbrains.annotations.NotNull;

import java.io.PrintStream;
import java.util.*;

public class VariantsInfoLog extends GenericVariantsInfoLog {

    private final static String ITEM_SEPARATOR = ",";


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
     * for variant calling and all tips
     * The maximum likelihood ADO state for each tip
     * [#patterns]
     */
    protected List<int[]>[] MLAdo;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public VariantsInfoLog(
            ScsAlignment scsData,
            TreeInterface tree,
            ScsSubstitutionModelBase substModel
    ) {
        this.ambiguityGenotypesStrategy = AmbiguityGenotypesStrategy.Default;
        this.scsData = scsData;
        this.tree = tree;
        this.substModel = substModel;

        final int numOfPatterns = scsData.getPatternCount();

        this.MLCategory = new ArrayList[numOfPatterns];
        this.MLGenotypesTips = new ArrayList[numOfPatterns];
        this.MLAdo = new ArrayList[numOfPatterns];

        this.numOfTips = this.tree.getLeafNodeCount();
        this.numOfNodes = this.tree.getNodeCount();

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
            throw new RuntimeException("ERROR! Duplicated cell names found.");

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
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

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
    public void addVariantsInfo(int index, GenericVariantsInfo.Base value) {
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

        if (this.MLAdo[patternIndex] == null)
            this.MLAdo[patternIndex] = new ArrayList<>();
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
    } // addMLGenotypesNodes

    @Override
    public void addMLGenotypesNodes(int patternIndex, final List<int[]> values, boolean overwrite) {
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
    } // addGenotypeLogLikelihoodConstantPattern

    @Override
    public void addGenotypeLogLikelihoodConstantPattern(int patternIndex, List<Double> values, boolean overwrite) {
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
    } // addGenotypeLogLikelihoodsTips

    @Override
    public void addGenotypeLogLikelihoodsTips(int patternIndex, final List<double[]> values, boolean overwrite) {
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
    } // addGenotypeLogLikelihoodsTips

    @Override
    public void addGenotypeLogLikelihoodsNodes(int patternIndex, final List<double[]> values, boolean overwrite) {
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
    } // addGenotypeLogLikelihoodsAll

    @Override
    public void addGenotypeLogLikelihoodsAll(int patternIndex, int tipIndex, final List<double[]> values, boolean overwrite) {
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
    } // addAdoLogLikelihoods

    /**
     * Not `Loggable`.
     *
     * @param out apparently
     */
    @Override
    public void log(@NotNull PrintStream out, String siteSeparator, String cellSeparator) {
        StringBuilder sb = new StringBuilder();

        for (int locusIndex = 0; locusIndex < scsData.getLociNr(); locusIndex++) {
            final int patternIndex = scsData.getPatternIndex(locusIndex);

            sb.append(getInfoAcrossCells(patternIndex, cellSeparator));

            if (locusIndex < scsData.getLociNr() - 1)
                sb.append(siteSeparator);
        }

        out.print(sb);
    } // log

    private String getInfoAcrossCells(int patternIndex, String separator) {
        StringBuilder sb = new StringBuilder();

        final int[] MLGenotypesTips = this.MLGenotypesTips[patternIndex].get(0);
        final int[] MLAdo = this.MLAdo[patternIndex].get(0);

        // loop through all cells
        for (int i = 0; i < this.numOfTips; i++) {
            // label of the cell on beast.tree
            final int tipIndex = this.originalTipNamesIndices[i];

            // get original genotype
            final String genotypeOri = this.substModel.getAlphabetGenotype(MLGenotypesTips[tipIndex]);

            sb.append(genotypeOri).append(ITEM_SEPARATOR).append(MLAdo[tipIndex]);

            if (i < this.numOfTips - 1)
                sb.append(separator);
        }

        return sb.toString();
    } // getInfoAcrossCells

}
