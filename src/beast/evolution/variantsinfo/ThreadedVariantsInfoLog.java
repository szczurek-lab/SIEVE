package beast.evolution.variantsinfo;

import org.jetbrains.annotations.NotNull;

import java.io.PrintStream;
import java.util.List;

public class ThreadedVariantsInfoLog extends GenericVariantsInfoLog {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    VariantsInfoLog[] variantsInfoLogs;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public ThreadedVariantsInfoLog(int numOfThreads) {
        this.ambiguityGenotypesStrategy = AmbiguityGenotypesStrategy.Default;

        assert numOfThreads > 0;
        this.variantsInfoLogs = new VariantsInfoLog[numOfThreads];
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void deeplyInitialise() {
        this.substModel = this.variantsInfoLogs[0].substModel;

        this.numOfTips = this.variantsInfoLogs[0].numOfTips;
        this.numOfNodes = this.variantsInfoLogs[0].numOfNodes;

        this.sortedTipNames = new String[this.numOfTips];
        this.originalTipNamesIndices = new int[this.numOfTips];

        System.arraycopy(this.variantsInfoLogs[0].sortedTipNames, 0, this.sortedTipNames, 0, this.numOfTips);
        System.arraycopy(this.variantsInfoLogs[0].originalTipNamesIndices, 0, this.originalTipNamesIndices, 0, this.numOfTips);
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
    public void addVariantsInfo(int index, GenericVariantsInfo.Base value) {
        assert value instanceof VariantsInfoLog;
        this.variantsInfoLogs[index] = (VariantsInfoLog) value;
    } // addVariantsInfo

    @Override
    public List<Integer> getMLCategory(int patternIndex) {
        throw new IllegalArgumentException("Unsupported function.");
    } // getMLCategory

    @Override
    public void initialiseAPattern(int patternIndex) {
        for (VariantsInfoLog i : variantsInfoLogs)
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
     * Not `Loggable`.
     *
     * @param out apparently
     */
    @Override
    public void log(@NotNull PrintStream out, String siteSeparator, String cellSeparator) {
        for (int i = 0; i < this.variantsInfoLogs.length; i++) {
            this.variantsInfoLogs[i].log(out, siteSeparator, cellSeparator);

            if (i < this.variantsInfoLogs.length - 1)
                out.print("\t");
        }
    } // log

    /**
     * Not `Loggable`.
     *
     * @param out apparently
     */
    @Override
    public void cloze(@NotNull PrintStream out) {
        for (VariantsInfoLog i : this.variantsInfoLogs) i.cloze(out);
    } // close

}
