package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.rawreadcountsmodel.RawReadCountsModelInterface;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.TreeInterface;
import beast.evolution.variantsinfo.GenericVariantsInfo;

import java.io.PrintStream;
import java.util.List;
import java.util.Random;

@Description("Abstract class of tree likelihood for SCS.")
public abstract class ScsGenericTreeLikelihood extends Distribution {


    //**********************************************
    //*                   Inputs                   *
    //**********************************************

    final public Input<ScsAlignment> scsDataInput = new Input<>("scsData", "single-cell sequence data " +
            "for the beast.tree", Input.Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Input.Validate.REQUIRED);

    final public Input<SiteModelInterface.Base> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Input.Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    public final Input<RawReadCountsModelInterface.Base> rawReadCountsModelInput = new Input<>("rawReadCountsModel", "a raw read counts model applied to leaves of a tree, accounting for technical errors", Input.Validate.REQUIRED);

    public final Input<Boolean> runTimeAnalysisInput = new Input<>("runTimeAnalysis", "analyze the running time in likelihood computation", Input.Validate.OPTIONAL);

    final public Input<Boolean> useLogPartialsInput = new Input<>("useLogPartials", "whether to use log-partials " +
            "when computing likelihood (default true)");

    final public Input<GenericVariantsInfo.Base> variantsInfoInput = new Input<>("variantsInfo",
            "collection of variants information (only in variant calling mode)", Input.Validate.OPTIONAL);

    final public Input<Boolean> useOnlyBranchLengthInput = new Input<>("useOnlyBranchLength",
            "in variant calling whether only using tree branch length or also using inferred rate " +
                    "(together with non-strict molecular clock model)",
            true, Input.Validate.OPTIONAL);

    final public Input<RealParameter> meanRateInput = new Input<>("meanRate",
            "mean clock rate (defaults to 1.0); should only be used in variant calling mode and be the " +
                    "same variable as that in the corresponding branch rate model.");

    final public Input<Boolean> traceMLGenotypesInput = new Input<>("traceMLGenotypes",
            "tracing maximum likelihood genotypes and ADO states during MCMC. A logger " +
                    "(`GenotypeAdoStateLogger`) to log genotypes and ADO states should also be defined.",
            false, Input.Validate.OPTIONAL);


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    protected boolean useLogPartials;

    /**
     * whether to correct ascertainment bias or not
     */
    protected boolean useAscBiasCorrection;

    protected double biasCorr;

    /**
     * whether is in variant calling mode or not
     */
    protected boolean inVariantCallingMode = false;

    protected boolean useOnlyBranchLength = true;

    /**
     * should only be used in variant calling mode and be the same variable as that in the corresponding branch rate model
     */
    protected double meanRate = 1.0;

    /**
     * for variant calling
     */
    protected GenericVariantsInfo.Base variantsInfo;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return null;
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return null;
    }

    /**
     * This method draws new values for the arguments conditional on the current value(s) of the conditionals.
     * <p/>
     * The new values are overwrite the argument values in the provided state.
     *
     * @param state  the state
     * @param random random number generator
     */
    @Override
    public void sample(State state, Random random) {
    }

    public boolean updateSeqCovModel() {
        return rawReadCountsModelInput.get().updateSeqCovModel();
    }

    /**
     * for the purpose of debug
     * must be overridden
     *
     * @param inDebugMode in debug mode or not
     */
    abstract public void updateInDebugMode(final boolean inDebugMode);

    /**
     * Calling variants.
     */
    abstract public void callVariants();

    /**
     * Is tracing maximum likelihood genotypes or not.
     *
     * @return yes or no
     */
    abstract public boolean isTraceMLGenotypes();

    /**
     * initialize header of logger
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    abstract public void initSiteHeader(int start, PrintStream out);

    abstract public String getSortedCellNamesFromVarInfo(String separator);


    //******************************************
    //*           Logger methods for           *
    //*        genotypes and ADO states        *
    //******************************************

    /**
     * log sampled allelic sequencing coverage and raw variance for each site
     *
     * @param out           apparently
     * @param siteSeparator apparently
     * @param cellSeparator apparently
     */
    public void logGenotypeAndAdo(PrintStream out, String siteSeparator, String cellSeparator) {
        callVariants();
        variantsInfo.log(out, siteSeparator, cellSeparator);
    } // logGenotypeAndAdo

    /**
     * close
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    abstract public void closeGenotypeAndAdo(int start, PrintStream out);


    //*******************************************
    //*           Logger methods for            *
    //*    allelic coverage and raw variance    *
    //*******************************************

    /**
     * log sampled allelic sequencing coverage and raw variance for each site
     *
     * @param out apparently
     */
    abstract public void logCovar(PrintStream out);

    /**
     * close
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    abstract public void closeCovar(int start, PrintStream out);

}
