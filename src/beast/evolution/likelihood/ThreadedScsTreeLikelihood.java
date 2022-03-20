package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.app.beauti.Beauti;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.FilteredScsAlignment;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.rawreadcountsmodel.RawReadCountsModelInterface;
import beast.evolution.rawreadcountsmodel.nucreadcountsmodel.NucReadCountsModelInterface;
import beast.evolution.rawreadcountsmodel.seqcovmodel.SeqCovModelInterface;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.TreeInterface;
import beast.evolution.variantsinfo.GenericVariantsInfoVCF;
import beast.evolution.variantsinfo.ThreadedVariantsInfoLog;
import beast.evolution.variantsinfo.ThreadedVariantsInfoVCF;
import beast.evolution.variantsinfo.VariantsInfoVCF;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import static beast.math.util.MathFunctions.logSumExp;

@Description("multithreading tree likelihood")
public class ThreadedScsTreeLikelihood extends ScsGenericTreeLikelihood {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads", "maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

    final public Input<String> proportionsInput = new Input<>("proportions", "specifies proportions of patterns used per thread as space "
            + "delimited string. This is useful when using a mixture of BEAGLE devices that run at different speeds, e.g GPU and CPU. "
            + "The string is duplicated if there are more threads than proportions specified. For example, "
            + "'1 2' as well as '33 66' with 2 threads specifies that the first thread gets a third of the patterns and the second "
            + "two thirds. With 3 threads, it is interpreted as '1 2 1' = 25%, 50%, 25% and with 7 threads it is "
            + "'1 2 1 2 1 2 1' = 10% 20% 10% 20% 10% 20% 10%. If not specified, all threads get the same proportion of patterns.");

    enum Scaling {none, always, _default}

    final public Input<Scaling> scalingInput = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(ThreadedScsTreeLikelihood.Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", ThreadedScsTreeLikelihood.Scaling._default, ThreadedScsTreeLikelihood.Scaling.values());

    /**
     * private list of likelihoods, to notify framework of ScsTreeLikelihoods being created in initAndValidate()
     **/
    final private Input<List<ScsTreeLikelihood>> likelihoodsInput = new Input<>("*", "", new ArrayList<>());

    @Override
    public List<Input<?>> listInputs() {
        List<Input<?>> list = super.listInputs();
        if (!Beauti.isInBeauti() && System.getProperty("beast.is.junit.testing") == null) {
            // do not expose internal likelihoods to BEAUti or junit tests
            list.add(likelihoodsInput);
        }
        return list;
    }

    /**
     * calculation engine
     **/
    private ScsTreeLikelihood[] treeLikelihood;

    private ExecutorService pool = null;
    private final List<Callable<Double>> likelihoodCallers = new ArrayList<>();
    private List<Callable<Double>> MLGenotypesCallers = null;

    /**
     * number of threads to use, changes when threading causes problems
     **/
    private int threadCount;
    private double[] logPByThread;
    private double[] constSumByThread;

    // specified a set ranges of patterns assigned to each thread
    // first patternPoints contains 0, then one point for each thread
    private int[] patternPoints;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        threadCount = BeastMCMC.m_nThreads;

        if (maxNrOfThreadsInput.get() > 0)
            threadCount = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);

        String instanceCount = System.getProperty("beast.instance.count");
        if (instanceCount != null && instanceCount.length() > 0)
            threadCount = Integer.parseInt(instanceCount);

        useAscBiasCorrection = scsDataInput.get().isAscBiasCorrection();

        logPByThread = new double[threadCount];
        constSumByThread = new double[threadCount];

        if (useLogPartialsInput.get() != null)
            useLogPartials = useLogPartialsInput.get();
        else
            useLogPartials = true;

        // sanity check: alignment should have same #taxa as tree
        if (scsDataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount())
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");

        // deeply initialize raw read counts model
        try {
            rawReadCountsModelInput.get().deeplyInitialize(
                    scsDataInput.get(),
                    siteModelInput.get().getCategoryCount(),
                    siteModelInput.get().substModelInput.get().getStateCount(),
                    ((ScsSubstitutionModelBase) (siteModelInput.get().substModelInput.get())).getModeledAlleles(),
                    useLogPartials
            );
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        if (this.rawReadCountsModelInput.get().isInVariantCallingMode()) {
            this.inVariantCallingMode = true;
            MLGenotypesCallers = new ArrayList<>();

            if (variantsInfoInput.get() == null)
                throw new IllegalArgumentException("'variantsInfo' is missing in variant calling mode.");

            variantsInfo = variantsInfoInput.get();

            if (!(variantsInfo instanceof ThreadedVariantsInfoVCF))
                throw new IllegalArgumentException("'variantsInfo' should be of type 'beast.evolution.variantsinfo.ThreadedVariantsInfoVCF'");

            variantsInfo.initialise(threadCount, null);
        } else if (this.traceMLGenotypesInput.get()) {
            MLGenotypesCallers = new ArrayList<>();
            variantsInfo = new ThreadedVariantsInfoLog(threadCount);
        }

        treeLikelihood = new ScsTreeLikelihood[threadCount];

        if (threadCount <= 1) {
            treeLikelihood[0] = new ScsTreeLikelihood();
            treeLikelihood[0].setID(getID() + "0");
            treeLikelihood[0].initByName(
                    "scsData", scsDataInput.get(),
                    "tree", treeInput.get(),
                    "siteModel", siteModelInput.get(),
                    "branchRateModel", branchRateModelInput.get(),
                    "rawReadCountsModel", rawReadCountsModelInput.get(),
                    "runTimeAnalysis", runTimeAnalysisInput.get(),
                    "useLogPartials", useLogPartials,
                    "variantsInfo", this.inVariantCallingMode ? createVariantsInfo(0, scsDataInput.get(), treeInput.get()) : null,
                    "useOnlyBranchLength", this.inVariantCallingMode ? useOnlyBranchLengthInput.get() : null,
                    "meanRate", this.inVariantCallingMode ? meanRateInput.get() : null,
                    "traceMLGenotypes", traceMLGenotypesInput.get(),
                    "scaling", scalingInput.get() + ""
            );
            treeLikelihood[0].getOutputs().add(this);
            likelihoodsInput.get().add(treeLikelihood[0]);

            if (traceMLGenotypesInput.get())
                variantsInfo.addVariantsInfo(0, treeLikelihood[0].getVariantsInfo());
        } else {
            pool = Executors.newFixedThreadPool(threadCount);

            calcPatternPoints(scsDataInput.get().getSiteCount());
            for (int i = 0; i < threadCount; i++) {
                String filterSpec = (patternPoints[i] + 1) + "-" + (patternPoints[i + 1]);

                treeLikelihood[i] = new ScsTreeLikelihood();
                treeLikelihood[i].setID(getID() + i);
                treeLikelihood[i].getOutputs().add(this);
                likelihoodsInput.get().add(treeLikelihood[i]);

                FilteredScsAlignment filter = new FilteredScsAlignment();
                filter.initByName(
                        "scsData", scsDataInput.get(),
                        "filter", filterSpec,
                        "dataType", scsDataInput.get().dataTypeInput.get(),
                        "ascertained", scsDataInput.get().ascBiasCorrectionInput.get(),
                        "meanAscBiasCorrection", scsDataInput.get().meanAscBiasCorrectionInput.get(),
                        "nrOfBackgroundSites", scsDataInput.get().getNrOfBackgroundSites()
                );

                treeLikelihood[i].initByName(
                        "scsData", filter,
                        "tree", treeInput.get(),
                        "siteModel", duplicate(siteModelInput.get(), filter, i),
                        "branchRateModel", duplicate(branchRateModelInput.get(), filter, i),
                        "rawReadCountsModel", duplicate(rawReadCountsModelInput.get(), filter, i),
                        "runTimeAnalysis", runTimeAnalysisInput.get(),
                        "useLogPartials", useLogPartials,
                        "variantsInfo", this.inVariantCallingMode ? createVariantsInfo(i, filter, treeInput.get()) : null,
                        "useOnlyBranchLength", this.inVariantCallingMode ? useOnlyBranchLengthInput.get() : null,
                        "meanRate", this.inVariantCallingMode ? meanRateInput.get() : null,
                        "traceMLGenotypes", traceMLGenotypesInput.get(),
                        "scaling", scalingInput.get() + ""
                );

                if (traceMLGenotypesInput.get())
                    variantsInfo.addVariantsInfo(i, treeLikelihood[i].getVariantsInfo());

                likelihoodCallers.add(new ThreadedScsTreeLikelihood.ScsTreeLikelihoodCaller(treeLikelihood[i], i));
                if (this.inVariantCallingMode || traceMLGenotypesInput.get())
                    MLGenotypesCallers.add(new ThreadedScsTreeLikelihood.MLGenotypeCaller(treeLikelihood[i], i));
            }
        }

        if (this.inVariantCallingMode || traceMLGenotypesInput.get())
            variantsInfo.deeplyInitialise();

    } // initAndValidate

    private void calcPatternPoints(int nSites) {
        patternPoints = new int[threadCount + 1];
        if (proportionsInput.get() == null) {
            int[] counts = new int[threadCount];
            Arrays.fill(counts, nSites / threadCount);

            // make sure the sites are splitting as evenly as possible
            for (int i = 0; i < nSites % threadCount; i++)
                counts[i]++;

            for (int i = 0; i < threadCount; i++)
                patternPoints[i + 1] = patternPoints[i] + counts[i];
        } else {
            String[] strs = proportionsInput.get().split("\\s+");
            double[] proportions = new double[threadCount];
            for (int i = 0; i < threadCount; i++)
                proportions[i] = Double.parseDouble(strs[i % strs.length]);

            // normalise
            double sum = 0;
            for (double d : proportions)
                sum += d;

            for (int i = 0; i < threadCount; i++)
                proportions[i] /= sum;

            // cumulative
            for (int i = 1; i < threadCount; i++)
                proportions[i] += proportions[i - 1];

            // calc ranges
            for (int i = 0; i < threadCount; i++)
                patternPoints[i + 1] = (int) (proportions[i] * nSites + 0.5);
        }

        // sanity check: there is no overlaps for any ranges
        for (int i = 0; i < patternPoints.length - 1; i++) {
            if (patternPoints[i] >= patternPoints[i + 1]) {
                throw new IllegalArgumentException("Error: overlaps found when partitioning sites. Redefine the " +
                        "'proportions' or just leave it blank.");
            }
        }
    } // calcPatternPoints

    private VariantsInfoVCF createVariantsInfo(
            int index,
            ScsAlignment data,
            TreeInterface tree
    ) {
        if (!(this.variantsInfo instanceof GenericVariantsInfoVCF)) return null;

        VariantsInfoVCF variantsInfo = new VariantsInfoVCF();
        variantsInfo.setID("variantsInfo" + index);
        variantsInfo.getOutputs().add(this.variantsInfo);

        variantsInfo.initByName(
                "scsData", data,
                "seqCovModel", ((GenericVariantsInfoVCF) this.variantsInfo).seqCovModelInput.get(),
                "tree", tree,
                "cellThreshold", ((GenericVariantsInfoVCF) this.variantsInfo).cellThresholdInput.get() == null ? null : ((GenericVariantsInfoVCF) this.variantsInfo).cellThresholdInput.get(),
                "ambiguity", ((GenericVariantsInfoVCF) this.variantsInfo).ambiguityGenotypesStrategyInput.get() == null ? null : ((GenericVariantsInfoVCF) this.variantsInfo).ambiguityGenotypesStrategyInput.get()
        );

        this.variantsInfo.addVariantsInfo(index, variantsInfo);

        return variantsInfo;
    } // createVariantsInfo

    /**
     * create new instance of src object, connecting all inputs from src object
     * Note if input is a SubstModel, it is duplicated as well.
     * Note if input is a raw read counts model, some critical variables are duplicated as well.
     *
     * @param src        object to be copied
     * @param targetData data corresponding to the target raw read counts model
     * @param i          index used to extend ID with.
     * @return copy of src object
     */
    private Object duplicate(BEASTInterface src, final ScsAlignment targetData, int i) {
        if (src == null) {
            return null;
        }

        BEASTInterface copy;
        try {
            copy = src.getClass().newInstance();
            copy.setID(src.getID() + "_" + i);
        } catch (InstantiationException | IllegalAccessException e) {
            e.printStackTrace();
            throw new RuntimeException("Programmer error: every object in the model should have a default constructor that is publicly accessible: " + src.getClass().getName());
        }

        for (Input<?> input : src.listInputs()) {
            if (input.get() != null) {
                if (input.get() instanceof List) {
                    // handle lists
                    for (Object o : (List<?>) input.get()) {
                        if (o instanceof BEASTInterface) {
                            // make sure it is not already in the list
                            copy.setInputValue(input.getName(), o);
                        }
                    }
                } else if (input.get() instanceof SubstitutionModel) {
                    // duplicate subst models
                    BEASTInterface substModel = (BEASTInterface) duplicate((BEASTInterface) input.get(), targetData, i);
                    copy.setInputValue(input.getName(), substModel);
                } else if (input.get() instanceof SeqCovModelInterface.Base) {
                    // duplicate sequencing coverage models
                    BEASTInterface seqCovModel = (BEASTInterface) duplicate((BEASTInterface) input.get(), targetData, i);
                    copy.setInputValue(input.getName(), seqCovModel);
                } else if (input.get() instanceof NucReadCountsModelInterface.Base) {
                    // duplicate nucleotide read counts models
                    BEASTInterface nucReadCountsModel = (BEASTInterface) duplicate((BEASTInterface) input.get(), targetData, i);
                    copy.setInputValue(input.getName(), nucReadCountsModel);
                } else {
                    // it is some other value
                    copy.setInputValue(input.getName(), input.get());
                }
            }
        }

        if (src instanceof RawReadCountsModelInterface.Base) {
            assert targetData instanceof FilteredScsAlignment;
            ((RawReadCountsModelInterface.Base) src).duplicate(
                    (RawReadCountsModelInterface.Base) copy,
                    targetData
            );
        }

        copy.initAndValidate();
        return copy;
    } // duplicate

    @Override
    public double calculateLogP() {
        logP = calculateLogPByBeagle();
        return logP;
    } // calculateLogP

    private double calculateLogPByBeagle() {
        try {
            if (threadCount > 1) {
                pool.invokeAll(likelihoodCallers);

                logP = 0;
                for (double f : logPByThread)
                    logP += f;

                if (useAscBiasCorrection) {
                    biasCorr = scsDataInput.get().getAscBiasCorrection(constSumByThread);

                    logP += biasCorr;
                }

            } else
                logP = treeLikelihood[0].calculateLogP();
        } catch (RejectedExecutionException | InterruptedException e) {
            e.printStackTrace();
            System.exit(0);
        }

        return logP;
    } // calculateLogPByBeagle

    /* return copy of pattern log likelihoods for each of the patterns in the alignment */
    public double[] getPatternLogLikelihoods() {
        double[] patternLogLikelihoods = new double[scsDataInput.get().getPatternCount()];
        int i = 0;
        for (ScsTreeLikelihood b : treeLikelihood) {
            double[] d = b.getPatternLogLikelihoods();
            System.arraycopy(d, 0, patternLogLikelihoods, i, d.length);
            i += d.length;
        }
        return patternLogLikelihoods;
    } // getPatternLogLikelihoods

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return Collections.singletonList(scsDataInput.get().getID());
    } // getArguments

    @Override
    public void callVariants() {
        try {
            if (threadCount > 1) {
                pool.invokeAll(MLGenotypesCallers);
            } else {
                treeLikelihood[0].callVariants();
            }
        } catch (RejectedExecutionException | InterruptedException e) {
            e.printStackTrace();
            System.exit(0);
        }
    } // callVariants


    //***********************************************
    //*           CalculationNode methods           *
    //***********************************************

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        boolean requiresRecalculation = false;
        for (ScsTreeLikelihood b : treeLikelihood) {
            requiresRecalculation |= b.requiresRecalculation();
        }
        return requiresRecalculation;
    } // requiresRecalculation

    @Override
    public void store() {
        super.store();
    } // store

    @Override
    public void restore() {
        super.restore();
    } // restore


    //**********************************************
    //*               Nested classes               *
    //**********************************************

    class ScsTreeLikelihoodCaller implements Callable<Double> {
        private final ScsTreeLikelihood likelihood;
        private final int threadNr;

        public ScsTreeLikelihoodCaller(ScsTreeLikelihood likelihood, int threadNr) {
            this.likelihood = likelihood;
            this.threadNr = threadNr;
        }

        public Double call() {
            try {
                logPByThread[threadNr] = likelihood.calculateLogPByThread();

                if (useAscBiasCorrection)
                    constSumByThread[threadNr] = likelihood.getBiasCorr();
            } catch (Exception e) {
                System.err.println("Something went wrong in thread " + threadNr);
                e.printStackTrace();
                System.exit(0);
            }
            return logPByThread[threadNr];
        }

    } // class ScsTreeLikelihoodCaller

    class MLGenotypeCaller implements Callable<Double> {
        private final ScsTreeLikelihood likelihood;
        private final int threadNr;

        public MLGenotypeCaller(ScsTreeLikelihood likelihood, int threadNr) {
            this.likelihood = likelihood;
            this.threadNr = threadNr;
        }

        @Override
        public Double call() {
            try {
                likelihood.callVariants();
            } catch (Exception e) {
                System.err.println("Something went wrong in thread " + threadNr);
                e.printStackTrace();
                System.exit(0);
            }

            return null;
        }
    } // class MLGenotypeCaller


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    /**
     * Is tracing maximum likelihood genotypes or not.
     *
     * @return yes or no
     */
    public boolean isTraceMLGenotypes() {
        return treeLikelihood[0].isTraceMLGenotypes();
    } // isTraceMLGenotypes

    /**
     * find in which likelihood core the provided pattern is being computed
     *
     * @param patterIndex_ which pattern is looking for?
     * @return the index of likelihood core
     */
    public int getNrOfLikelihood(int patterIndex_) {
        if (patterIndex_ < patternPoints[0] || patterIndex_ >= patternPoints[patternPoints.length - 1]) {
            throw new IllegalArgumentException("The pattern index being looked for (" + patterIndex_ + ") is out " +
                    "of bound [" + patternPoints[0] + ", " + patternPoints[patternPoints.length - 1] + "). ");
        } else {
            for (int i = 0; i < patternPoints.length - 1; i++) {
                if (patterIndex_ >= patternPoints[i] && patterIndex_ < patternPoints[i + 1]) {
                    return i;
                }
            }

            return -1;
        }
    }

    /**
     * for the purpose of debug
     * must be overridden
     *
     * @param inDebugMode in debug mode or not
     */
    @Override
    public void updateInDebugMode(final boolean inDebugMode) {
        for (ScsTreeLikelihood t : treeLikelihood) {
            t.updateInDebugMode(inDebugMode);
        }
    } // updateInDebugMode

    @Override
    public String getSortedCellNamesFromVarInfo(String separator) {
        return treeLikelihood[0].getSortedCellNamesFromVarInfo(separator);
    } // getSortedCellNamesFromVarInfo

    /**
     * initialize header of logger
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    @Override
    public void initSiteHeader(int start, PrintStream out) {
        for (int i = 0; i != treeLikelihood.length; i++) {
            treeLikelihood[i].initSiteHeader(start, out);

            if (i < treeLikelihood.length - 1)
                out.print("\t");

            start += treeLikelihood[i].scsDataInput.get().getLociNr();
        }
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
        for (ScsTreeLikelihood i : treeLikelihood) {
            i.closeGenotypeAndAdo(start, out);
            start += i.scsDataInput.get().getLociNr();
        }
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
        for (int i = 0; i != treeLikelihood.length; i++) {
            treeLikelihood[i].logCovar(out);

            if (i < treeLikelihood.length - 1)
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
        for (ScsTreeLikelihood i : treeLikelihood) {
            i.closeCovar(start, out);
            start += i.scsDataInput.get().getLociNr();
        }
    } // closeCovar

}
