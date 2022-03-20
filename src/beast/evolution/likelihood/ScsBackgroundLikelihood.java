package beast.evolution.likelihood;


import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.rawreadcountsmodel.RawReadCountsModelInterface;
import beast.math.distributions.BetaBinomial;
import beast.math.distributions.DirichletMultinomial;

import java.util.List;
import java.util.Random;

import static org.apache.commons.math.special.Gamma.logGamma;

@Description("Compute the likelihood of background sites with presumed homozygous reference genotype (0/0) for " +
        "whole genome sequencing data or whole exome sequencing data.")
public class ScsBackgroundLikelihood extends Distribution {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    final public Input<ScsAlignment> scsDataInput = new Input<>("scsData", "single-cell sequence data " +
            "for the beast.tree", Input.Validate.REQUIRED);

    public final Input<RawReadCountsModelInterface.Base> rawReadCountsModelInput = new Input<>("rawReadCountsModel",
            "an error model applied to leaves of a tree, accounting for the sequencing error",
            Input.Validate.REQUIRED);

    public final Input<Boolean> runAnalysisInput = new Input<>("runAnalysis", "for debug purpose",
            Input.Validate.OPTIONAL);

    /**
     * the parameters of wild type genotype situation
     */
    protected double[] wildTypeNucReadCountsModelParams;

    protected int nrOfTaxa;
    protected long nrOfBackgroundSites;

    protected List<List<long[]>> backgroundInfo;

    protected RawReadCountsModelInterface.Base rawReadCountsModel;

    private boolean runAnalysis;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        nrOfTaxa = scsDataInput.get().getTaxonCount();
        nrOfBackgroundSites = scsDataInput.get().getNrOfBackgroundSites();
        backgroundInfo = scsDataInput.get().getBackgroundInfo();

        rawReadCountsModel = rawReadCountsModelInput.get();

        wildTypeNucReadCountsModelParams = new double[rawReadCountsModel.getNrOfParamsInReadCountsModel()];

        if (runAnalysisInput.get() == null)
            runAnalysis = false;
        else
            runAnalysis = runAnalysisInput.get();
    } // initAndValidate


    //**********************************************
    //*            Distribution methods            *
    //**********************************************

    /**
     * compute the log likelihood of background information
     * alpha <-> delta, beta <-> tau
     *
     * @return log likelihood of background information
     */
    @Override
    public double calculateLogP() {
        // get the latest parameters of beta-binomial or Dirichlet-multinomial distribution for situation 1
        rawReadCountsModel.getWildTypeNucReadCountModelParams(this.wildTypeNucReadCountsModelParams);

        if (this.rawReadCountsModel.getNucReadCountsModelDistribution() instanceof BetaBinomial)
            calculateLogPBetaBinomial();
        else if (this.rawReadCountsModel.getNucReadCountsModelDistribution() instanceof DirichletMultinomial)
            calculateLogPDirichletMultinomial();
        else
            throw new IllegalArgumentException("Error! Cannot process such a distribution: " +
                    this.rawReadCountsModel.getNucReadCountsModelDistribution().getClass().getName());

        return this.logP;
    } // calculateLogP

    private void calculateLogPBetaBinomial() {
        // note that in the computation, combinations with reads being 0 are included.

        assert this.backgroundInfo != null;

        double tmp = 0.0;
        this.logP = 0.0;

        final long startTime = System.currentTimeMillis();
        if (this.runAnalysis)
            System.out.println("Background likelihood computation:");

        /*
         * first, compute:
         * nrOfBackgroundSites * nrOfTaxa * (logGamma(alpha + beta) - logGamma(alpha) - logGamma(beta))
         * which only depends on the parameter in the error model (beta-binomial distribution)
         */
        this.logP = this.nrOfBackgroundSites * this.nrOfTaxa * (logGamma(this.wildTypeNucReadCountsModelParams[0] + this.wildTypeNucReadCountsModelParams[1]) - logGamma(this.wildTypeNucReadCountsModelParams[0]) - logGamma(this.wildTypeNucReadCountsModelParams[1]));
        if (runAnalysis) {
            tmp = this.logP;
            System.out.println("  Part 1 (likelihood of error model parameters): " + String.format("%.5f", tmp));
        }

        /*
         * second, compute:
         * sum(m) (N_m * logGamma(m + alpha))
         * regarding variant
         */
        for (long[] var : this.backgroundInfo.get(0)) {
            this.logP += var[1] * logGamma(var[0] + this.wildTypeNucReadCountsModelParams[0]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 2 (variant likelihood): " + String.format("%.5f", tmp));
        }

        /*
         * third, compute:
         * sum(c - m) (N_{c - m} * logGamma(c - m + beta))
         * regarding normal
         */
        for (long[] nor : backgroundInfo.get(1)) {
            this.logP += nor[1] * logGamma(nor[0] + this.wildTypeNucReadCountsModelParams[1]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 3 (normal likelihood): " + String.format("%.5f", tmp));
        }

        /*
         * fourth, compute:
         * - sum(c) (N_{c} * logGamma(c + alpha + beta))
         * regarding coverage
         */
        for (long[] cov : backgroundInfo.get(2)) {
            this.logP -= cov[1] * logGamma(cov[0] + this.wildTypeNucReadCountsModelParams[0] + this.wildTypeNucReadCountsModelParams[1]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 4 (coverage likelihood): " + String.format("%.5f", tmp));
            System.out.println("  Final likelihood: " + String.format("%.5f", this.logP));
        }

        final long endTime = System.currentTimeMillis();

        if (this.runAnalysis)
            System.out.println("  Run time: " + (endTime - startTime) + " milliseconds.\n");
    } // calculateLogPBetaBinom

    private void calculateLogPDirichletMultinomial() {
        // note that in the computation, combinations with reads being 0 are excluded.

        assert this.backgroundInfo != null;

        double tmp = 0.0;
        this.logP = 0.0;

        final long startTime = System.currentTimeMillis();
        if (this.runAnalysis)
            System.out.println("Background likelihood computation:");

        /*
         * first, compute:
         * nrOfBackgroundSites * nrOfTaxa * (logGamma(alpha_0) - logGamma(alpha_1) - logGamma(alpha_2) - logGamma(alpha_3) - logGamma(alpha_4))
         * which only depends on the parameter in the model of nucleotide read counts (Dirichlet-multinomial distribution)
         */
        this.logP = this.nrOfBackgroundSites * this.nrOfTaxa * (
                logGamma(this.wildTypeNucReadCountsModelParams[4])
                        - logGamma(this.wildTypeNucReadCountsModelParams[0])
                        - logGamma(this.wildTypeNucReadCountsModelParams[1])
                        - logGamma(this.wildTypeNucReadCountsModelParams[2])
                        - logGamma(this.wildTypeNucReadCountsModelParams[3])
        );
        if (this.runAnalysis) {
            tmp = this.logP;
            System.out.println("  Part 1 (likelihood of error model parameters): " + String.format("%.5f", tmp));
        }

        /*
         * second, compute (occurrence of 0 will be excluded):
         * - sum(c) (N_{c} * logGamma(alpha_0 + c))
         * regarding coverage
         */
        for (final long[] item : this.backgroundInfo.get(4)) {
            this.logP -= item[1] * logGamma(this.wildTypeNucReadCountsModelParams[4] + item[0]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 2 (coverage likelihood): " + String.format("%.5f", tmp));
        }

        /*
         * third, compute:
         * sum(m_1) (N_{m_1} * logGamma(alpha_1 + m_1))
         * regarding variant1
         */
        for (final long[] item : this.backgroundInfo.get(0)) {
            this.logP += item[1] * logGamma(this.wildTypeNucReadCountsModelParams[0] + item[0]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 3 (variant1 likelihood): " + String.format("%.5f", tmp));
        }

        /*
         * fourth, compute:
         * sum(m_2) (N_{m_2} * logGamma(alpha_2 + m_2))
         * regarding variant2
         */
        for (final long[] item : this.backgroundInfo.get(1)) {
            this.logP += item[1] * logGamma(this.wildTypeNucReadCountsModelParams[1] + item[0]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 4 (variant2 likelihood): " + String.format("%.5f", tmp));
        }

        /*
         * fifth, compute:
         * sum(m_3) (N_{m_3} * logGamma(alpha_3 + m_3))
         * regarding variant3
         */
        for (final long[] item : this.backgroundInfo.get(2)) {
            this.logP += item[1] * logGamma(this.wildTypeNucReadCountsModelParams[2] + item[0]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 5 (variant3 likelihood): " + String.format("%.5f", tmp));
        }

        /*
         * sixth, compute:
         * sum(m_4) N_{m_4} * (logGamma(alpha_4 + m_4))
         * regarding normal
         */
        for (final long[] item : this.backgroundInfo.get(3)) {
            this.logP += item[1] * logGamma(this.wildTypeNucReadCountsModelParams[3] + item[0]);
        }
        if (this.runAnalysis) {
            tmp = this.logP - tmp;
            System.out.println("  Part 6 (normal likelihood): " + String.format("%.5f", tmp));
            System.out.println("  Final likelihood: " + String.format("%.5f", this.logP));
        }

        final long endTime = System.currentTimeMillis();

        if (this.runAnalysis)
            System.out.println("  Run time: " + (endTime - startTime) + " milliseconds.\n");
    } // calculateLogPDirichletMultinomial

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

}
