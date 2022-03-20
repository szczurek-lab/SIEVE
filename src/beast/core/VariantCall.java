package beast.core;

import beast.core.util.CompoundDistribution;
import beast.core.variantslogger.VariantsLogger;
import beast.evolution.likelihood.ScsGenericTreeLikelihood;

@Description("For the purpose of variant calling")
public class VariantCall extends Runnable {


    //**********************************************
    //*                   Inputs                   *
    //**********************************************

    final public Input<State> startStateInput =
            new Input<>("state", "elements of the state space", Input.Validate.REQUIRED);

    final public Input<Distribution> posteriorInput =
            new Input<>("distribution", "probability distribution to sample over (e.g. a posterior)",
                    Input.Validate.REQUIRED);

    final public Input<VariantsLogger> variantsLoggerInput =
            new Input<>("variantsLogger", "a logger for variant calling results", Input.Validate.REQUIRED);


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    /**
     * The state that takes care of managing StateNodes.
     */
    protected State state;

    protected double logPosterior;
    protected Distribution posterior;

    protected VariantsLogger variantsLogger;


    //***********************************************
    //*           Initialize and validate           *
    //***********************************************

    public void initAndValidate() {
        this.state = startStateInput.get();

        this.state.initialise();
        this.state.setPosterior(posteriorInput.get());

        this.variantsLogger = variantsLoggerInput.get();
    } // initAndValidate


    //*********************************************
    //*                  Methods                  *
    //*********************************************

    @Override
    public void run() {
        state.initAndValidate();
        state.setEverythingDirty(true);

        posterior = posteriorInput.get();

        logPosterior = state.robustlyCalcPosterior(posterior);

        state.storeCalculationNodes();
    } // run

    public double getLogPosterior() {
        return this.logPosterior;
    } // getLogPosterior

    public double getLogTreeLikelihood() {
        for (Distribution dist1 : ((CompoundDistribution) posterior).pDistributions.get()) {
            if (dist1 instanceof CompoundDistribution) {
                for (Distribution dist2 : ((CompoundDistribution) dist1).pDistributions.get()) {
                    if (dist2 instanceof ScsGenericTreeLikelihood) {
                        return dist2.getArrayValue();
                    }
                }
            }
        }

        return 0.0;
    } // getLogTreeLikelihood

    public void callVariantsAndGenerateVCF() {
        callVariants();
        generateVCF();
    } // callVariantsAndGenerateVCF

    protected void callVariants() {
        for (Distribution dist1 : ((CompoundDistribution) posterior).pDistributions.get()) {
            if (dist1 instanceof CompoundDistribution) {
                for (Distribution dist2 : ((CompoundDistribution) dist1).pDistributions.get()) {
                    if (dist2 instanceof ScsGenericTreeLikelihood) {
                        ((ScsGenericTreeLikelihood) dist2).callVariants();
                    }
                }
            }
        }
    } // callVariant

    protected void generateVCF() {
        variantsLogger.init();
        variantsLogger.log();
        variantsLogger.close();
    } // generateVCF

}
