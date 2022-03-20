package beast.core;

import beast.core.util.CompoundDistribution;
import beast.evolution.likelihood.ScsGenericTreeLikelihood;
import beast.evolution.likelihood.ScsTreeLikelihood;
import beast.evolution.likelihood.ThreadedScsTreeLikelihood;

@Description("An extension of State class")
public class ScsState extends State {

    @Override
    public void initAndValidate() {
    }

    @Override
    public double robustlyCalcPosterior(final Distribution posterior) {
        updateInDebugMode(posterior, true);
        store(-1);
        setEverythingDirty(true);
        checkCalculationNodesDirtiness();
        final double logLikelihood = posterior.calculateLogP();
        setEverythingDirty(false);
        acceptCalculationNodes();
        updateInDebugMode(posterior, false);

        return logLikelihood;
    } // robustlyCalcPosterior

    /**
     * note that what we are saving here is the most recent
     *
     * @param posterior
     */
    @Override
    public double robustlyCalcNonStochasticPosterior(Distribution posterior) {
        updateInDebugMode(posterior, true);

        store(-1);
        setEverythingDirty(true);
        storeCalculationNodes();
        checkCalculationNodesDirtiness();
        final double logLikelihood = posterior.getNonStochasticLogP();
        setEverythingDirty(false);
        acceptCalculationNodes();

        updateInDebugMode(posterior, false);

        return logLikelihood;
    } // robustlyCalcNonStochasticPosterior

    private void updateInDebugMode(final Distribution posterior, final boolean inDebugMode) {
        if (posterior instanceof CompoundDistribution) {

            for (Distribution dists1 : ((CompoundDistribution) posterior).pDistributions.get()) {
                if (dists1 instanceof CompoundDistribution) {
                    for (Distribution dists2 : ((CompoundDistribution) dists1).pDistributions.get()) {
                        if (dists2 instanceof ScsGenericTreeLikelihood) {
                            if (dists2 instanceof ScsTreeLikelihood) {
                                ((ScsTreeLikelihood) dists2).updateInDebugMode(inDebugMode);
                            } else if (dists2 instanceof ThreadedScsTreeLikelihood) {
                                ((ThreadedScsTreeLikelihood) dists2).updateInDebugMode(inDebugMode);
                            }
                        }
                    }
                }
            }
        }
    } // updateInDebugMode

}
