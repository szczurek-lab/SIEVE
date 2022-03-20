package beast.evolution.tree.coalescent;

import beast.core.Description;
import beast.core.Input;
import beast.math.Binomial;

@Description("Coalescent tree prior for tumor trees with a trunk and a MRCA root.")
public class ScsCoalescent extends Coalescent {

    final public Input<Double> expectedNrOfMuInput = new Input<>("expectedNrOfMu",
            "the expected number of accumulated mutations per site for a healthy cell becoming the tumor " +
                    "founding cell (1 over lambda)");

    protected double expectedNrOfMu;

    @Override
    public void initAndValidate() {
        if (expectedNrOfMuInput.get() != null) {
            if (expectedNrOfMuInput.get() <= 0) {
                throw new IllegalArgumentException("the expected number of accumulated mutations per site must be " +
                        "positive, but " + expectedNrOfMuInput.get() + "is found (" + this.getClass().getName() + ")");
            } else {
                expectedNrOfMu = expectedNrOfMuInput.get();
            }
        }

        super.initAndValidate();
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a population size function.
     * <p>
     * The tree is assumed having a trunk, connecting TMRCA and MRCA
     *
     * @param intervals       the intervals whose likelihood is computed
     * @param popSizeFunction the population size function
     * @param threshold       the minimum allowable coalescent interval size; negative infinity will be returned if
     *                        any non-zero intervals are smaller than this
     * @return the log likelihood of the intervals given the population size function
     */
    public double calculateLogLikelihood(IntervalList intervals, PopulationFunction popSizeFunction, double threshold) {

        double logL = 0.0;

        double startTime = 0.0;
        final int n = intervals.getIntervalCount();
        for (int i = 0; i < n; i++) {

            final double duration = intervals.getInterval(i);
            final double finishTime = startTime + duration;

            final double intervalArea = popSizeFunction.getIntegral(startTime, finishTime);
            if (intervalArea == 0 && duration > 1e-10) {
                /* the above test used to be duration != 0, but that leads to numerical issues on resume
                 * (https://github.com/CompEvol/beast2/issues/329) */
                return Double.NEGATIVE_INFINITY;
            }
            final int lineageCount = intervals.getLineageCount(i);

            final double kChoose2 = Binomial.choose2(lineageCount);
            // common part
            logL += -kChoose2 * intervalArea;

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                final double demographicAtCoalPoint = popSizeFunction.getPopSize(finishTime);

                // if value at end is many orders of magnitude different than mean over interval reject the interval
                // This is protection against cases where ridiculous infinitesimal
                // population size at the end of a linear interval drive coalescent values to infinity.

                if (duration == 0.0 || demographicAtCoalPoint * (intervalArea / duration) >= threshold) {
                    //                if( duration == 0.0 || demographicAtCoalPoint >= threshold * (duration/intervalArea) ) {
                    logL -= Math.log(demographicAtCoalPoint);
                } else {
                    // remove this at some stage
                    //  System.err.println("Warning: " + i + " " + demographicAtCoalPoint + " " + (intervalArea/duration) );
                    return Double.NEGATIVE_INFINITY;
                }
            } else if (intervals.getIntervalType(i) == IntervalType.NOTHING && i == n - 1) {
                if (expectedNrOfMuInput.get() != null) {
                    /* reach the root of the tree
                     * as the tree root has the largest height, therefore it must be at the end of the intervals.
                     * assuming it complies with exponential distribution
                     */
                    logL += Math.log(1 / expectedNrOfMu) - (finishTime - startTime) / expectedNrOfMu;
                }
            }
            startTime = finishTime;
        }

        return logL;
    } // calculateLogLikelihood

}
