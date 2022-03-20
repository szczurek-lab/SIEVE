package beast.math.distributions;

import static org.apache.commons.math.special.Beta.logBeta;

public class DirichletMultinomial extends ScsParametricDistribution {

    /**
     * Compute probability density of Dirichlet-multinomial distribution.
     *
     * @param counts    An array of counts, where the last element is the total count.
     * @param fractions An array of fractions matching each count in {@param counts}, respectively.
     * @return probability density
     */
    public double density(int[] counts, double[] fractions) {
        return Math.exp(logDensity(counts, fractions));
    } // density

    /**
     * Compute log probability density of Dirichlet-multinomial distribution.
     *
     * @param counts    An array of counts, where the last element is the total count.
     * @param fractions An array of fractions matching each count in {@param counts}, respectively.
     * @return probability density
     */
    public double logDensity(int[] counts, double[] fractions) {
        if (counts.length != fractions.length)
            throw new IllegalArgumentException("Error! Counts should have the same length as fractions does.");

        if (counts[counts.length - 1] == 0)
            return 1.0;

        double logValue = 0;
        int i = 0;
        for (; i < counts.length - 1; i++)
            logValue += logPartial(counts[i], fractions[i]);

        return logPartial(counts[i], fractions[i]) - logValue;
    } // logDensity

    /**
     * Helper method to compute part of the probability.
     * Logarithm of : {@param count} * Beta({@param fraction}, {@param count})
     *
     * @param count    apparently
     * @param fraction apparently
     * @return log({ @ param count }) + logBeta({@param fraction}, {@param count})
     */
    private double logPartial(double count, double fraction) {
        return count <= 0 ? 0 : Math.log(count) + logBeta(fraction, count);
    } // logPartial

    /**
     * Return the probability density for a particular point.
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    @Override
    public double density(double x) {
        throw new IllegalArgumentException("Error! This method is not implemented.");
    }

    @Override
    public double logDensity(double x) {
        throw new IllegalArgumentException("Error! This method is not implemented.");
    }

}
