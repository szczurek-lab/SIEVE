package beast.math.distributions;

import beast.core.Description;

import static org.apache.commons.math.special.Gamma.logGamma;

@Description("implement beta-binomial distribution")
public class BetaBinomial extends ScsParametricDistribution {

    private long nrOfTrials;
    private double alpha;
    private double beta;

    /**
     * Return the probability density for a particular point.
     *
     * @param x          The point at which the density should be computed.
     * @param nrOfTrials The number of trials.
     * @param alpha      alpha parameter
     * @param beta       beta parameter
     * @return The pdf at point x.
     */
    public double density(double x, long nrOfTrials, double alpha, double beta) {
        if (nrOfTrials < 0 || alpha <= 0 || beta <= 0) {
            throw new IllegalArgumentException("nrOfTrials, alpha and beta should be positive (" +
                    this.getClass().getName() + ")");
        }

        setNrOfTrials(nrOfTrials);
        setAlpha(alpha);
        setBeta(beta);
        return density(x);
    }

    /**
     * Return the log probability density for a particular point.
     *
     * @param x          The point at which the density should be computed.
     * @param nrOfTrials The number of trials.
     * @param alpha      alpha parameter
     * @param beta       beta parameter
     * @return The pdf at point x.
     */
    public double logDensity(double x, int nrOfTrials, double alpha, double beta) {
        if (nrOfTrials < 0 || alpha <= 0 || beta <= 0) {
            throw new IllegalArgumentException("nrOfTrials, alpha and beta should be positive (" +
                    this.getClass().getName() + ")");
        }

        setNrOfTrials(nrOfTrials);
        setAlpha(alpha);
        setBeta(beta);
        return logDensity(x);
    }

    /**
     * Return the probability density for a particular point.
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    public double density(double x) {
        return Math.exp(logDensity(x));
    }

    /**
     * Return the log probability density for a particular point.
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    public double logDensity(double x) {
        return logGamma(nrOfTrials + 1) +
                logGamma(x + alpha) +
                logGamma(nrOfTrials - x + beta) +
                logGamma(alpha + beta) -
                logGamma(x + 1) -
                logGamma(nrOfTrials - x + 1) -
                logGamma(nrOfTrials + alpha + beta) -
                logGamma(alpha) -
                logGamma(beta);
    }


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    public long getNrOfTrials() {
        return nrOfTrials;
    }

    public void setNrOfTrials(long nrOfTrials) {
        this.nrOfTrials = nrOfTrials;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

}
