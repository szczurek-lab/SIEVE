package beast.math.distributions;

import beast.core.Description;

import static org.apache.commons.math.special.Gamma.logGamma;

@Description("parameterized as modelling the number of failures before reaching a given number of successes with " +
        "a given successful probability")
public class NegativeBinomial extends ScsParametricDistribution {

    private double nrOfSuccess;
    private double propOfSuccess;

    /**
     * Return the probability density for a particular point.
     *
     * @param x             The point at which the density should be computed.
     * @param nrOfSuccess   The number of successful events.
     * @param propOfSuccess The probability of a successful event.
     * @return The pdf at point x.
     */
    public double density(double x, double nrOfSuccess, double propOfSuccess) {
        if (nrOfSuccess < 0 || propOfSuccess < 0) {
            throw new IllegalArgumentException("Both nrOfSuccess and propOfSuccess should be positive (" +
                    this.getClass().getName() + ")");
        }

        setNrOfSuccess(nrOfSuccess);
        setPropOfSuccess(propOfSuccess);
        return density(x);
    } // density

    /**
     * Return the log probability density for a particular point.
     *
     * @param x             The point at which the density should be computed.
     * @param nrOfSuccess   The number of successful events.
     * @param propOfSuccess The probability of a successful event.
     * @return The pdf at point x.
     */
    public double logDensity(double x, double nrOfSuccess, double propOfSuccess) {
        if (nrOfSuccess < 0 || propOfSuccess < 0) {
            throw new IllegalArgumentException("Both nrOfSuccess and propOfSuccess should be positive (" +
                    this.getClass().getName() + ")");
        }

        setNrOfSuccess(nrOfSuccess);
        setPropOfSuccess(propOfSuccess);
        return logDensity(x);
    } // logDensity

    /**
     * Return the probability density for a particular point.
     * refer to https://www.johndcook.com/negative_binomial.pdf for gamma parameterization
     * gamma distribution is parameterized with shape parameter 'alpha' and scale parameter 'beta'
     * alpha = nrOfSuccess
     * beta = 1 / propOfSuccess - 1
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    public double density(double x) {
        return Math.exp(logDensity(x));
    } // density

    /**
     * Return the probability density for a particular point.
     * refer to https://www.johndcook.com/negative_binomial.pdf for gamma parameterization
     * gamma distribution is parameterized with shape parameter 'alpha' and scale parameter 'beta'
     * alpha = nrOfSuccess
     * beta = 1.0 / propOfSuccess - 1.0
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    public double logDensity(double x) {
        double beta = 1.0 / propOfSuccess - 1.0;
        return logGamma(nrOfSuccess + x) + (nrOfSuccess + x) * Math.log(1.0 - propOfSuccess) - logGamma(x + 1.0) - logGamma(nrOfSuccess) - nrOfSuccess * Math.log(beta);
    } // logDensity


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    public double getNrOfSuccess() {
        return nrOfSuccess;
    }

    public void setNrOfSuccess(double nrOfSuccess) {
        this.nrOfSuccess = nrOfSuccess;
    }

    public double getPropOfSuccess() {
        return propOfSuccess;
    }

    public void setPropOfSuccess(double propOfSuccess) {
        this.propOfSuccess = propOfSuccess;
    }

}
