package beast.math.distributions;

import beast.core.CalculationNode;
import beast.core.Description;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;

@Description("")
public abstract class ScsParametricDistribution extends CalculationNode implements ContinuousDistribution {

    @Override
    public void initAndValidate() {

    }


    //***********************************************
    //*               Not implemented               *
    //***********************************************

    /**
     * For this distribution, X, this method returns x such that P(X &lt; x) = p.
     *
     * @param p the cumulative probability.
     * @return x.
     * @throws MathException if the inverse cumulative probability can not be
     *                       computed due to convergence or other numerical errors.
     */
    @Override
    @Deprecated
    public double inverseCumulativeProbability(double p) throws MathException {
        return 0;
    }

    /**
     * For a random variable X whose values are distributed according
     * to this distribution, this method returns P(X &le; x).  In other words,
     * this method represents the  (cumulative) distribution function, or
     * CDF, for this distribution.
     *
     * @param x the value at which the distribution function is evaluated.
     * @return the probability that a random variable with this
     * distribution takes a value less than or equal to <code>x</code>
     * @throws MathException if the cumulative probability can not be
     *                       computed due to convergence or other numerical errors.
     */
    @Override
    @Deprecated
    public double cumulativeProbability(double x) throws MathException {
        return 0;
    }

    /**
     * For a random variable X whose values are distributed according
     * to this distribution, this method returns P(x0 &le; X &le; x1).
     *
     * @param x0 the (inclusive) lower bound
     * @param x1 the (inclusive) upper bound
     * @return the probability that a random variable with this distribution
     * will take a value between <code>x0</code> and <code>x1</code>,
     * including the endpoints
     * @throws MathException            if the cumulative probability can not be
     *                                  computed due to convergence or other numerical errors.
     * @throws IllegalArgumentException if <code>x0 > x1</code>
     */
    @Override
    @Deprecated
    public double cumulativeProbability(double x0, double x1) throws MathException {
        return 0;
    }

}
