package beast.math.util;

import beast.app.variantcaller.EstimatesTypeCollection;
import beast.math.statistic.DiscreteStatistics;
import org.apache.commons.lang3.ArrayUtils;
import org.jetbrains.annotations.NotNull;
import smile.stat.distribution.KernelDensity;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class MathFunctions {


    private MathFunctions() {}


    //**********************************************
    //*              Static Variables              *
    //**********************************************

    /**
     * Scale a log likelihood to phred-scaled score by multiplication
     * log10(number) = log(number) * log10(e)
     */
    public static final double PHRED_SCORE_COEFFICIENT = -10 * Math.log10(Math.E);


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static double[] expArr(double @NotNull [] x) {
        double[] y = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            y[i] = Math.exp(x[i]);
        }

        return y;
    } // expArr

    /**
     * compute mean (ruling out NaN, positive infinity and negative infinity)
     *
     * @param x           list of numbers
     * @param includeZero whether to include zero or not
     * @return mean
     */
    public static double mean(double[] x, boolean includeZero) {
        double m = 0;
        int count = x.length;
        for (double aX : x) {
            if (Double.isNaN(aX) || Double.isInfinite(aX) || (!includeZero && aX == 0)) {
                count--;
            } else {
                m += aX;
            }
        }

        return m / count;
    } // mean

    /**
     * compute the geometric mean.
     *
     * @param arr   data
     * @param start start position of arr
     * @param end   end position of arr (exclusive)
     * @return geometric mean
     */
    public static double geometricMean(
            final int[] arr,
            final int start,
            final int end
    ) {
        if (start >= end || start >= arr.length || end > arr.length)
            throw new IllegalArgumentException("Error: Illegal bounds for computing geometric mean.");

        double sum = 0.0;
        int len = 0;

        for (int i = start; i < end; i++) {
            if (arr[i] > 0) {
                ++len;
                sum += Math.log(arr[i]);
            }
        }

        return len == 0 ? 0.0 : Math.exp(sum / len);
    } // geometricMean

    /**
     * compute the logarithm of sum
     * https://statmodeling.stat.columbia.edu/2016/06/11/log-sum-of-exponentials/
     *
     * @param logElements logarithm elements
     * @return obviously
     */
    public static double logSumExp(final double[] logElements) {
        double result = 0.0;

        // find the maximum element
        double max = DiscreteStatistics.max(logElements);

        for (double i : logElements) {
            result += Math.exp(i - max);
        }

        return Math.log(result) + max;
    } // logSumExp

    public static int sum(final int[] values) {
        return Arrays.stream(values).reduce(0, Integer::sum);
    } // sum

    public static double sum(final double[] values) {
        return Arrays.stream(values).reduce(0, Double::sum);
    } // sum

    public static <T extends Comparable<T>> T min(final T[] values) {
        T min = values[0];
        for (int i = 1; i < values.length; i++) {
            if (values[i].compareTo(min) < 0) min = values[i];
        }
        return min;
    } // min

    public static <T extends Comparable<T>> T max(final T[] values) {
        T max = values[0];
        for (int i = 1; i < values.length; i++) {
            if (values[i].compareTo(max) > 0) max = values[i];
        }
        return max;
    } // max

    public static <T extends Comparable<T>> int minIndex(final T[] values) {
        int id = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i].compareTo(values[id]) < 0) id = i;
        }
        return id;
    } // minIndex

    public static <T extends Comparable<T>> List<Integer> minIndices(final T[] values) {
        T min = values[0];
        List<Integer> id = Stream.of(0).collect(Collectors.toList());
        for (int i = 1; i < values.length; i++) {
            if (values[i].compareTo(min) < 0) {
                min = values[i];
                id.clear();
                id.add(i);
            } else if (values[i].compareTo(min) == 0)
                id.add(i);
        }
        return id;
    } // minIndices

    public static <T extends Comparable<T>> int maxIndex(final T[] values) {
        int id = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i].compareTo(values[id]) > 0) id = i;
        }
        return id;
    } // maxIndex

    public static <T extends Comparable<T>> List<Integer> maxIndices(final T[] values) {
        T max = values[0];
        List<Integer> id = Stream.of(0).collect(Collectors.toList());
        for (int i = 1; i < values.length; i++) {
            if (values[i].compareTo(max) > 0) {
                max = values[i];
                id.clear();
                id.add(i);
            } else if (values[i].compareTo(max) == 0)
                id.add(i);
        }
        return id;
    } // maxIndices

    /**
     * get estimates
     * <p>
     * output added in the following order:
     * 0 - mean
     * 1 - median
     * mode (MAP), different kernel density estimate function:
     * 2 - "gaussian"
     *
     * @param estimatesType    which kind of estimates to get?
     * @param modeKDEType      which kind of KDE distribution is used for mode estimates?
     * @param burninPercentage apparently
     * @param samples          input
     * @param samplesEstimates output
     */
    public static void getEstimates(EstimatesTypeCollection.EstimatesType estimatesType,
                                    EstimatesTypeCollection.ModeKDEType modeKDEType,
                                    int burninPercentage,
                                    Map<String, List<Double>> samples,
                                    Map<String, Map<String, Double>> samplesEstimates) throws NullPointerException {
        for (String key : samples.keySet()) {
            List<Double> origin = samples.get(key);
            Map<String, Double> estimates = new HashMap<>();

            int burninCount = Math.max(0, (burninPercentage * origin.size()) / 100);

            double[] filtered = ArrayUtils.toPrimitive(origin.subList(burninCount, origin.size()).toArray(new Double[0]));

            getEstimates(estimatesType, modeKDEType, filtered, estimates);

            samplesEstimates.put(key, estimates);
        }
    } // getEstimates

    /**
     * get estimates
     * <p>
     * output added in the following order:
     * 0 - mean
     * 1 - median
     * 2 - mode (MAP), using "gaussian" kernel density estimate function
     *
     * @param estimatesType    which kind of estimates to get?
     * @param burninPercentage apparently
     * @param samples          input
     * @param samplesEstimates output
     */
    public static void getArrEstimates(
            EstimatesTypeCollection.EstimatesType estimatesType,
            int burninPercentage,
            Map<String, List<Double>[][]> samples,
            Map<String, Map<String, double[][]>> samplesEstimates
    ) throws NullPointerException {
        int i, j;
        for (String key : samples.keySet()) {
            i = 0;

            Map<String, double[][]> estimates = new HashMap<>();

            for (List<Double>[] matrix : samples.get(key)) {
                j = 0;

                for (List<Double> origin : matrix) {
                    int burninCount = Math.max(0, (burninPercentage * origin.size()) / 100);

                    double[] filtered = ArrayUtils.toPrimitive(origin.subList(burninCount, origin.size()).toArray(new Double[0]));

                    getArrEstimates(estimatesType, filtered, estimates,
                            samples.get(key).length, i, matrix.length, j);

                    j++;
                }

                i++;
            }

            samplesEstimates.put(key, estimates);
        }
    } // getArrEstimates

    /**
     * output added in the following order:
     * 0 - mean
     * 1 - median
     * mode (MAP), different kernel density estimate function:
     * 2 - "gaussian"
     *
     * @param estimatesType which kind of estimates to get?
     * @param modeKDEType   which kind of KDE distribution is used for mode estimates?
     * @param data          input
     * @param estimates     output
     */
    public static void getEstimates(EstimatesTypeCollection.EstimatesType estimatesType,
                                    EstimatesTypeCollection.ModeKDEType modeKDEType,
                                    final double[] data,
                                    Map<String, Double> estimates) throws NullPointerException {
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE && modeKDEType == null) {
            throw new NullPointerException("Unspecified KDE distribution for mode estimates.");
        }

        if (estimates == null) {
            estimates = new HashMap<>();
        }

        // mean
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MEAN ||
                estimatesType == EstimatesTypeCollection.EstimatesType.ALL) {
            estimates.put(EstimatesTypeCollection.EstimatesType.MEAN.toString().toLowerCase(), DiscreteStatistics.mean(data));
        }

        // median
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MEDIAN ||
                estimatesType == EstimatesTypeCollection.EstimatesType.ALL) {
            estimates.put(EstimatesTypeCollection.EstimatesType.MEDIAN.toString().toLowerCase(), DiscreteStatistics.median(data));
        }

        // mode with different KDE distributions
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE) {
            estimates.put(modeKDEType.toString().toLowerCase(), getModeEstimates(data));
        }

        if (estimatesType == EstimatesTypeCollection.EstimatesType.ALL) {
            for (EstimatesTypeCollection.ModeKDEType type : EstimatesTypeCollection.ModeKDEType.values()) {
                estimates.put(type.toString().toLowerCase(), getModeEstimates(data));
            }
        }
    } // getEstimates

    /**
     * output added in the following order:
     * 0 - mean
     * 1 - median
     * 2 - mode (MAP), using "gaussian" kernel density estimate function
     *
     * @param estimatesType which kind of estimates to get?
     * @param data          input
     * @param estimates     output
     * @param numOfMatrix   the number of matrices
     * @param index1        which matrix it is?
     * @param len           number of estimates for a key
     * @param index2        allelic sequencing coverage or raw variance?
     */
    public static void getArrEstimates(
            EstimatesTypeCollection.EstimatesType estimatesType,
            final double[] data,
            Map<String, double[][]> estimates,
            int numOfMatrix,
            int index1,
            int len,
            int index2
    ) {
        if (estimates == null) {
            estimates = new HashMap<>();
        }

        // mean
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MEAN ||
                estimatesType == EstimatesTypeCollection.EstimatesType.ALL) {
            estimates.computeIfAbsent(EstimatesTypeCollection.EstimatesType.MEAN.toString().toLowerCase(), k -> new double[numOfMatrix][len]);
            estimates.get(EstimatesTypeCollection.EstimatesType.MEAN.toString().toLowerCase())[index1][index2] = DiscreteStatistics.mean(data);
        }

        // median
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MEDIAN ||
                estimatesType == EstimatesTypeCollection.EstimatesType.ALL) {
            estimates.computeIfAbsent(EstimatesTypeCollection.EstimatesType.MEDIAN.toString().toLowerCase(), k -> new double[numOfMatrix][len]);
            estimates.get(EstimatesTypeCollection.EstimatesType.MEDIAN.toString().toLowerCase())[index1][index2] = DiscreteStatistics.median(data);
        }

        // mode
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE ||
                estimatesType == EstimatesTypeCollection.EstimatesType.ALL) {
            getModeArrEstimates(data, estimates, numOfMatrix, index1, len, index2);
        }
    } // getArrEstimates

    /**
     * get estimates for mode
     *
     * @param data input
     * @return mode value
     */
    public static double getModeEstimates(final double[] data) {
        KernelDensity kde = new KernelDensity(data);

        double[] probabilities = new double[data.length];
        double max = 0;
        int maxIndex = 0;

        for (int i = 0; i < data.length; i++) {
            probabilities[i] = kde.p(data[i]);
            if (probabilities[i] > max) {
                max = probabilities[i];
                maxIndex = i;
            }
        }

        return data[maxIndex];
    } // getModeEstimates

    /**
     * get estimates for mode
     *
     * @param data        input
     * @param estimates   output
     * @param len         number of estimates for a key
     * @param numOfMatrix the number of matrices
     * @param index1      which matrix it is?
     * @param len         number of estimates for a key
     * @param index2      allelic sequencing coverage or raw variance?
     */
    public static void getModeArrEstimates(
            final double[] data,
            Map<String, double[][]> estimates,
            int numOfMatrix,
            int index1,
            int len,
            int index2) {
        estimates.computeIfAbsent(EstimatesTypeCollection.ModeKDEType.GAUSSIAN.toString().toLowerCase(), k -> new double[numOfMatrix][len]);
        estimates.get(EstimatesTypeCollection.ModeKDEType.GAUSSIAN.toString().toLowerCase())[index1][index2] = getModeEstimates(data);
    } // getModeArrEstimates

    /**
     * Convert a logE value to a rounded phred scaled value.
     *
     * @param value log value
     * @return phread scaled value
     */
    public static long convertLogE2RoundedPhredScaled(double value) {
        return Math.round(PHRED_SCORE_COEFFICIENT * value);
    } // convertLogE2RoundedPhredScaled

    /**
     * Convert a logE value to an unrounded phred scaled value.
     *
     * @param value log value
     * @return phread scaled value
     */
    public static double convertLogE2PhredScaled(double value) {
        return PHRED_SCORE_COEFFICIENT * value;
    } // convertLogE2PhredScaled

}
