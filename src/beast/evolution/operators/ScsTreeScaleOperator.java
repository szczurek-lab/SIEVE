package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.ScsTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.text.DecimalFormat;

@Description("Modified based on ScaleOperator; should only be applied to ScsTree")
public class ScsTreeScaleOperator extends Operator {

    public final Input<Tree> treeInput = new Input<>("tree", "if specified, all beast.tree divergence times are scaled", Input.Validate.REQUIRED);

    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);

    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "scale root of a tree only, ignored if tree is not specified (default false)", false);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);

    private double m_fScaleFactor;
    private double upper, lower;

    @Override
    public void initAndValidate() {
        m_fScaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
    }

    protected double getScaler() {
        return (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
    }

    @Override
    @SuppressWarnings("deprecation")
    public double proposal() {
        final double scale = getScaler();
        final Tree tree = treeInput.get(this);
        assert tree instanceof ScsTree; // input tree should be an instance of ScsTree

        try {
            if (rootOnlyInput.get()) {
                final Node root = tree.getRoot();
                final double newHeight = root.getHeight() * scale;

                if (newHeight < Math.max(root.getLeft().getHeight(), (root.getRight() == null ? 0 : root.getRight().getHeight()))) {
                    return Double.NEGATIVE_INFINITY;
                }
                root.setHeight(newHeight);

                ((ScsTree) tree).updateRootTMRCA();

                return -Math.log(scale);
            } else {
                // scale the beast.tree
                final int internalNodes = tree.scale(scale);

                ((ScsTree) tree).updateRootTMRCA();

                // TODO: check the degree of freedom
                return Math.log(scale) * (internalNodes - 2);
            }
        } catch (Exception e) {
            return Double.NEGATIVE_INFINITY;
        }
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / m_fScaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return m_fScaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        m_fScaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(m_fScaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

}
