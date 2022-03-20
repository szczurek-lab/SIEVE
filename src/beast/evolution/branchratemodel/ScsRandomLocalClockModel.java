package beast.evolution.branchratemodel;

import beast.core.Citation;
import beast.core.Description;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;

@Description("Enhanced random Local Clock Model.")
@Citation(value =
        "Drummond AJ, Suchard MA (2010) Bayesian random local clocks, or one rate to rule them all. BMC biology 8, 114.",
        DOI = "10.1186/1741-7007-8-114",
        year = 2010, firstAuthorSurname = "drummond")
public class ScsRandomLocalClockModel extends RandomLocalClockModel {

    private boolean recompute = true;

    /**
     * This is a recursive function that does the work of
     * calculating the unscaled branch rates across the tree
     * taking into account the indicator variables.
     *
     * @param node the node
     * @param rate the rate of the parent node
     */
    private void calculateUnscaledBranchRates(Node node, double rate, BooleanParameter indicators, RealParameter rates) {

        int nodeNumber = getNr(node);

        if (!node.isRoot()) {
            if (indicators.getValue(nodeNumber)) {
                if (ratesAreMultipliers) {
                    rate *= rates.getValue(nodeNumber);
                } else {
                    rate = rates.getValue(nodeNumber);
                }
            }
        }
        unscaledBranchRates[nodeNumber] = rate;

        if (!node.isLeaf()) {
            for (Node i : node.getChildren()) {
                calculateUnscaledBranchRates(i, rate, indicators, rates);
            }
        }
    }

    private void recalculateScaleFactor() {
        BooleanParameter indicators = indicatorParamInput.get();
        RealParameter rates = rateParamInput.get();

        double rootRate = 1.0;
        if (includeRootInput.get()) rootRate = rates.getValue(tree.getRoot().getNr());

        calculateUnscaledBranchRates(tree.getRoot(), rootRate, indicators, rates);

        if (scaling) {

            double timeTotal = 0.0;
            double branchTotal = 0.0;

            for (int i = 0; i < tree.getNodeCount(); i++) {
                Node node = tree.getNode(i);
                if (!node.isRoot()) {

                    double branchInTime = node.getParent().getHeight() - node.getHeight();

                    double branchLength = branchInTime * unscaledBranchRates[node.getNr()];

                    timeTotal += branchInTime;
                    branchTotal += branchLength;
                }
            }

            scaleFactor = timeTotal / branchTotal;

            scaleFactor *= meanRate.getValue();
        } else {
            scaleFactor = 1.0;
        }
    }

    @Override
    public double getRateForBranch(Node node) {
        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads
        synchronized (this) {
            if (recompute) {
                recalculateScaleFactor();
                recompute = false;
            }
        }

        return unscaledBranchRates[getNr(node)] * scaleFactor;
    }

    private int getNr(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr > tree.getRoot().getNr()) {
            nodeNr--;
        }
        return nodeNr;
    }

    public boolean isScaling() {
        return scaling;
    }

}
