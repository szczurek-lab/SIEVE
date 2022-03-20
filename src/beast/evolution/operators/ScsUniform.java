package beast.evolution.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.evolution.tree.ScsTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Randomly selects true internal scsTree node (i.e. not the root) and move node height uniformly in interval " +
        "restricted by the nodes parent and children.")
public class ScsUniform extends Uniform {

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    @SuppressWarnings("deprecation")
    public double proposal() {
        // a scsTree instance should be the input; otherwise use Uniform operator
        final Tree tree = treeInput.get(this);
        assert tree instanceof ScsTree; // input tree should be an instance of ScsTree

        // randomly select internal node
        final int nodeCount = tree.getNodeCount();

        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount() == 1)
            return Double.NEGATIVE_INFINITY;

        Node node;
        do {
            final int nodeNr = nodeCount / 2 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf());
        final double upper = node.getParent().getHeight();
        final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
        node.setHeight(newValue);

        ((ScsTree) tree).updateRootTMRCA();

        return 0.0;
    }

}
