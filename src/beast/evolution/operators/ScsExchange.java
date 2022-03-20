package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.ScsNode;
import beast.evolution.tree.ScsTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Modified based on Exchange move")
public class ScsExchange extends TreeOperator {

    final public Input<Boolean> isNarrowInput = new Input<>("isNarrow", "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);

    @Override
    public void initAndValidate() {
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        assert tree instanceof ScsTree; // input tree should be an instance of ScsTree

        double logHastingsRatio = 0;

        if (isNarrowInput.get()) {
            logHastingsRatio = narrow(tree);
        } else {
            logHastingsRatio = wide(tree);
        }

        return logHastingsRatio;
    }

    /**
     * whether an ancestor node is a valid grandpa?
     *
     * @param n node
     * @return valid or not
     */
    @SuppressWarnings("deprecation")
    private int isAncValidGP(final Node n) {
        return (n.isRoot() || (n.getLeft().isLeaf() && n.getRight().isLeaf())) ? 0 : 1;
    }

    /**
     * whether a node (leaf or ancestor) is a valid grandpa?
     *
     * @param n node
     * @return valid or not
     */
    private int isValidGP(final Node n) {
        return n.isLeaf() ? 0 : isAncValidGP(n);
    }

    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     */
    @SuppressWarnings("deprecation")
    public double narrow(final Tree tree) {
        final int internalNodes = tree.getInternalNodeCount();
        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        Node grandParent;
        do {
            grandParent = tree.getNode(internalNodes + Randomizer.nextInt(internalNodes));
        } while (grandParent.isRoot() || (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()));

        // one of the parent's children is about to exchange with its uncle
        Node parentIndex = grandParent.getLeft();
        Node uncle = grandParent.getRight();
        if (parentIndex.getHeight() < uncle.getHeight()) {
            parentIndex = grandParent.getRight();
            uncle = grandParent.getLeft();
        }

        if (parentIndex.isLeaf()) {
            // tree with dated tips
            return Double.NEGATIVE_INFINITY;
        }

        int validGP = 0;
        {
            for (int i = internalNodes; i < 2 * internalNodes; ++i) {
                validGP += isAncValidGP(tree.getNode(i));
            }
        }

        final int c2 = isValidGP(parentIndex) + isValidGP(uncle);

        final Node i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());
        exchangeNodes(i, uncle, parentIndex, grandParent);

        final int validGPafter = validGP - c2 + isValidGP(parentIndex) + isValidGP(uncle);

        ((ScsTree) tree).updateRootTMRCA();

        return Math.log((float) validGP / validGPafter);
    } // narrow

    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     *
     * @param tree
     */
    public double wide(final Tree tree) {
        final int nodeCount = tree.getNodeCount();

        Node i;
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot() || ((ScsNode) i).isTMRCA());

        Node j;
        do {
            j = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (j.getNr() == i.getNr() || j.isRoot() || ((ScsNode) j).isTMRCA());

        final Node p = i.getParent();
        final Node jP = j.getParent();

        if ((p != jP) && (i != jP) && (j != p)
                && (j.getHeight() < p.getHeight())
                && (i.getHeight() < jP.getHeight())) {
            exchangeNodes(i, j, p, jP);

            // All the nodes on the path from i/j to the common ancestor of i/j parents had a topology change,
            // so they need to be marked FILTHY.
            if (markCladesInput.get()) {
                Node iup = p;
                Node jup = jP;
                while (iup != jup) {
                    if (iup.getHeight() < jup.getHeight()) {
                        assert !iup.isRoot() && !((ScsNode) iup).isTMRCA();
                        iup = iup.getParent();
                        iup.makeDirty(Tree.IS_FILTHY);
                    } else {
                        assert !jup.isRoot() && !((ScsNode) jup).isTMRCA();
                        jup = jup.getParent();
                        jup.makeDirty(Tree.IS_FILTHY);
                    }
                }
            }

            ((ScsTree) tree).updateRootTMRCA();

            return 0;
        }

        // Randomly selected nodes i and j are not valid candidates for a wide exchange.
        // reject instead of counting (like we do for narrow).
        return Double.NEGATIVE_INFINITY;
    } // wide

    /* exchange sub-trees whose root are i and j */

    protected void exchangeNodes(Node i, Node j,
                                 Node p, Node jP) {
        // precondition p -> i & jP -> j
        replace(p, i, j);
        replace(jP, j, i);
        // postcondition p -> j & p -> i
    }

}
