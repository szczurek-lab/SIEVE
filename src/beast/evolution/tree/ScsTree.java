package beast.evolution.tree;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.util.Log;
import beast.util.NexusParser;
import beast.util.NotSingleException;
import beast.util.ScsNexusParser;
import beast.util.ScsTreeParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

@Description("The tree built according to the finite-sites mutation, deletion, and insertion model using single-cell " +
        "DNA sequencing data has a root which rules the tumor subtree. In other words, the only difference between " +
        "our tree and the beast tree is that we have one more node with a trunk connecting the beast tree.")
public class ScsTree extends Tree {

    final static private String MEAN_RATE_NAME = "rate";
    final static private String MEDIAN_RATE_NAME = "rate_median";


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    final public Input<String> treeFileNameInput = new Input<>("treeFileName", "the full path of a file " +
            "containing a tree (nexus or newick) to initialize from");

    final public Input<Boolean> prioritizeMedianRateInput = new Input<>("prioritizeMedianRate",
            "prioritizing median rate value over mean rate value for non-strict molecular clock model",
            false, Input.Validate.OPTIONAL);

    /**
     * the most recent common ancestor of all sampled tumor cells
     * the only daughter of thr root of the tree (MRCA)
     */
    protected Node rootTMRCA;

    private boolean prioritizeMedianRate = false;

    /**
     * branch rates from the input tree file
     * [number of nodes]
     */
    private double[] branchRates;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public ScsTree() {

    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        if (m_initial.get() != null && !(this instanceof StateNodeInitialiser)) {
            throw new RuntimeException("initial-input should be specified for tree that is not a StateNodeInitialiser");
        }

        if (treeFileNameInput.get() != null) {
            try {
                final Tree other = readFromFile(treeFileNameInput.get());
                root = other.root.copy();
                rootTMRCA = root.getChild(0);
                nodeCount = other.nodeCount;
                internalNodeCount = other.internalNodeCount;
                leafNodeCount = other.leafNodeCount;
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }

            if (prioritizeMedianRateInput.get() != null)
                prioritizeMedianRate = prioritizeMedianRateInput.get();

            branchRates = new double[nodeCount];
            initialiseBranchRates(root);
        }

        if (nodeCount < 0) {
            if (m_taxonset.get() != null) {
                makeCaterpillar(0, 1, false);
            } else {
                throw new IllegalArgumentException("Empty data used to build a tree (" + this.getClass().getName() + ")");
            }
        }

        if (nodeCount >= 0) {
            initArrays();
        }

        // process trail time information, if any provided
        processTraits(m_traitList.get());

        // ensure tree is compatible with time trait
        if (timeTraitSet != null)
            adjustTreeNodeHeights(root);

        // ensure all nodes have their taxon names set up
        String[] taxa = getTaxaNames();
        for (int i = 0; i < getNodeCount() && i < taxa.length; i++) {
            if (taxa[i] != null) {
                if (m_nodes[i] == null) {
                    Log.warning("WARNING: Expected a node for taxon " + taxa[i] + " but did not find one in " +
                            "the expected location in the m_nodes array");
                } else if (m_nodes[i].getID() == null) {
                    m_nodes[i].setID(taxa[i]);
                }
            }
        }
    } // initAndValidate

    /**
     * Read a tree from a file.
     * The file should only contain one tree.
     * Support NEXUS and NEWICK.
     *
     * @param inputFileName apparently
     */
    public Tree readFromFile(String inputFileName) throws IOException, NotSingleException {
        List<Tree> parsedTrees;
        boolean isNexus = true;
        BufferedReader fin;
        String str;

        // Nexus tree or Newick tree?
        fin = new BufferedReader(new FileReader(inputFileName));
        if (!fin.ready()) {
            throw new IOException(inputFileName + " appears empty.");
        }
        str = fin.readLine();
        if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
            // Newick tree found
            isNexus = false;
        }

        // Parse tree
        if (isNexus) {
            NexusParser nexusParser = new ScsNexusParser(nodeTypeInput.get());
            nexusParser.parseFile(new File(inputFileName));
            parsedTrees = nexusParser.trees;
        } else {
            parsedTrees = new ArrayList<>();

            while (true) {
                Tree thisTree;

                try {
                    thisTree = new ScsTreeParser(null, str, 0, false, nodeTypeInput.get());
                } catch (ArrayIndexOutOfBoundsException e) {
                    thisTree = new ScsTreeParser(null, str, 1, false, nodeTypeInput.get());
                }

                parsedTrees.add(thisTree);

                if (fin.ready())
                    str = fin.readLine().trim();
                else
                    break;
            }
        }
        fin.close();

        // Verify; only one input tree is allowed
        if (parsedTrees.size() != 1) {
            throw new NotSingleException("Only one input tree is expected, but " + parsedTrees.size() + " provided.");
        }

        return parsedTrees.get(0);
    }

    /**
     * Initialize branch rates from the input tree.
     *
     * @param node the node currently working on
     */
    private void initialiseBranchRates(Node node) {
        Object rate;
        rate = prioritizeMedianRate ? node.getMetaData(MEDIAN_RATE_NAME) : node.getMetaData(MEAN_RATE_NAME);
        if (rate == null)
            rate = node.getMetaData(MEAN_RATE_NAME);
        if (rate == null)
            branchRates[node.getNr()] = 1.0;
        else
            branchRates[node.getNr()] = (Double) rate;

        if (!node.isLeaf()) {
            for (Node i : node.getChildren()) {
                initialiseBranchRates(i);
            }
        }
    }

    /**
     * Get the branch rate of node {@param node}.
     *
     * @param node which node to check?
     * @return branch rate
     */
    public double getBranchRate(Node node) {
        return branchRates[node.getNr()];
    }

    @Override
    @SuppressWarnings("deprecation")
    public void makeCaterpillar(final double minInternalHeight, final double step, final boolean finalize) {
        // make a caterpillar
        final List<String> taxa = m_taxonset.get().asStringList();

        // make a tumor subtree
        Node left = newNode();
        left.labelNr = 0;
        left.height = 0;
        left.setID(taxa.get(0));
        for (int i = 1; i < taxa.size(); i++) {
            final Node right = newNode();
            right.labelNr = i;
            right.height = 0;
            right.setID(taxa.get(i));
            final Node parent = newNode();
            parent.labelNr = taxa.size() + i - 1;
            parent.height = minInternalHeight + i * step;
            left.parent = parent;
            parent.setLeft(left);
            right.parent = parent;
            parent.setRight(right);
            left = parent;
        }
        rootTMRCA = left;

        // add the root and trunk
        root = newNode();
        root.labelNr = 2 * taxa.size() - 1;
        root.height = minInternalHeight + taxa.size() * step;
        rootTMRCA.parent = root;
        root.setLeft(rootTMRCA);

        // some statistics
        leafNodeCount = taxa.size();
        nodeCount = leafNodeCount * 2;
        internalNodeCount = leafNodeCount;

        if (finalize) {
            initArrays();
        }
    } // makeCaterpillar

    @Override
    public int scale(final double scale) {
        assert root instanceof ScsNode;

        return root.scale(scale);
    }

    /**
     * update rootTMRCA to be the child of root
     */
    @SuppressWarnings("deprecation")
    public void updateRootTMRCA() {
        if (rootTMRCA != root.getLeft()) {
            rootTMRCA = root.getLeft();
        }
    } // updateRootTMRCA

    /**
     * deep copy, returns a completely new tree
     *
     * @return a deep copy of this beast.tree.
     */
    @Override
    @SuppressWarnings("deprecation")
    public Tree copy() {
        ScsTree tree = new ScsTree();
        tree.setID(getID());
        tree.index = index;
        tree.root = root.copy();
        tree.rootTMRCA = tree.root.getLeft();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        return tree;
    }

    /**
     * copy of all values into existing tree *
     */
    @Override
    @SuppressWarnings("deprecation")
    public void assignTo(final StateNode other) {
        final ScsTree tree = (ScsTree) other;
        final Node[] nodes = new ScsNode[nodeCount];
        listNodes(tree.root, nodes);
        tree.setID(getID());
        //tree.index = index;
        root.assignTo(nodes);
        tree.root = nodes[root.getNr()];
        tree.rootTMRCA = tree.root.getLeft();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    @SuppressWarnings("deprecation")
    public void assignFrom(final StateNode other) {
        final ScsTree tree = (ScsTree) other;
        final Node[] nodes = new ScsNode[tree.getNodeCount()];//tree.getNodesAsArray();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = newNode();
        }
        setID(tree.getID());
        //index = tree.index;
        root = nodes[tree.root.getNr()];
        root.assignFrom(nodes, tree.root);
        root.parent = null;
        rootTMRCA = root.getLeft();
        nodeCount = tree.nodeCount;
        internalNodeCount = tree.internalNodeCount;
        leafNodeCount = tree.leafNodeCount;
        initArrays();
    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    @SuppressWarnings("deprecation")
    public void assignFromFragile(final StateNode other) {
        // invalidate cache
        postCache = null;

        final ScsTree tree = (ScsTree) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[tree.root.getNr()];
        final Node[] otherNodes = tree.m_nodes;
        final int rootNr = root.getNr();
        assignFrom(0, rootNr, otherNodes);
        root.height = otherNodes[rootNr].height;
        root.parent = null;
        if (otherNodes[rootNr].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[rootNr].getLeft().getNr()]);
            rootTMRCA = root.getLeft();
        } else {
            root.setLeft(null);
        }
        if (otherNodes[rootNr].getRight() != null) {
            root.setRight(m_nodes[otherNodes[rootNr].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        assignFrom(rootNr + 1, nodeCount, otherNodes);
    }

    /**
     * helper to assignFromFragile *
     */
    @SuppressWarnings("deprecation")
    private void assignFrom(final int start, final int end, final Node[] otherNodes) {
        for (int i = start; i < end; i++) {
            Node sink = m_nodes[i];
            Node src = otherNodes[i];
            sink.height = src.height;
            sink.parent = m_nodes[src.parent.getNr()];
            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }

    /**
     * StateNode implementation *
     */
    @Override
    protected void store() {

        // this condition can only be true for sampled ancestor trees
        if (m_storedNodes.length != nodeCount) {
            final Node[] tmp = new ScsNode[nodeCount];
            System.arraycopy(m_storedNodes, 0, tmp, 0, m_storedNodes.length - 1);
            if (nodeCount > m_storedNodes.length) {
                tmp[m_storedNodes.length - 1] = m_storedNodes[m_storedNodes.length - 1];
                tmp[nodeCount - 1] = newNode();
                tmp[nodeCount - 1].setNr(nodeCount - 1);
            }
            m_storedNodes = tmp;
        }

        storeNodes(0, nodeCount);
        storedRoot = m_storedNodes[root.getNr()];
    }

    /**
     * Stores nodes with index i, for start <= i < end
     * (i.e. including start but not including end)
     *
     * @param start the first index to be stored
     * @param end   nodes are stored up to but not including this index
     */
    private void storeNodes(final int start, final int end) {
        // Use direct members for speed (we are talking 5-7% or more from total time for large trees :)
        for (int i = start; i < end; i++) {
            final Node sink = m_storedNodes[i];
            final Node src = m_nodes[i];
            sink.height = src.height;

            if (src.parent != null) {
                sink.parent = m_storedNodes[src.parent.getNr()];
            } else {
                // currently only called in the case of sampled ancestor trees
                // where root node is not always last in the list
                sink.parent = null;
            }

            final List<Node> children = sink.children;
            final List<Node> srcChildren = src.children;

            if (children.size() == srcChildren.size()) {
                // shave some more time by avoiding list clear and add
                for (int k = 0; k < children.size(); ++k) {
                    final Node srcChild = srcChildren.get(k);
                    // don't call addChild, which calls  setParent(..., true);
                    final Node c = m_storedNodes[srcChild.getNr()];
                    c.parent = sink;
                    children.set(k, c);
                }
            } else {
                children.clear();
                //sink.removeAllChildren(false);
                for (final Node srcChild : srcChildren) {
                    // don't call addChild, which calls  setParent(..., true);
                    final Node c = m_storedNodes[srcChild.getNr()];
                    c.parent = sink;
                    children.add(c);
                    //sink.addChild(c);
                }
            }
        }
    }

    @Override
    @SuppressWarnings("deprecation")
    public void restore() {

        // necessary for sampled ancestor trees
        nodeCount = m_storedNodes.length;

        final Node[] tmp = m_storedNodes;
        m_storedNodes = m_nodes;
        m_nodes = tmp;
        root = m_nodes[storedRoot.getNr()];
        rootTMRCA = root.getLeft();

        // necessary for sampled ancestor trees,
        // we have the nodes, no need for expensive recursion
        leafNodeCount = 0;
        for (Node n : m_nodes) {
            leafNodeCount += n.isLeaf() ? 1 : 0;
        }

        //leafNodeCount = root.getLeafNodeCount();

        hasStartedEditing = false;

        for (Node n : m_nodes) {
            n.isDirty = Tree.IS_CLEAN;
        }

        postCache = null;
    }

}
