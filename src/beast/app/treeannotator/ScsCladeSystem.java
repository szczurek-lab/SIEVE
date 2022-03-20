package beast.app.treeannotator;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;

public class ScsCladeSystem extends CladeSystem {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private boolean processSA = true;

    private boolean hasTrunk = true;

    /**
     * only not null when the tree has a trunk
     */
    private List<Object[]> rootAttributeValues = null;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public ScsCladeSystem(boolean hasTrunk) {
        this.hasTrunk = hasTrunk;
        if (hasTrunk) {
            rootAttributeValues = new ArrayList<>();
        }
    }

    public ScsCladeSystem(boolean hasTrunk, Tree targetTree) {
        this.hasTrunk = hasTrunk;
        if (hasTrunk) {
            rootAttributeValues = new ArrayList<>();
        }

        add(targetTree, true);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * adds all the clades in the tree
     */
    @Override
    public void add(Tree tree, boolean includeTips) {
        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).
        addClades(tree.getRoot(), includeTips);
    }

    private BitSet addClades(Node node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2 * index);

            if (includeTips) {
                addClade(bits);
            }

        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                bits.or(addClades(node1, includeTips));
            }

            for (int i = 1; i < bits.length(); i = i + 2) {
                bits.set(i, false);
            }
            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits.set(2 * index + 1);
            }

            if (!node.isRoot() || (node.isRoot() && !hasTrunk)) {
                addClade(bits);
            }
        }

        return bits;
    }

    private void addClade(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade == null) {
            clade = new Clade(bits);
            cladeMap.put(bits, clade);
        }
        clade.setCount(clade.getCount() + 1);
    }

    @Override
    public double getSumCladeCredibility(Node node, BitSet bits) {

        double sum = 0.0;

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2 * index);
        } else {

            BitSet bits2 = new BitSet();
            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                sum += getSumCladeCredibility(node1, bits2);
            }

            for (int i = 1; i < bits2.length(); i = i + 2) {
                bits2.set(i, false);
            }

            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits2.set(2 * index + 1);
            }

            if (!node.isRoot() || (node.isRoot() && !hasTrunk)) {
                sum += getCladeCredibility(bits2);
            }

            if (bits != null) {
                bits.or(bits2);
            }
        }

        return sum;
    }

    @Override
    public double getLogCladeCredibility(Node node, BitSet bits) {

        double logCladeCredibility = 0.0;

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2 * index);
        } else {

            BitSet bits2 = new BitSet();
            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                logCladeCredibility += getLogCladeCredibility(node1, bits2);
            }

            for (int i = 1; i < bits2.length(); i = i + 2) {
                bits2.set(i, false);
            }

            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits2.set(2 * index + 1);
            }

            if (!node.isRoot() || (node.isRoot() && !hasTrunk)) {
                logCladeCredibility += Math.log(getCladeCredibility(bits2));
            }

            if (bits != null) {
                bits.or(bits2);
            }
        }

        return logCladeCredibility;
    }

    private double getCladeCredibility(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade == null) {
            return 0.0;
        }
        return clade.getCredibility();
    }

    @Override
    public void collectAttributes(Tree tree, Set<String> attributeNames) {
        collectAttributes(tree.getRoot(), attributeNames);
    }

    private BitSet collectAttributes(Node node, Set<String> attributeNames) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            if (index < 0) {
                throw new IllegalArgumentException("Taxon, " + node.getID() + ", not found in target tree");
            }
            bits.set(2 * index);

        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                bits.or(collectAttributes(node1, attributeNames));
            }

            for (int i = 1; i < bits.length(); i = i + 2) {
                bits.set(i, false);
            }
            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits.set(2 * index + 1);
            }
        }

        if (node.isRoot() && hasTrunk) {
            collectAttributesForRoot(node, attributeNames);
        } else {
            collectAttributesForClade(bits, node, attributeNames);
        }

        return bits;
    }

    private void collectAttributesForClade(BitSet bits, Node node, Set<String> attributeNames) {
        Clade clade = cladeMap.get(bits);
        if (clade != null) {

            if (clade.attributeValues == null) {
                clade.attributeValues = new ArrayList<>();
            }

            int i = 0;
            Object[] values = new Object[attributeNames.size()];
            for (String attributeName : attributeNames) {

                Object value;
                switch (attributeName) {
                    case "height":
                        value = node.getHeight();
                        break;
                    case "length":
                        value = getBranchLength(node);
                        break;
                    default:
                        value = node.getMetaData(attributeName);
                        if (value instanceof String && ((String) value).startsWith("\"")) {
                            value = ((String) value).replaceAll("\"", "");
                        }
                        break;
                }

                values[i] = value;

                i++;
            }
            clade.attributeValues.add(values);

            clade.setCount(clade.getCount() + 1);
        }
    }

    /**
     * collect attributes for root if the tree has a trunk
     *
     * @param node           tree root
     * @param attributeNames apparently
     */
    private void collectAttributesForRoot(Node node, Set<String> attributeNames) {
        assert rootAttributeValues != null : "clade system is not properly initialized";

        int i = 0;
        Object[] values = new Object[attributeNames.size()];
        for (String attributeName : attributeNames) {
            Object value;

            switch (attributeName) {
                case "height":
                    value = node.getHeight();
                    break;
                case "length":
                    value = 0.0;
                    break;
                default:
                    value = node.getMetaData(attributeName);
                    if (value instanceof String && ((String) value).startsWith("\"")) {
                        value = ((String) value).replaceAll("\"", "");
                    }
                    break;
            }

            values[i] = value;
            i++;
        }

        rootAttributeValues.add(values);
    }

    private Object getBranchLength(Node node) {
        if (node.isRoot()) {
            return 0;
        }
        return node.getParent().getHeight() - node.getHeight();
    }

    @Override
    public BitSet removeClades(Node node, boolean includeTips) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = getTaxonIndex(node);
            bits.set(2 * index);

            if (includeTips) {
                removeClade(bits);
            }

        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                bits.or(removeClades(node1, includeTips));
            }

            for (int i = 1; i < bits.length(); i = i + 2) {
                bits.set(i, false);
            }
            if (node.isFake() && processSA) {
                int index = getTaxonIndex(node.getDirectAncestorChild());
                bits.set(2 * index + 1);
            }

            if (!node.isRoot() || (node.isRoot() && !hasTrunk)) {
                removeClade(bits);
            }
        }

        return bits;
    }

    private void removeClade(BitSet bits) {
        Clade clade = cladeMap.get(bits);
        if (clade != null) {
            clade.setCount(clade.getCount() - 1);
        }

    }

    // Get tree clades as bitSets on target taxa
    // codes is an array of existing BitSet objects, which are reused
    @Override
    void getTreeCladeCodes(Tree tree, BitSet[] codes) {
        if (hasTrunk) {
            assert tree.getRoot().getChildCount() == 1;
            getTreeCladeCodes(tree.getRoot().getChild(0), codes);
        } else {
            getTreeCladeCodes(tree.getRoot(), codes);
        }
    }

    public int getTaxonIndex(Node node) {
        return node.getNr();
    }

    public void setProcessSA(boolean processSA) {
        this.processSA = processSA;
    }

    public List<Object[]> getRootAttributeValues() {
        return rootAttributeValues;
    }

}
