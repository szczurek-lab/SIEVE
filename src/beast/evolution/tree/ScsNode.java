package beast.evolution.tree;

import beast.core.Description;

import java.util.TreeMap;

@Description("")
public class ScsNode extends Node {

    /**
     * figure out whether current node is the most recent common ancestor of all tumor cells (TMRCA) or not
     * if yes, it is the only child of its parent (MRCA, tree root)
     *
     * @return true if current node is TMRCA
     */
    public boolean isTMRCA() {
        return parent.isRoot() && parent.children.size() == 1;
    } // isTMRCA

    /**
     * @return (deep) copy of node
     */
    @Override
    public Node copy() {
        final Node node = new ScsNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.lengthMetaDataString = lengthMetaDataString;
        node.metaData = new TreeMap<>(metaData);
        node.lengthMetaData = new TreeMap<>(lengthMetaData);
        node.parent = null;
        node.setID(getID());

        for (final Node child : getChildren()) {
            node.addChild(child.copy());
        }
        return node;
    } // copy

    /**
     * scale height of this node and all its descendants
     *
     * @param scale scale factor
     * @return degrees of freedom scaled (used for HR calculations)
     */
    @Override
    @SuppressWarnings("deprecation")
    public int scale(final double scale) {
        startEditing();

        int dof = 0;

        isDirty |= Tree.IS_DIRTY;
        if (!isLeaf() && !isFake()) {
            height *= scale;

            if (isRoot() || parent.getHeight() != getHeight())
                dof += 1;
        }
        if (!isLeaf()) {
            dof += getLeft().scale(scale);
            if (getRight() != null) {
                dof += getRight().scale(scale);
            }
            if (height < getLeft().height || (getRight() != null && height < getRight().height)) {
                throw new IllegalArgumentException("Scale gives negative branch length");
            }
        }

        return dof;
    }

}
