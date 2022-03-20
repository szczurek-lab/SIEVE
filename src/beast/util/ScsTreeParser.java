package beast.util;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;

import java.util.List;

public class ScsTreeParser extends TreeParser {

    /**
     * Create a tree from the given newick format
     *
     * @param taxaNames                             a list of taxa names to use, or null.
     *                                              If null then IsLabelledNewick will be set to true
     * @param newick                                the newick of the tree
     * @param offset                                the offset to map node numbers in newick format to indices in taxaNames.
     *                                              so, name(node with nodeNumber) = taxaNames[nodeNumber-offset]
     * @param adjustTipHeightsWhenMissingDateTraits true if tip heights should be adjusted to zero
     */
    public ScsTreeParser(
            final List<String> taxaNames,
            final String newick,
            final int offset,
            final boolean adjustTipHeightsWhenMissingDateTraits,
            final String nodeType
    ) {
        if (taxaNames == null) {
            isLabelledNewickInput.setValue(true, this);
        } else {
            m_taxonset.setValue(new TaxonSet(Taxon.createTaxonList(taxaNames)), this);
        }
        newickInput.setValue(newick, this);
        offsetInput.setValue(offset, this);
        adjustTipHeightsInput.setValue(adjustTipHeightsWhenMissingDateTraits, this);
        nodeTypeInput.setValue(nodeType, this);
        labels = taxaNames;
        initAndValidate();
    }

}
