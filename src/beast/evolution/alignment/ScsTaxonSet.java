package beast.evolution.alignment;

import beast.core.Description;
import beast.core.Input;

import java.util.*;

@Description("A TaxonSet is an ordered set of taxa. The order on the taxa is provided at the time of construction " +
        "either from a list of taxon objects or an alignment.")
public class ScsTaxonSet extends TaxonSet {
    final public Input<ScsAlignment> scsAlignmentInput = new Input<>("scsAlignment", "single cell " +
            "sequencing data alignment where each sequence represents a taxon");

    public ScsTaxonSet() {
    }

    public ScsTaxonSet(final List<Taxon> taxa) {
        taxonsetInput.setValue(taxa, this);
        initAndValidate();
    }

    public ScsTaxonSet(final ScsAlignment alignment) {
        scsAlignmentInput.setValue(alignment, this);
        initAndValidate();
    }


    public ScsTaxonSet(final String id, final List<Taxon> taxa) {
        setID(id);
        taxonsetInput.setValue(taxa, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        if (alignmentInput.get() != null) {
            throw new IllegalArgumentException("The alignment attribute should not be used for single cell " +
                    "sequencing data (in class " + this.getClass().getName() + ").");
        }

        taxonList = taxonsetInput.get();
        if (scsAlignmentInput.get() != null) {
            if (taxonList.size() > 0) {
                throw new IllegalArgumentException("Only one of taxon and scsAlignmentInput should be specified, not both (id=" + getID() + ").");
            }
            taxaNames = scsAlignmentInput.get().taxaNames;
        } else {
            if (taxonList.size() == 0) {
                throw new IllegalArgumentException(getID() + ": Either taxon or scsAlignmentInput should be specified (id=" + getID() + ").");
            }
            taxaNames = new ArrayList<>();
            for (final Taxon taxon : taxonList) {
                taxaNames.add(taxon.getID());
            }
        }
    } // initAndValidate

}
