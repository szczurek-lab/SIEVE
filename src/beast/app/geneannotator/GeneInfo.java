package beast.app.geneannotator;

import beast.evolution.alignment.VariantSiteInfo;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase.EvolutionaryEventType;
import org.jetbrains.annotations.NotNull;

import java.util.Arrays;
import java.util.Comparator;

public class GeneInfo implements Comparable<GeneInfo> {

    private final VariantSiteInfo snv;
    private final EvolutionaryEventType[] evoEventType;
    private final String gene;

    public GeneInfo(
            VariantSiteInfo snv,
            EvolutionaryEventType[] evoEventType,
            String gene
    ) {
        this.snv = snv;
        this.evoEventType = new EvolutionaryEventType[evoEventType.length];
        System.arraycopy(evoEventType, 0, this.evoEventType, 0, evoEventType.length);
        this.gene = gene;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        GeneInfo geneInfo = (GeneInfo) o;

        if (!snv.equals(geneInfo.snv)) return false;
//        if (!Arrays.equals(evoEventType, geneInfo.evoEventType)) return false;
        return gene.equals(geneInfo.gene);
    }

    @Override
    public int hashCode() {
        int result = snv.hashCode();
        for (EvolutionaryEventType i : evoEventType) {
            result = 31 * result + i.toString().hashCode();
        }
        result = 31 * result + gene.hashCode();
        return result;
    }

    public VariantSiteInfo getSnv() {
        return this.snv;
    }

    public EvolutionaryEventType[] getEvoEventType() {
        return this.evoEventType;
    }

    public String getGene() {
        return this.gene;
    }

    public boolean violateISA() {
        if (this.evoEventType.length > 1) return true;
        return this.evoEventType[0].violateISA();
    }

    @Override
    public int compareTo(@NotNull GeneInfo g) {
        final int v1 = this.snv.compareTo(g.snv);
        if (v1 != 0) return v1;

        if (this.evoEventType.length > g.evoEventType.length)
            return 1;
        else if (this.evoEventType.length < g.evoEventType.length)
            return -1;
        else {
            for (int i = 0; i < this.evoEventType.length; i++) {
                final int v2 = this.evoEventType[i].compare(g.evoEventType[i]);
                if (v2 != 0) return v2;
            }
        }

        return this.gene.compareTo(g.gene);
    }

    public static class GeneInfoComparator implements Comparator<GeneInfo> {

        @Override
        public int compare(GeneInfo g1, GeneInfo g2) {
            return g1.compareTo(g2);
        }

    }

    public static class GeneNamesAlphaComparator implements Comparator<GeneInfo> {

        @Override
        public int compare(GeneInfo g1, GeneInfo g2) {
            return g1.gene.compareTo(g2.gene);
        }

    }

}
