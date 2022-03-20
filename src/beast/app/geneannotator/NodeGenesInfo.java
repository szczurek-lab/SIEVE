package beast.app.geneannotator;

import beast.evolution.alignment.VariantSiteInfo;

import java.util.ArrayList;
import java.util.List;

import static beast.evolution.substitutionmodel.ScsSubstitutionModelBase.EvolutionaryEventType.arrayToString;

public class NodeGenesInfo {

    private List<GeneInfo> fullGenes;
    private List<GeneInfo> ISAGenes;
    private List<GeneInfo> FSAGenes;

    private List<String> fullChr;
    private List<Long> fullPos;
    private List<String> fullRefNuc;
    private List<String> fullAltNuc;
    private List<String> fullEventType;
    private List<String> fullGeneName;

    private List<String> ISAChr;
    private List<Long> ISAPos;
    private List<String> ISARefNuc;
    private List<String> ISAAltNuc;
    private List<String> ISAEventType;
    private List<String> ISAGeneName;

    private List<String> FSAChr;
    private List<Long> FSAPos;
    private List<String> FSARefNuc;
    private List<String> FSAAltNuc;
    private List<String> FSAEventType;
    private List<String> FSAGeneName;

    public NodeGenesInfo() {
        this.fullGenes = new ArrayList<>();
        this.ISAGenes = new ArrayList<>();
        this.FSAGenes = new ArrayList<>();
    }

    public void addGeneInfo(GeneInfo g) {
        if (!this.fullGenes.contains(g))
            this.fullGenes.add(g);
    }

    public void addISAGeneInfo(GeneInfo g) {
        if (!this.ISAGenes.contains(g))
            this.ISAGenes.add(g);
    }

    public void addFSAGeneInfo(GeneInfo g) {
        if (!this.FSAGenes.contains(g))
            this.FSAGenes.add(g);
    }

    public void convertGeneInfo() {
        if (this.fullGenes.size() == 0)
            return;

        final GeneInfo.GeneNamesAlphaComparator comp = new GeneInfo.GeneNamesAlphaComparator();

        this.fullGenes.sort(comp);

        this.fullChr = new ArrayList<>();
        this.fullPos = new ArrayList<>();
        this.fullRefNuc = new ArrayList<>();
        this.fullAltNuc = new ArrayList<>();
        this.fullEventType = new ArrayList<>();
        this.fullGeneName = new ArrayList<>();
        for (GeneInfo i : this.fullGenes) {
            final VariantSiteInfo j = i.getSnv();
            this.fullChr.add(j.getChromosome());
            this.fullPos.add(j.getPosition());
            this.fullRefNuc.add(String.valueOf(j.getRefNuc()));
            this.fullAltNuc.add(j.getAltNucsAsString(""));

            this.fullEventType.add(arrayToString(i.getEvoEventType()));
            this.fullGeneName.add(i.getGene());
        }

        if (this.ISAGenes.size() != 0) {
            this.ISAGenes.sort(comp);

            this.ISAChr = new ArrayList<>();
            this.ISAPos = new ArrayList<>();
            this.ISARefNuc = new ArrayList<>();
            this.ISAAltNuc = new ArrayList<>();
            this.ISAEventType = new ArrayList<>();
            this.ISAGeneName = new ArrayList<>();
            for (GeneInfo i : this.ISAGenes) {
                final VariantSiteInfo j = i.getSnv();
                this.ISAChr.add(j.getChromosome());
                this.ISAPos.add(j.getPosition());
                this.ISARefNuc.add(String.valueOf(j.getRefNuc()));
                this.ISAAltNuc.add(j.getAltNucsAsString(""));

                this.ISAEventType.add(arrayToString(i.getEvoEventType()));
                this.ISAGeneName.add(i.getGene());
            }
        }

        if (this.FSAGenes.size() != 0) {
            this.FSAGenes.sort(comp);

            this.FSAChr = new ArrayList<>();
            this.FSAPos = new ArrayList<>();
            this.FSARefNuc = new ArrayList<>();
            this.FSAAltNuc = new ArrayList<>();
            this.FSAEventType = new ArrayList<>();
            this.FSAGeneName = new ArrayList<>();
            for (GeneInfo i : this.FSAGenes) {
                final VariantSiteInfo j = i.getSnv();
                this.FSAChr.add(j.getChromosome());
                this.FSAPos.add(j.getPosition());
                this.FSARefNuc.add(String.valueOf(j.getRefNuc()));
                this.FSAAltNuc.add(j.getAltNucsAsString(""));

                this.FSAEventType.add(arrayToString(i.getEvoEventType()));
                this.FSAGeneName.add(i.getGene());
            }
        }
    }

    public List<GeneInfo> getFullGenes() {
        return fullGenes;
    }

    public List<GeneInfo> getISAGenes() {
        return ISAGenes;
    }

    public List<GeneInfo> getFSAGenes() {
        return FSAGenes;
    }

    public List<String> getFullChr() {
        return fullChr;
    }

    public List<Long> getFullPos() {
        return fullPos;
    }

    public List<String> getFullRefNuc() {
        return fullRefNuc;
    }

    public List<String> getFullAltNuc() {
        return fullAltNuc;
    }

    public List<String> getFullEventType() {
        return fullEventType;
    }

    public List<String> getFullGeneName() {
        return fullGeneName;
    }

    public List<String> getISAChr() {
        return ISAChr;
    }

    public List<Long> getISAPos() {
        return ISAPos;
    }

    public List<String> getISARefNuc() {
        return ISARefNuc;
    }

    public List<String> getISAAltNuc() {
        return ISAAltNuc;
    }

    public List<String> getISAEventType() {
        return ISAEventType;
    }

    public List<String> getISAGeneName() {
        return ISAGeneName;
    }

    public List<String> getFSAChr() {
        return FSAChr;
    }

    public List<Long> getFSAPos() {
        return FSAPos;
    }

    public List<String> getFSARefNuc() {
        return FSARefNuc;
    }

    public List<String> getFSAAltNuc() {
        return FSAAltNuc;
    }

    public List<String> getFSAEventType() {
        return FSAEventType;
    }

    public List<String> getFSAGeneName() {
        return FSAGeneName;
    }

}
