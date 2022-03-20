package beast.evolution.variantsinfo.vcfentry;

import beast.evolution.alignment.VariantSiteInfo;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;

public class VCFEntry {


    //**********************************************
    //*             Constant variables             *
    //**********************************************

    private final static String ERROR_OCCUPIED = "Error. Reassigning a non-null element.";

    private final static String CLASS_SEPARATOR = "\t";

    final static String FLOAT_FORMAT = "%.3f";

    private final static String CHROMOSOME_HEADER = "CHROM";
    private final static String POS_HEADER = "POS";
    private final static String ID_HEADER = "ID";
    private final static String REF_HEADER = "REF";
    private final static String ALT_HEADER = "ALT";
    private final static String QUAL_HEADER = "QUAL";
    private final static String FILTER_HEADER = "FILTER";
    private final static String INFO_HEADER = "INFO";
    private final static String FORMAT_HEADER = "FORMAT";
    private static final List<String> cellNames = new ArrayList<>();


    //*********************************************
    //*                 Variables                 *
    //*********************************************

    private final VariantSiteInfo varInfo;
    private final String id;
    private final long qual;
    private final String filter;
    private final VCFInfo info;
    private final CellLocusGenotype[] cellLocusGenotype;


    //********************************************
    //*               Constructors               *
    //********************************************

    public VCFEntry(
            @NotNull VariantSiteInfo varInfo,
            String id,
            long qual,
            String filter,
            VCFInfo info,
            final CellLocusGenotype[] cellLocusGenotype
    ) {
        this.varInfo = varInfo;
        this.id = id;
        this.qual = qual;
        this.filter = filter;
        this.info = info;
        this.cellLocusGenotype = cellLocusGenotype;
    }


    //*********************************************
    //*                  Methods                  *
    //*********************************************

    public static void addCellNames(List<String> names, boolean overwrite) {
        if (cellNames.size() != 0) {
            if (overwrite)
                cellNames.clear();
            else
                return;
        }

        cellNames.addAll(names);
    } // addCellNames

    public static List<String> getCellNames() {
        return cellNames;
    } // getCellNames

    public String getAdoStateAsString(String separator) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.cellLocusGenotype.length; i++) {
            sb.append(this.cellLocusGenotype[i].getAdoState());

            if (i < this.cellLocusGenotype.length - 1)
                sb.append(separator);
        }

        return sb.toString();
    } // getAdoStateAsString

    public String getCellLocusOriGenotypeAsString(String separator) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.cellLocusGenotype.length; i++) {
            sb.append(this.cellLocusGenotype[i].getGenotypeOri());

            if (i < this.cellLocusGenotype.length - 1)
                sb.append(separator);
        }

        return sb.toString();
    } // getCellLocusOriGenotypeAsString

    public String getCellLocusTernaryGenotypeAsString(String separator) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.cellLocusGenotype.length; i++) {
            sb.append(this.cellLocusGenotype[i].getTernary());

            if (i < this.cellLocusGenotype.length - 1)
                sb.append(separator);
        }

        return sb.toString();
    } // getCellLocusTernaryGenotypeAsString

    public String getCellLocusGenotypeLikelihoodAsString(String classSeparator, String itemSeparator) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.cellLocusGenotype.length; i++) {
            sb.append(this.cellLocusGenotype[i].getGenotypeLikelihoodsAsString(itemSeparator, FLOAT_FORMAT));

            if (i < this.cellLocusGenotype.length - 1)
                sb.append(classSeparator);
        }

        return sb.toString();
    } // getCellLocusGenotypeLikelihoodAsString

    public String getLociInfoAsString(String classSeparator, String itemSeparator) {
        StringBuilder sb = new StringBuilder();

        sb.append(this.varInfo.getChromosome()).append(classSeparator);
        sb.append(this.varInfo.getPosition()).append(classSeparator);
        sb.append(this.varInfo.getRefNuc()).append(classSeparator);
        sb.append(this.varInfo.getAltNucsAsString(itemSeparator));

        return sb.toString();
    } // getLociInfoAsString

    public static String getCombinedHeader() {
        if (cellNames.size() == 0)
            throw new IllegalArgumentException("Error! No cell names set.");

        StringBuilder sb = new StringBuilder();

        sb.append("#");
        sb.append(CHROMOSOME_HEADER).append(CLASS_SEPARATOR);
        sb.append(POS_HEADER).append(CLASS_SEPARATOR);
        sb.append(ID_HEADER).append(CLASS_SEPARATOR);
        sb.append(REF_HEADER).append(CLASS_SEPARATOR);
        sb.append(ALT_HEADER).append(CLASS_SEPARATOR);
        sb.append(QUAL_HEADER).append(CLASS_SEPARATOR);
        sb.append(FILTER_HEADER).append(CLASS_SEPARATOR);
        sb.append(INFO_HEADER).append(CLASS_SEPARATOR);
        sb.append(FORMAT_HEADER);

        for (String name : cellNames) {
            sb.append(CLASS_SEPARATOR).append(name);
        }

        return sb.toString();
    } // getCombinedHeader

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        // chrom
        sb.append(this.varInfo.getChromosome()).append(CLASS_SEPARATOR);

        // pos
        sb.append(this.varInfo.getPosition()).append(CLASS_SEPARATOR);

        // ID
        sb.append(this.id).append(CLASS_SEPARATOR);

        // ref nuc
        sb.append(this.varInfo.getRefNuc()).append(CLASS_SEPARATOR);

        // alt nucs
        sb.append(this.varInfo.getAltNucsAsString(",")).append(CLASS_SEPARATOR);

        // qual
        sb.append(this.qual).append(CLASS_SEPARATOR);

        // filter
        sb.append(this.filter).append(CLASS_SEPARATOR);

        // info
        sb.append(this.info.toString()).append(CLASS_SEPARATOR);

        // format
        sb.append(CellLocusGenotype.getHeader());

        // cell locus genotype
        for (CellLocusGenotype i : this.cellLocusGenotype)
            sb.append(CLASS_SEPARATOR).append(i.toString());

        return sb.toString();
    } // toString

}
