package beast.evolution.variantsinfo.vcfentry;

import beast.math.util.MathFunctions;
import org.jetbrains.annotations.NotNull;

import java.util.Arrays;

import static beast.math.util.MathFunctions.convertLogE2RoundedPhredScaled;

public class CellLocusGenotype {

    //**********************************************
    //*             Constant variables             *
    //**********************************************

    private final static String ITEM_SEPARATOR = ",";
    private final static String CLASS_SEPARATOR = ":";

    private final static String GENOTYPE_HEADER = "GT";
    private final static String ALT_ALLELE_DEPTH_HEADER = "AD";
    private final static String READ_DEPTH_HEADER = "DP";
    private final static String GENOTYPE_QUALITY_HEADER = "GQ";
    private final static String GENOTYPE_PHRED_SCALED_HEADER = "PL";
    private final static String ADO_STATE_HEADER = "ADO";
    private final static String ADO_QUALITY_HEADER = "ADOQ";
    private final static String ADO_PHRED_SCALED_HEADER = "ADOP";

    private final static String COMBINED_HEADER;

    static {
        COMBINED_HEADER = getCombinedHeader();
    }


    //*********************************************
    //*                 Variables                 *
    //*********************************************

    private final String genotypeOri;
    private final String genotypeAdapted;
    private final int ternary;
    private int[] alleleDepth;
    private final int readDepth;
    private final double[] genotypeLogLikelihoods;
    private final long genotypeQual;
    private final long[] genotypePhredScaledLikelihoods;
    private final int adoState;
    private final double[] adoLogLikelihoods;
    private final long adoQual;
    private final long[] adoPhredScaledLikelihoods;


    //********************************************
    //*               Constructors               *
    //********************************************

    public CellLocusGenotype(
            @NotNull String genotypeOri,
            @NotNull String genotypeAdapted,
            int ternary,
            int[] alleleDepth,
            int readDepth,
            double[] genotypeLogLikelihoods,
            int adoState,
            double[] adoLogLikelihoods
    ) {
        this.genotypeOri = genotypeOri;
        this.genotypeAdapted = genotypeAdapted;
        this.ternary = ternary;
        this.alleleDepth = new int[alleleDepth.length];
        this.readDepth = readDepth;
        this.genotypeLogLikelihoods = new double[genotypeLogLikelihoods.length];
        this.genotypePhredScaledLikelihoods = new long[genotypeLogLikelihoods.length];
        this.adoState = adoState;
        this.adoLogLikelihoods = new double[adoLogLikelihoods.length];
        this.adoPhredScaledLikelihoods = new long[adoLogLikelihoods.length];
        long[] genotypePLTmp = new long[genotypeLogLikelihoods.length];
        long[] adoPLTmp = new long[adoLogLikelihoods.length];

        System.arraycopy(alleleDepth, 0, this.alleleDepth, 0, alleleDepth.length);
        System.arraycopy(genotypeLogLikelihoods, 0, this.genotypeLogLikelihoods, 0, genotypeLogLikelihoods.length);
        System.arraycopy(adoLogLikelihoods, 0, this.adoLogLikelihoods, 0, adoLogLikelihoods.length);

        for (int i = 0; i < genotypeLogLikelihoods.length; i++)
            this.genotypePhredScaledLikelihoods[i] = convertLogE2RoundedPhredScaled(genotypeLogLikelihoods[i]);

        long genotypePLMin = MathFunctions.min(
                Arrays.stream(this.genotypePhredScaledLikelihoods).boxed().toArray(Long[]::new)
        );

        for (int i = 0; i < genotypeLogLikelihoods.length; i++)
            this.genotypePhredScaledLikelihoods[i] -= genotypePLMin;

        System.arraycopy(this.genotypePhredScaledLikelihoods, 0, genotypePLTmp, 0, genotypeLogLikelihoods.length);
        Arrays.sort(genotypePLTmp);
        this.genotypeQual = genotypePLTmp[1] > 99 ? 99 : genotypePLTmp[1];

        for (int i = 0; i < adoLogLikelihoods.length; i++)
            this.adoPhredScaledLikelihoods[i] = convertLogE2RoundedPhredScaled(adoLogLikelihoods[i]);

        long adoPLMin = MathFunctions.min(
                Arrays.stream(this.adoPhredScaledLikelihoods).boxed().toArray(Long[]::new)
        );

        for (int i = 0; i < adoLogLikelihoods.length; i++)
            this.adoPhredScaledLikelihoods[i] -= adoPLMin;

        System.arraycopy(this.adoPhredScaledLikelihoods, 0, adoPLTmp, 0, adoLogLikelihoods.length);
        Arrays.sort(adoPLTmp);
        this.adoQual = adoPLTmp[1] > 99 ? 99 : adoPLTmp[1];
    }


    //*********************************************
    //*                  Methods                  *
    //*********************************************

    private static String getCombinedHeader() {
        StringBuilder str = new StringBuilder();

        str.append(GENOTYPE_HEADER).append(CLASS_SEPARATOR);
        str.append(ALT_ALLELE_DEPTH_HEADER).append(CLASS_SEPARATOR);
        str.append(READ_DEPTH_HEADER).append(CLASS_SEPARATOR);
        str.append(GENOTYPE_QUALITY_HEADER).append(CLASS_SEPARATOR);
        str.append(GENOTYPE_PHRED_SCALED_HEADER).append(CLASS_SEPARATOR);
        str.append(ADO_STATE_HEADER).append(CLASS_SEPARATOR);
        str.append(ADO_QUALITY_HEADER).append(CLASS_SEPARATOR);
        str.append(ADO_PHRED_SCALED_HEADER);

        return str.toString();
    } // getCombinedHeader

    public static String getHeader() {
        return COMBINED_HEADER;
    } // getHeader

    public void setAlleleDepth(int[] alleleDepth) {
        this.alleleDepth = new int[alleleDepth.length];
        System.arraycopy(alleleDepth, 0, this.alleleDepth, 0, alleleDepth.length);
    } // alleleDepth

    public int getAlleleDepthLength() {
        return this.alleleDepth.length;
    } // getAlleleDepthLength

    public int getAdoState() {
        return this.adoState;
    } // getAdoState

    public String getGenotypeOri() {
        return this.genotypeOri;
    } // getGenotypeOri

    public int getTernary() {
        return this.ternary;
    } // getTernary

    public String getGenotypeLikelihoodsAsString(String separator, String floatFormat) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.genotypeLogLikelihoods.length; i++) {
            sb.append(String.format(floatFormat, Math.exp(this.genotypeLogLikelihoods[i])));

            if (i < this.genotypeLogLikelihoods.length - 1)
                sb.append(separator);
        }

        return sb.toString();
    } // getGenotypeLikelihoodsAsString

    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();

        str.append(this.genotypeAdapted).append(CLASS_SEPARATOR);
        str.append(
                String.join(
                        ITEM_SEPARATOR,
                        Arrays.stream(this.alleleDepth).mapToObj(String::valueOf).toArray(String[]::new)
                )
        ).append(CLASS_SEPARATOR);
        str.append(this.readDepth).append(CLASS_SEPARATOR);
        str.append(this.genotypeQual).append(CLASS_SEPARATOR);
        str.append(
                String.join(
                        ITEM_SEPARATOR,
                        Arrays.stream(this.genotypePhredScaledLikelihoods).mapToObj(String::valueOf).toArray(String[]::new)
                )
        ).append(CLASS_SEPARATOR);
        str.append(this.adoState).append(CLASS_SEPARATOR);
        str.append(this.adoQual).append(CLASS_SEPARATOR);
        str.append(
                String.join(
                        ITEM_SEPARATOR,
                        Arrays.stream(this.adoPhredScaledLikelihoods).mapToObj(String::valueOf).toArray(String[]::new)
                )
        );

        return str.toString();
    } // toString

}
