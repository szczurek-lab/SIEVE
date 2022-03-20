package beast.evolution.variantsinfo.vcfentry;

import java.util.ArrayList;
import java.util.List;

import static beast.evolution.variantsinfo.vcfentry.VCFEntry.FLOAT_FORMAT;

public class VCFInfo {


    //**********************************************
    //*             Constant variables             *
    //**********************************************

    private final static String ITEM_SEPARATOR = ",";
    private final static String CLASS_SEPARATOR = ";";

    private final static String ALLELE_COUNT_HEADER = "AC";
    private final static String ALLELE_FREQUENCY_HEADER = "AF";
    private final static String TOTAL_ALLELE_NUMS_HEADER = "AN";
    private final static String TOTAL_READ_DEPTH_HEADER = "DP";


    //*********************************************
    //*                 Variables                 *
    //*********************************************

    private final int[] alleleCount;
    private final double[] alleleFrequency;
    private final int totalAlleleNums;
    private final int totalReadDepth;


    //********************************************
    //*               Constructors               *
    //********************************************

    public VCFInfo(
            int[] alleleCount,
            int totalAlleleNums,
            int totalReadDepth
    ) {
        this.alleleCount = new int[alleleCount.length];
        this.alleleFrequency = new double[alleleCount.length];

        System.arraycopy(alleleCount, 0, this.alleleCount, 0, alleleCount.length);

        for (int i = 0; i < alleleCount.length; i++) {
            this.alleleFrequency[i] = ((double) this.alleleCount[i]) / totalAlleleNums;
        }

        this.totalAlleleNums = totalAlleleNums;
        this.totalReadDepth = totalReadDepth;
    }


    //*********************************************
    //*                  Methods                  *
    //*********************************************

    @Override
    public String toString() {
        StringBuilder str = new StringBuilder();
        List<String> tmp = new ArrayList<>();

        // AC
        str.append(ALLELE_COUNT_HEADER).append("=");
        for (int ac : this.alleleCount) {
            tmp.add(String.valueOf(ac));
        }
        str.append(String.join(ITEM_SEPARATOR, tmp.toArray(new String[0]))).append(CLASS_SEPARATOR);

        tmp.clear();

        // AF
        str.append(ALLELE_FREQUENCY_HEADER).append("=");
        for (double af : this.alleleFrequency) {
            tmp.add(String.format(FLOAT_FORMAT, af));
        }
        str.append(String.join(ITEM_SEPARATOR, tmp.toArray(new String[0]))).append(CLASS_SEPARATOR);

        // AN
        str.append(TOTAL_ALLELE_NUMS_HEADER).append("=").append(this.totalAlleleNums).append(CLASS_SEPARATOR);

        // DP
        str.append(TOTAL_READ_DEPTH_HEADER).append("=").append(this.totalReadDepth);

        return str.toString();
    } // toString

}
