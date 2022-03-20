package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;

import java.lang.reflect.InvocationTargetException;

@Description(value = "A finite mutation, deletion and insertion substitution model with 15 genotypes")
public class ScsFiniteMuDelInsModel extends ScsSubstitutionModelBase {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public final Input<RealParameter> deletionRateInput = new Input<>("deletionRate", "deletion rate " +
            "(non-negative; default 0); measured relatively to mutation rate", Input.Validate.REQUIRED);
    public final Input<RealParameter> insertionRateInput = new Input<>("insertionRate", "insertion rate " +
            "(non-negative; default 0); measured relatively to mutation rate", Input.Validate.REQUIRED);

    protected RealParameter deletionRate;
    protected RealParameter insertionRate;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * constructor for testing purpose
     */
    public ScsFiniteMuDelInsModel() {
        super();
        initAndValidate();
    } // constructor

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // model 0, 1, 2, 3 alleles
        modeledAlleles = new int[]{0, 1, 2, 3};

        // fix the number of genotypes
        nrOfStates = 15;

        // fix the genotype of the tree root
        rootGenotype = 6;

        // fix the genotype of the constant site
        constGenotype = 6;

        // initialize the nrOfExistingAlleles
        nrOfExistingAlleles = new int[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            if (i <= 5) {
                nrOfExistingAlleles[i] = 3;
            } else if (i <= 11) {
                nrOfExistingAlleles[i] = 2;
            } else if (i <= 13) {
                nrOfExistingAlleles[i] = 1;
            } else {
                nrOfExistingAlleles[i] = 0;
            }
        }

        nrOfAltAlleles = new int[]{0, 1, 1, 2, 2, 3, 0, 0, 1, 1, 2, 2, 0, 1};

        ternaryCodes = new int[]{0, 1, 1, 1, 1, 2, 0, 0, 1, 1, 2, 2, 0, 2, 3};

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException
                | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)]; // maybe not of any use

        // parameter sanity check
        if (getParameter(deletionRateInput).getArrayValue() < 0.0
                || getParameter(insertionRateInput).getArrayValue() < 0.0) {
            throw new IllegalArgumentException("At least one of deletionRate and insertionRate is out of bound, " +
                    "which should not be smaller than 0 (" + this.getClass().getName() + ")");
        }
        deletionRate = getParameter(deletionRateInput);
        insertionRate = getParameter(insertionRateInput);
    } // initAndValidate

    @Override
    protected void setupRelativeRates() {
        /*
         * 0/00 0/01 1/00 0/11 1/01 1/11 0/0 -/00 0/1 -/01 1/1 -/11 0/- 1/- -
         */

        /* 0/00 */
        relativeRates[0] = 2.0 / 3; // -> 0/01
        relativeRates[1] = 1.0 / 3; // -> 1/00
        relativeRates[2] = 0.0; // -> 0/11
        relativeRates[3] = 0.0; // -> 1/01
        relativeRates[4] = 0.0; // -> 1/11
        relativeRates[5] = 2 * deletionRate.getArrayValue() / 3; // -> 0/0
        relativeRates[6] = deletionRate.getArrayValue() / 3; // -> -/00
        relativeRates[7] = 0.0; // -> 0/1
        relativeRates[8] = 0.0; // -> -/01
        relativeRates[9] = 0.0; // -> 1/1
        relativeRates[10] = 0.0; // -> -/11
        relativeRates[11] = 0.0; // -> 0/-
        relativeRates[12] = 0.0; // -> 1/-
        relativeRates[13] = 0.0; // -> -

        /* 0/01 */
        relativeRates[14] = 1.0 / 9; // -> 0/00
        relativeRates[15] = 0.0; // -> 1/00
        relativeRates[16] = 1.0 / 3; // -> 0/11
        relativeRates[17] = 1.0 / 3; // -> 1/01
        relativeRates[18] = 0.0; // -> 1/11
        relativeRates[19] = deletionRate.getArrayValue() / 3; // -> 0/0
        relativeRates[10] = 0.0; // -> -/00
        relativeRates[21] = deletionRate.getArrayValue() / 3; // -> 0/1
        relativeRates[22] = deletionRate.getArrayValue() / 3; // -> -/01
        relativeRates[23] = 0.0; // -> 1/1
        relativeRates[24] = 0.0; // -> -/11
        relativeRates[25] = 0.0; // -> 0/-
        relativeRates[26] = 0.0; // -> 1/-
        relativeRates[27] = 0.0; // -> -

        /* 1/00 */
        relativeRates[28] = 1.0 / 9; // -> 0/00
        relativeRates[29] = 0.0; // -> 0/01
        relativeRates[30] = 0.0; // -> 0/11
        relativeRates[31] = 2.0 / 3; // -> 1/01
        relativeRates[32] = 0.0; // -> 1/11
        relativeRates[33] = 0.0; // -> 0/0
        relativeRates[34] = deletionRate.getArrayValue() / 3; // -> -/00
        relativeRates[35] = 2 * deletionRate.getArrayValue() / 3; // -> 0/1
        relativeRates[36] = 0.0; // -> -/01
        relativeRates[37] = 0.0; // -> 1/1
        relativeRates[38] = 0.0; // -> -/11
        relativeRates[39] = 0.0; // -> 0/-
        relativeRates[40] = 0.0; // -> 1/-
        relativeRates[41] = 0.0; // -> -

        /* 0/11 */
        relativeRates[42] = 0.0; // -> 0/00
        relativeRates[43] = 2.0 / 9; // -> 0/01
        relativeRates[44] = 0.0; // -> 1/00
        relativeRates[45] = 0.0; // -> 1/01
        relativeRates[46] = 1.0 / 3; // -> 1/11
        relativeRates[47] = 0.0; // -> 0/0
        relativeRates[48] = 0.0; // -> -/00
        relativeRates[49] = 2 * deletionRate.getArrayValue() / 3; // -> 0/1
        relativeRates[50] = 0.0; // -> -/01
        relativeRates[51] = 0.0; // -> 1/1
        relativeRates[52] = deletionRate.getArrayValue() / 3; // -> -/11
        relativeRates[53] = 0.0; // -> 0/-
        relativeRates[54] = 0.0; // -> 1/-
        relativeRates[55] = 0.0; // -> -

        /* 1/01 */
        relativeRates[56] = 0.0; // -> 0/00
        relativeRates[57] = 1.0 / 9; // -> 0/01
        relativeRates[58] = 1.0 / 9; // -> 1/00
        relativeRates[59] = 0.0; // -> 0/11
        relativeRates[60] = 1.0 / 3; // -> 1/11
        relativeRates[61] = 0.0; // -> 0/0
        relativeRates[62] = 0.0; // -> -/00
        relativeRates[63] = deletionRate.getArrayValue() / 3; // -> 0/1
        relativeRates[64] = deletionRate.getArrayValue() / 3; // -> -/01
        relativeRates[65] = deletionRate.getArrayValue() / 3; // -> 1/1
        relativeRates[66] = 0.0; // -> -/11
        relativeRates[67] = 0.0; // -> 0/-
        relativeRates[68] = 0.0; // -> 1/-
        relativeRates[69] = 0.0; // -> -

        /* 1/11 */
        relativeRates[70] = 0.0; // -> 0/00
        relativeRates[71] = 0.0; // -> 0/01
        relativeRates[72] = 0.0; // -> 1/00
        relativeRates[73] = 1.0 / 9; // -> 0/11
        relativeRates[74] = 2.0 / 9; // -> 1/01
        relativeRates[75] = 0.0; // -> 0/0
        relativeRates[76] = 0.0; // -> -/00
        relativeRates[77] = 0.0; // -> 0/1
        relativeRates[78] = 0.0; // -> -/01
        relativeRates[79] = 2 * deletionRate.getArrayValue() / 3; // -> 1/1
        relativeRates[80] = deletionRate.getArrayValue() / 3; // -> -/11
        relativeRates[81] = 0.0; // -> 0/-
        relativeRates[82] = 0.0; // -> 1/-
        relativeRates[83] = 0.0; // -> -

        /* 0/0 */
        relativeRates[84] = 2.0 * insertionRate.getArrayValue() / 3; // -> 0/00
        relativeRates[85] = 0.0; // -> 0/01
        relativeRates[86] = 0.0; // -> 1/00
        relativeRates[87] = 0.0; // -> 0/11
        relativeRates[88] = 0.0; // -> 1/01
        relativeRates[89] = 0.0; // -> 1/11
        relativeRates[90] = 0.0; // -> -/00
        relativeRates[91] = 2.0 / 3; // -> 0/1
        relativeRates[92] = 0.0; // -> -/01
        relativeRates[93] = 0.0; // -> 1/1
        relativeRates[94] = 0.0; // -> -/11
        relativeRates[95] = 2.0 * deletionRate.getArrayValue() / 3; // -> 0/-
        relativeRates[96] = 0.0; // -> 1/-
        relativeRates[97] = 0.0; // -> -

        /* -/00 */
        relativeRates[98] = 0.0; // -> 0/00
        relativeRates[99] = 0.0; // -> 0/01
        relativeRates[100] = 0.0; // -> 1/00
        relativeRates[101] = 0.0; // -> 0/11
        relativeRates[102] = 0.0; // -> 1/01
        relativeRates[103] = 0.0; // -> 1/11
        relativeRates[104] = 0.0; // -> 0/0
        relativeRates[105] = 0.0; // -> 0/1
        relativeRates[106] = 2.0 / 3; // -> -/01
        relativeRates[107] = 0.0; // -> 1/1
        relativeRates[108] = 0.0; // -> -/11
        relativeRates[109] = deletionRate.getArrayValue() / 3; // -> 0/-
        relativeRates[110] = 0.0; // -> 1/-
        relativeRates[111] = 0.0; // -> -

        /* 0/1 */
        relativeRates[112] = 0.0; // -> 0/00
        relativeRates[113] = 0.0; // -> 0/01
        relativeRates[114] = insertionRate.getArrayValue() / 3; // -> 1/00
        relativeRates[115] = insertionRate.getArrayValue() / 3; // -> 0/11
        relativeRates[116] = 0.0; // -> 1/01
        relativeRates[117] = 0.0; // -> 1/11
        relativeRates[118] = 1.0 / 9; // -> 0/0
        relativeRates[119] = 0.0; // -> -/00
        relativeRates[120] = 0.0; // -> -/01
        relativeRates[121] = 1.0 / 3; // -> 1/1
        relativeRates[122] = 0.0; // -> -/11
        relativeRates[123] = deletionRate.getArrayValue() / 3; // -> 0/-
        relativeRates[124] = deletionRate.getArrayValue() / 3; // -> 1/-
        relativeRates[125] = 0.0; // -> -

        /* -/01 */
        relativeRates[126] = 0.0; // -> 0/00
        relativeRates[127] = 0.0; // -> 0/01
        relativeRates[128] = 0.0; // -> 1/00
        relativeRates[129] = 0.0; // -> 0/11
        relativeRates[130] = 0.0; // -> 1/01
        relativeRates[131] = 0.0; // -> 1/11
        relativeRates[132] = 0.0; // -> 0/0
        relativeRates[133] = 1.0 / 9; // -> -/00
        relativeRates[134] = 0.0; // -> 0/1
        relativeRates[135] = 0.0; // -> 1/1
        relativeRates[136] = 1.0 / 3; // -> -/11
        relativeRates[137] = deletionRate.getArrayValue() / 3; // -> 0/-
        relativeRates[138] = deletionRate.getArrayValue() / 3; // -> 1/-
        relativeRates[139] = 0.0; // -> -

        /* 1/1 */
        relativeRates[140] = 0.0; // -> 0/00
        relativeRates[141] = 0.0; // -> 0/01
        relativeRates[142] = 0.0; // -> 1/00
        relativeRates[143] = 0.0; // -> 0/11
        relativeRates[144] = 0.0; // -> 1/01
        relativeRates[145] = 2 * insertionRate.getArrayValue() / 3; // -> 1/11
        relativeRates[146] = 0.0; // -> 0/0
        relativeRates[147] = 0.0; // -> -/00
        relativeRates[148] = 2.0 / 9; // -> 0/1
        relativeRates[149] = 0.0; // -> -/01
        relativeRates[150] = 0.0; // -> -/11
        relativeRates[151] = 0.0; // -> 0/-
        relativeRates[152] = 2 * deletionRate.getArrayValue() / 3; // -> 1/-
        relativeRates[153] = 0.0; // -> -

        /* -/11 */
        relativeRates[154] = 0.0; // -> 0/00
        relativeRates[155] = 0.0; // -> 0/01
        relativeRates[156] = 0.0; // -> 1/00
        relativeRates[157] = 0.0; // -> 0/11
        relativeRates[158] = 0.0; // -> 1/01
        relativeRates[159] = 0.0; // -> 1/11
        relativeRates[160] = 0.0; // -> 0/0
        relativeRates[161] = 0.0; // -> -/00
        relativeRates[162] = 0.0; // -> 0/1
        relativeRates[163] = 2.0 / 9; // -> -/01
        relativeRates[164] = 0.0; // -> 1/1
        relativeRates[165] = 0.0; // -> 0/-
        relativeRates[166] = 2 * deletionRate.getArrayValue() / 3; // -> 1/-
        relativeRates[167] = 0.0; // -> -

        /* 0/- */
        relativeRates[168] = 0.0; // -> 0/00
        relativeRates[169] = 0.0; // -> 0/01
        relativeRates[170] = 0.0; // -> 1/00
        relativeRates[171] = 0.0; // -> 0/11
        relativeRates[172] = 0.0; // -> 1/01
        relativeRates[173] = 0.0; // -> 1/11
        relativeRates[174] = 0.0; // -> 0/0
        relativeRates[175] = insertionRate.getArrayValue() / 3; // -> -/00
        relativeRates[176] = 0.0; // -> 0/1
        relativeRates[177] = 0.0; // -> -/01
        relativeRates[178] = 0.0; // -> 1/1
        relativeRates[179] = 0.0; // -> -/11
        relativeRates[180] = 1.0 / 3; // -> 1/-
        relativeRates[181] = deletionRate.getArrayValue() / 3; // -> -

        /* 1/- */
        relativeRates[182] = 0.0; // -> 0/00
        relativeRates[183] = 0.0; // -> 0/01
        relativeRates[184] = 0.0; // -> 1/00
        relativeRates[185] = 0.0; // -> 0/11
        relativeRates[186] = 0.0; // -> 1/01
        relativeRates[187] = 0.0; // -> 1/11
        relativeRates[188] = 0.0; // -> 0/0
        relativeRates[189] = 0.0; // -> -/00
        relativeRates[190] = 0.0; // -> 0/1
        relativeRates[191] = 0.0; // -> -/01
        relativeRates[192] = 0.0; // -> 1/1
        relativeRates[193] = insertionRate.getArrayValue() / 3; // -> -/11
        relativeRates[194] = 1.0 / 9; // -> 0/-
        relativeRates[195] = deletionRate.getArrayValue() / 3; // -> -

        /* - */
        relativeRates[196] = 0.0; // -> 0/00
        relativeRates[197] = 0.0; // -> 0/01
        relativeRates[198] = 0.0; // -> 1/00
        relativeRates[199] = 0.0; // -> 0/11
        relativeRates[200] = 0.0; // -> 1/01
        relativeRates[201] = 0.0; // -> 1/11
        relativeRates[202] = 0.0; // -> 0/0
        relativeRates[203] = 0.0; // -> -/00
        relativeRates[204] = 0.0; // -> 0/1
        relativeRates[205] = 0.0; // -> -/01
        relativeRates[206] = 0.0; // -> 1/1
        relativeRates[207] = 0.0; // -> -/11
        relativeRates[208] = 0.0; // -> 0/-
        relativeRates[209] = 0.0; // -> 1/-
    } // setupRelativeRates

    @Override
    public String getAlphabetGenotype(final int index) {
        if (index < 0 || index > nrOfStates - 1) {
            throw new IllegalArgumentException("Index exceeds the boundary (0 - " + (nrOfStates - 1) + ")");
        }

        if (index == 0) {
            return "0/00";
        } else if (index == 1) {
            return "0/01";
        } else if (index == 2) {
            return "1/00";
        } else if (index == 3) {
            return "0/11";
        } else if (index == 4) {
            return "1/01";
        } else if (index == 5) {
            return "1/11";
        } else if (index == 6) {
            return "0/0";
        } else if (index == 7) {
            return "./00";
        } else if (index == 8) {
            return "0/1";
        } else if (index == 9) {
            return "./01";
        } else if (index == 10) {
            return "1/1";
        } else if (index == 11) {
            return "./11";
        } else if (index == 12) {
            return "./0";
        } else if (index == 13) {
            return "./1";
        } else {
            return "./.";
        }
    } // getAlphabetGenotype

    @Override
    public boolean canHandleDataType(ReadCounts readCountsDataType) {
        return readCountsDataType instanceof CovSup || readCountsDataType instanceof FullSupsCov;
    } // canHandleDataType


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    public RealParameter getDeletionRate() {
        return deletionRate;
    } // getDeletionRate

    public RealParameter getInsertionRate() {
        return insertionRate;
    } // getInsertionRate

    public int[] getNrOfExistingAlleles() {
        return nrOfExistingAlleles;
    } // getNrOfAlleles

    @Override
    public int getNrOfAlleles(final int genotypeIndex) {
        return nrOfExistingAlleles[genotypeIndex];
    } // getNrOfAlleles

    /**
     * Get the evolutionary events during the process of a parent evolving to a child.
     *
     * @param parentGenotype genotype of parent
     * @param childGenotype  genotype of a child
     * @return an array of {@link EvolutionaryEventType}
     */
    @Override
    public EvolutionaryEventType[] getEvolutionaryEvents(int parentGenotype, int childGenotype) {
        // Some very special cases
        if (parentGenotype == 8) {
            if (childGenotype == 7)
                return new EvolutionaryEventType[]{EvolutionaryEventType.DELETION, EvolutionaryEventType.INSERTION};

            if (childGenotype == 9)
                return new EvolutionaryEventType[]{EvolutionaryEventType.DELETION, EvolutionaryEventType.INSERTION, EvolutionaryEventType.SINGLE_MUTATION};
        }

        return super.getEvolutionaryEvents(parentGenotype, childGenotype);
    } // getEvolutionaryEvents

} // ScsFiniteMuDelInsModel
