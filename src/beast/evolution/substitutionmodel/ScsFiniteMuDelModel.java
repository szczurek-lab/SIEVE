package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;

import java.lang.reflect.InvocationTargetException;

@Description(value = "A finite mutation and deletion substitution model with 6 genotypes")
public class ScsFiniteMuDelModel extends ScsSubstitutionModelBase {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public final Input<RealParameter> deletionRateInput = new Input<>("deletionRate", "deletion rate " +
            "(non-negative; default 0); measured relatively to mutation rate", Input.Validate.REQUIRED);

    protected RealParameter deletionRate;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * constructor for testing purpose
     */
    public ScsFiniteMuDelModel() {
        super();
        initAndValidate();
    } // constructor

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // model 0, 1, 2 alleles
        modeledAlleles = new int[]{0, 1, 2};

        // fix the number of genotypes
        nrOfStates = 6;

        // fix the genotype of the tree root
        rootGenotype = 0;

        // fix the genotype of the constant site
        constGenotype = 0;

        // initialize the nrOfExistingAlleles
        nrOfExistingAlleles = new int[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            if (i <= 2) {
                nrOfExistingAlleles[i] = 2;
            } else if (i <= 4) {
                nrOfExistingAlleles[i] = 1;
            } else {
                nrOfExistingAlleles[i] = 0;
            }
        }

        nrOfAltAlleles = new int[]{0, 1, 2, 0, 1, 0};

        ternaryCodes = new int[]{0, 1, 2, 0, 2, 3};

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
        if (getParameter(deletionRateInput).getArrayValue() < 0.0) {
            throw new IllegalArgumentException("deletionRate is out of bound, which should not be smaller than 0 (" +
                    this.getClass().getName() + ")");
        }
        deletionRate = getParameter(deletionRateInput);
    } // initAndValidate

    @Override
    protected void setupRelativeRates() {
        /*
         * 0/0 0/1 1/1 0/- 1/- -
         */

        /* 0/0 */
        relativeRates[0] = 1.0; // -> 0/1
        relativeRates[1] = 0.0; // -> 1/1
        relativeRates[2] = deletionRate.getArrayValue(); // -> 0/-
        relativeRates[3] = 0.0; // -> 1/-
        relativeRates[4] = 0.0; // -> -

        /* 0/1 */
        relativeRates[5] = 1.0 / 6; // -> 0/0
        relativeRates[6] = 1.0 / 2; // -> 1/1
        relativeRates[7] = deletionRate.getArrayValue() / 2; // -> 0/-
        relativeRates[8] = deletionRate.getArrayValue() / 2; // -> 1/-
        relativeRates[9] = 0.0; // -> -

        /* 1/1 */
        relativeRates[10] = 0.0; // -> 0/0
        relativeRates[11] = 1.0 / 3; // -> 0/1
        relativeRates[12] = 0.0; // -> 0/-
        relativeRates[13] = deletionRate.getArrayValue(); // -> 1/-
        relativeRates[14] = 0.0; // -> -

        /* 0/- */
        relativeRates[15] = 0.0; // -> 0/0
        relativeRates[16] = 0.0; // -> 0/1
        relativeRates[17] = 0.0; // -> 1/1
        relativeRates[18] = 1.0 / 2; // -> 1/-
        relativeRates[19] = deletionRate.getArrayValue() / 2; // -> -

        /* 1/- */
        relativeRates[20] = 0.0; // -> 0/0
        relativeRates[21] = 0.0; // -> 0/1
        relativeRates[22] = 0.0; // -> 1/1
        relativeRates[23] = 1.0 / 6; // -> 0/-
        relativeRates[24] = deletionRate.getArrayValue() / 2; // -> -

        /* - */
        relativeRates[25] = 0.0; // -> 0/0
        relativeRates[26] = 0.0; // -> 0/1
        relativeRates[27] = 0.0; // -> 1/1
        relativeRates[28] = 0.0; // -> 0/-
        relativeRates[29] = 0.0; // -> 1/-
    } // setupRelativeRates

    @Override
    public String getAlphabetGenotype(final int index) {
        if (index < 0 || index > nrOfStates - 1) {
            throw new IllegalArgumentException("Index exceeds the boundary (0 - " + (nrOfStates - 1) + ")");
        }

        if (index == 0) {
            return "0/0";
        } else if (index == 1) {
            return "0/1";
        } else if (index == 2) {
            return "1/1";
        } else if (index == 3) {
            return "0/-";
        } else if (index == 4) {
            return "1/-";
        } else {
            return "-/-";
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

    public int[] getNrOfExistingAlleles() {
        return nrOfExistingAlleles;
    } // getNrOfAlleles

    @Override
    public int getNrOfAlleles(final int genotypeIndex) {
        return nrOfExistingAlleles[genotypeIndex];
    } // getNrOfAlleles

} // ScsFiniteMuDelModel
