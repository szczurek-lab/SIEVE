package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;

import java.lang.reflect.InvocationTargetException;

@Description(value = "A finite mutation substitution model with 3 genotypes")
public class ScsFiniteMuModel extends ScsSubstitutionModelBase {


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * constructor for testing purpose
     */
    public ScsFiniteMuModel() {
        super();
        initAndValidate();
    } // constructor

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // only 2 alleles modeled
        modeledAlleles = new int[]{2};

        // fix the number of genotypes
        nrOfStates = 3;

        // fix the genotype of the tree root
        rootGenotype = 0;

        // fix the genotype of the constant site
        constGenotype = 0;

        // initialize the nrOfExistingAlleles
        nrOfExistingAlleles = new int[1];
        nrOfExistingAlleles[0] = 2;

        nrOfAltAlleles = new int[]{0, 1, 2};

        ternaryCodes = new int[]{0, 1, 2};

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException
                | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)]; // maybe not of any use
    } // initAndValidate

    @Override
    protected void setupRelativeRates() {
        /*
         * 0/0 0/1 1/1
         */

        /* 0/0 */
        relativeRates[0] = 1.0; // -> 0/1
        relativeRates[1] = 0.0; // -> 1/1

        /* 0/1 */
        relativeRates[2] = 1.0 / 6; // -> 0/0
        relativeRates[3] = 1.0 / 2; // -> 1/1

        /* 1/1 */
        relativeRates[4] = 0.0; // -> 0/0
        relativeRates[5] = 1.0 / 3; // -> 0/1
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
        } else {
            return "1/1";
        }
    } // getAlphabetGenotype

    @Override
    public boolean canHandleDataType(ReadCounts readCountsDataType) {
        return readCountsDataType instanceof CovSup || readCountsDataType instanceof FullSupsCov;
    } // canHandleDataType


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    @Override
    public int getNrOfAlleles(final int genotypeIndex) {
        return nrOfExistingAlleles[0];
    } // getNrOfAlleles

} // ScsFiniteMuModel
