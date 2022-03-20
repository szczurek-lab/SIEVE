package beast.evolution.substitutionmodel;

import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;
import beast.evolution.variantsinfo.vcfentry.CandidateAltNuc;

import java.lang.reflect.InvocationTargetException;
import java.util.List;

public class ScsFiniteMuExtendedModel extends ScsSubstitutionModelBase {

    /**
     * constructor for testing purpose
     */
    public ScsFiniteMuExtendedModel() {
        super();
        initAndValidate();
    } // constructor


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // only 2 alleles modeled
        modeledAlleles = new int[]{2};

        // fix the number of genotypes
        nrOfStates = 4;

        // fix the genotype of the tree root
        rootGenotype = 0;

        // fix the genotype of the constant site
        constGenotype = 0;

        // initialize the nrOfExistingAlleles
        nrOfExistingAlleles = new int[1];
        nrOfExistingAlleles[0] = 2;

        nrOfAltAlleles = new int[]{0, 1, 2, 2};

        ternaryCodes = new int[]{0, 1, 2, 2};

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
         * 0/0 0/1 1/1 1/2
         */

        /* 0/0 */
        relativeRates[0] = 1.0; // -> 0/1
        relativeRates[1] = 0.0; // -> 1/1
        relativeRates[2] = 0.0; // -> 1/2

        /* 0/1 */
        relativeRates[3] = 1.0 / 6; // -> 0/0
        relativeRates[4] = 1.0 / 6; // -> 1/1
        relativeRates[5] = 1.0 / 3; // -> 1/2

        /* 1/1 */
        relativeRates[6] = 0.0; // -> 0/0
        relativeRates[7] = 1.0 / 3; // -> 0/1
        relativeRates[8] = 2.0 / 3; // -> 1/2

        /* 1/2 */
        relativeRates[9] = 0.0; // -> 0/0
        relativeRates[10] = 1.0 / 3; // -> 0/1
        relativeRates[11] = 1.0 / 3; // -> 1/1
    } // setupRelativeRates

    @Override
    public String getAlphabetGenotype(int index) {
        if (index < 0 || index > nrOfStates - 1) {
            throw new IllegalArgumentException("Index exceeds the boundary (0 - " + (nrOfStates - 1) + ")");
        }

        if (index == 0) {
            return "0/0";
        } else if (index == 1) {
            return "0/1";
        } else if (index == 2) {
            return "1/1";
        } else {
            return "1/2";
        }
    } // getAlphabetGenotype

    /**
     * Adjust an allele.
     * <p/>
     * This function only handles these characters: '0', '1', and '2'.
     * <p/>
     * Any substitution model not matching the above rule should override this function.
     *
     * @param allele          allele to be adjusted
     * @param locusAltNucs    locus-wise alternative nucleotides
     * @param cellAltNucs     cell-wise alternative nucleotides in descending order according to read counts
     * @param missingAllele   character to represent a missing allele
     * @param existingAltNucs existing alt nucs in a cell
     * @return an adjusted allele
     */
    @Override
    protected char adaptAllele(
            char allele,
            List<CandidateAltNuc> locusAltNucs,
            char[] cellAltNucs,
            char missingAllele,
            List<Character> existingAltNucs
    ) {
        switch (allele) {
            case '0':
                return '0';
            case '1':
                return (cellAltNucs == null || cellAltNucs.length == 0) ?
                        adaptAltAllele(locusAltNucs, 'N', existingAltNucs) :
                        adaptAltAllele(locusAltNucs, cellAltNucs[0], existingAltNucs);
            case '2':
                return (cellAltNucs == null || cellAltNucs.length == 0) ?
                        adaptAltAllele(locusAltNucs, 'N', existingAltNucs) :
                        adaptAltAllele(locusAltNucs, cellAltNucs[1], existingAltNucs);
            default:
                throw new IllegalArgumentException("Error! Unsupported character: " + allele +
                        ". Only '0', '1', and '2' are allowed.");
        }
    } // adaptAllele

    @Override
    public boolean canHandleDataType(ReadCounts readCountsDataType) {
        return readCountsDataType instanceof CovSup || readCountsDataType instanceof FullSupsCov;
    } // canHandleDataType

    @Override
    public EvolutionaryEventType[] getEvolutionaryEvents(int parentGenotype, int childGenotype) {
        switch (parentGenotype) {
            case 0:
                switch (childGenotype) {
                    case 0:
                        return null;
                    case 1:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_MUTATION};
                    case 2:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HOMO_SIMU_DOUBLE_MUTATION};
                    case 3:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HETERO_SIMU_DOUBLE_MUTATION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 1:
                switch (childGenotype) {
                    case 0:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_BACK_MUTATION};
                    case 1:
                        return null;
                    case 2:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HOMO_SINGLE_MUTATION_ADDITION};
                    case 3:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HETERO_SINGLE_MUTATION_ADDITION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 2:
                switch (childGenotype) {
                    case 0:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.DOUBLE_BACK_MUTATION};
                    case 1:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_BACK_MUTATION};
                    case 2:
                        return null;
                    case 3:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HETERO_SUBST_SINGLE_MUTATION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 3:
                switch (childGenotype) {
                    case 0:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.DOUBLE_BACK_MUTATION};
                    case 1:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_BACK_MUTATION};
                    case 2:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HOMO_SUBST_SINGLE_MUTATION};
                    case 3:
                        return null;
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            default:
                throw new IllegalArgumentException("Error! Unsupported genotype for the parent: " + parentGenotype);
        }
    } // getEvolutionaryEvents

    /**
     * Match an evolutionary event. Only supporting '0', '1' and '2'.
     *
     * @param from from which allele?
     * @param to   to which allele?
     * @return matched evolutionary event
     */
    @Override
    @Deprecated
    protected EvolutionaryEventType matchEvolutionaryEvent(char from, char to) {
        if (!Character.isDigit(from) || !Character.isDigit(to))
            throw new IllegalStateException(INVALID_TRANSITION);

        if (from == to)
            return null;

        final String query = from + String.valueOf(to);
        switch (query) {
            case "01":
            case "02":
            case "12":
            case "21":
                return EvolutionaryEventType.SINGLE_MUTATION;
            case "10":
            case "20":
                return EvolutionaryEventType.SINGLE_BACK_MUTATION;
            default:
                throw new IllegalStateException("Unexpected value: " + query);
        }
    }

    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    @Override
    public int getNrOfAlleles(int genotypeIndex) {
        return nrOfExistingAlleles[0];
    } // getNrOfAlleles

}
