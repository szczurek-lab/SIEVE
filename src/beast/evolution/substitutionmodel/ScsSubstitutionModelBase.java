package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.ReadCounts;
import beast.evolution.variantsinfo.vcfentry.CandidateAltNuc;
import com.google.common.primitives.Chars;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Each extended class of ScsSubstitutionModelBase should implement:
 * canHandleDataType(ReadCounts readCountsDataType)
 * setupRelativeRates()
 * setupRateMatrix()
 */
@Description(value = "Base implementation of a substitution model for single cell DNA sequencing read counts data")
public abstract class ScsSubstitutionModelBase extends GeneralSubstitutionModel {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    protected final static String INVALID_TRANSITION = "Transition of genotypes from the first to the second is illegal: ";

    public enum EvolutionaryEventType {

        // single mutation (0/0 -> 0/1)
        SINGLE_MUTATION("SM"),

        // homozygous simultaneous double mutation (0/0 -> 1/1)
        HOMO_SIMU_DOUBLE_MUTATION("HoSDM"),

        // heterozygous simultaneous double mutation (0/0 -> 1/1')
        HETERO_SIMU_DOUBLE_MUTATION("HeSDM"),

        // single back mutation (0/1 -> 0/0)
        SINGLE_BACK_MUTATION("SB"),

        // double back mutation (1/1 -> 0/0)
        DOUBLE_BACK_MUTATION("DB"),

        // homozygous single mutation addition (0/1 -> 1/1)
        HOMO_SINGLE_MUTATION_ADDITION("HoSMA"),

        // heterozygous single mutation addition (0/1 -> 1/1')
        HETERO_SINGLE_MUTATION_ADDITION("HeSMA"),

        // homozygous substitute single mutation (1/1' -> 1/1)
        HOMO_SUBST_SINGLE_MUTATION("HoSSM"),

        // heterozygous substitute single mutation (1/1 -> 1/1')
        HETERO_SUBST_SINGLE_MUTATION("HeSSM"),

        DELETION("D"),
        INSERTION("I");

        String desc;

        EvolutionaryEventType(String s) {
            desc = s;
        }

        @Override
        public String toString() {
            return desc;
        }

        public static String arrayToString(final EvolutionaryEventType[] arr) {
            return arrayToString(arr, "");
        } // arrayToString

        public static String arrayToString(final EvolutionaryEventType[] arr, final String separator) {
            StringBuilder stringBuilder = new StringBuilder();

            for (int index = 0; index < arr.length; index++) {
                stringBuilder.append(arr[index].toString());

                if (index < arr.length - 1)
                    stringBuilder.append(separator);
            }

            return stringBuilder.toString();
        } // arrayToString

        /**
         * Whether the evolutionary event violates infinite-sites assumption or not?
         *
         * @return yes or no
         */
        public boolean violateISA() {
            return this != EvolutionaryEventType.SINGLE_MUTATION;
        }

        public int compare(@NotNull EvolutionaryEventType e) {
            return this.desc.compareTo(e.desc);
        }

        public static class EventsComparator implements Comparator<EvolutionaryEventType> {

            @Override
            public int compare(EvolutionaryEventType e1, EvolutionaryEventType e2) {
                return e1.compare(e2);
            }

        }

    }

    /**
     * genotype of the root
     */
    protected int rootGenotype;

    /**
     * genotype of the constant site
     */
    protected int constGenotype;

    /**
     * the number of alternative alleles each genotype contains
     */
    protected int[] nrOfAltAlleles;

    /**
     * corresponding number of existing alleles for each genotype
     */
    protected int[] nrOfExistingAlleles;

    /**
     * the number of existing alleles this model studies
     */
    protected int[] modeledAlleles;

    /**
     * ternary codes corresponding to each genotype
     * 0: homogeneous reference
     * 1: heterogeneous alternative
     * 2: homogeneous alternative
     * 3: missing data
     */
    protected int[] ternaryCodes;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * constructor for testing purpose
     */
    public ScsSubstitutionModelBase() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        try {
            ratesInput.setValue(null, this);
        } catch (Exception e) {
            throw new IllegalArgumentException(e.getMessage());
        }

        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        try {
            frequenciesInput.setValue(null, this);
        } catch (Exception e) {
            throw new IllegalArgumentException(e.getMessage());
        }
    } // constructor

    /**
     * ratesInput and frequenciesInput should not be specified in the configuration xml document as
     * they are not applied to the derived classes of ScsSubstitutionModelBase
     */
    @Override
    public void initAndValidate() {
        if (ratesInput.get() != null) {
            throw new IllegalArgumentException("the rates attribute should not be used for the selected substitution model (" + this.getClass().getName() + ")");
        }
        if (frequenciesInput.get() != null) {
            throw new IllegalArgumentException("the frequencies attribute should not be used for the selected substitution model (" + this.getClass().getName() + ")");
        }
    } // initAndValidate

    /**
     * @param parameterInput a parameter of type Input<Function>
     * @return either the input value or default 0.0
     */
    protected RealParameter getParameter(Input<RealParameter> parameterInput) {
        if (parameterInput.get() != null) {
            return parameterInput.get();
        }
        return new RealParameter("0.0");
    } // getParameter

    @Override
    protected void setupRateMatrix() {
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0.0;
            for (int j = 0; j < i; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j];
            }
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j - 1];
            }
        }

        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double sum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    sum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -sum;
        }

        // normalise rate matrix to one expected substitution per unit time
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            subst += -rateMatrix[i][i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / subst;
            }
        }
    } // setupRateMatrix

    @Override
    protected abstract void setupRelativeRates();

    public abstract String getAlphabetGenotype(final int index);

    public abstract boolean canHandleDataType(ReadCounts readCountsDataType);

    /**
     * Get genotype compatible with VCF standard.
     *
     * @param genotype      genotype code
     * @param locusAltNucs  locus-wise alternative nucleotides
     * @param cellAltNucs   cell-wise alternative nucleotides in descending order according to read counts
     * @param missingAllele character to represent a missing allele
     * @return adjusted genotype
     */
    public String getGenotypeForVCF(
            int genotype,
            List<CandidateAltNuc> locusAltNucs,
            char[] cellAltNucs,
            char missingAllele
    ) {
        List<Character> existingAltNucs = new ArrayList<>();

        final String[] allelesOri = getAlphabetGenotype(genotype).split("/");

        String chrom1Str = String.copyValueOf(
                adaptAllelesChrom(
                        allelesOri[0],
                        locusAltNucs,
                        cellAltNucs,
                        missingAllele,
                        existingAltNucs
                )
        );
        String chrom2Str = String.copyValueOf(
                adaptAllelesChrom(
                        allelesOri[1],
                        locusAltNucs,
                        cellAltNucs,
                        missingAllele,
                        existingAltNucs
                )
        );

        return String.join("/", new String[]{chrom1Str, chrom2Str});
    } // getGenotype

    /**
     * Get alleles compatible with VCF standard per chromosome.
     *
     * @param alleles         alleles on a strand of chromosome
     * @param locusAltNucs    locus-wise alternative nucleotides
     * @param cellAltNucs     cell-wise alternative nucleotides in descending order according to read counts
     * @param missingAllele   character to represent a missing allele
     * @param existingAltNucs existing alt nucs in a cell
     * @return adjusted alleles
     */
    protected char[] adaptAllelesChrom(
            String alleles,
            List<CandidateAltNuc> locusAltNucs,
            char[] cellAltNucs,
            char missingAllele,
            List<Character> existingAltNucs
    ) {
        List<Character> results = new ArrayList<>();

        for (int i = 0; i < alleles.length(); i++) {
            results.add(
                    adaptAllele(
                            alleles.charAt(i),
                            locusAltNucs,
                            cellAltNucs,
                            missingAllele,
                            existingAltNucs
                    )
            );
        }

        return Chars.toArray(results);
    } // adaptAllelesChrom

    /**
     * Adjust an allele.
     * <p/>
     * This function only handles these characters: '0', '1', and '-'.
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
            case '-':
                return missingAllele;
            case '1':
                return (cellAltNucs == null || cellAltNucs.length == 0) ?
                        adaptAltAllele(locusAltNucs, 'N', existingAltNucs) :
                        adaptAltAllele(locusAltNucs, cellAltNucs[0], existingAltNucs);
            default:
                throw new IllegalArgumentException("Error! Unsupported character: " + allele +
                        ". Only '0', '1', and '-' are allowed.");
        }
    } // adaptAllele

    /**
     * Adjust the alternative allele: {@param cellAltNuc}.
     * <p/>
     * e.g. {@param allele} is '1', {@param cellAltNuc} is 'C', {@param locusAltNucs} is {'T', 'C', 'G'},
     * then the index of 'C' in {@param locusAltNucs} plus 1 is returned, i.e., 2.
     *
     * @param locusAltNucs    locus-wise alternative nucleotides
     * @param cellAltNuc      cell-wise alternative nucleotides in descending order according to read counts
     * @param existingAltNucs existing alt nucs in a cell
     * @return an adjusted allele
     */
    protected char adaptAltAllele(
            List<CandidateAltNuc> locusAltNucs,
            char cellAltNuc,
            List<Character> existingAltNucs
    ) {
        if (Character.toUpperCase(cellAltNuc) == 'N')
            return '1';

        for (int i = 0; i < locusAltNucs.size(); i++) {
            if (Character.toUpperCase(locusAltNucs.get(i).getNuc()) == Character.toUpperCase(cellAltNuc)) {
                locusAltNucs.get(i).addCount(1);
                if (!existingAltNucs.contains(cellAltNuc)) {
                    existingAltNucs.add(cellAltNuc);
                    locusAltNucs.get(i).addNumOfCells(1);
                }

                return Character.forDigit(i + 1, 10);
            }
        }

        locusAltNucs.add(
                new CandidateAltNuc(cellAltNuc)
        );
        locusAltNucs.get(locusAltNucs.size() - 1).addCount(1);
        existingAltNucs.add(cellAltNuc);
        locusAltNucs.get(locusAltNucs.size() - 1).addNumOfCells(1);
        return Character.forDigit(locusAltNucs.size(), 10);
    } // adaptAltAllele

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return false;
    }

    public int getRootGenotype() {
        return rootGenotype;
    } // getRootGenotype

    public int getConstGenotype() {
        return constGenotype;
    } // getConstGenotype

    public abstract int getNrOfAlleles(final int genotypeIndex); // getNrOfAlleles

    public int getModeledAllelesSize() {
        return modeledAlleles.length;
    } // getModeledAllelesSize

    public int[] getModeledAlleles() {
        return modeledAlleles;
    } // getModeledAlleles

    public void getModeledAlleles(int[] out) {
        System.arraycopy(modeledAlleles, 0, out, 0, modeledAlleles.length);
    } // getModeledAlleles

    public String getAllGenotypes(String delimiter) {
        String[] genotypes = new String[this.nrOfStates];

        for (int i = 0; i < nrOfStates; i++) {
            genotypes[i] = getAlphabetGenotype(i);
        }

        return String.join(delimiter, genotypes);
    } // getAllGenotypes

    public int getNrOfAltAlleles(int genotypeIndex) {
        return this.nrOfAltAlleles[genotypeIndex];
    } // getNrOfAltAlleles

    public int getTernaryCode(int genotypeIndex) {
        return this.ternaryCodes[genotypeIndex];
    } // getTernaryCode

    public boolean isVariant(int genotypeIndex) {
        return genotypeIndex != this.rootGenotype;
    } // isVariant

    /**
     * Get the evolutionary events during the process of a parent evolving to a child.
     *
     * @param parentGenotype genotype of parent
     * @param childGenotype  genotype of a child
     * @return an array of {@link EvolutionaryEventType}
     */
    public EvolutionaryEventType[] getEvolutionaryEvents(final int parentGenotype, final int childGenotype) {
        if (parentGenotype == childGenotype)
            return null;

        final String pGT = getAlphabetGenotype(parentGenotype);
        final String cGT = getAlphabetGenotype(childGenotype);

        final String[] pAlleles = pGT.split("/");
        final String[] cAlleles = cGT.split("/");

        assert pAlleles.length == 2;
        assert cAlleles.length == 2;

        List<EvolutionaryEventType> events = new ArrayList<>();

        if (pGT.length() > 3 && cGT.length() > 3) {
            // MuDelIns model with insertions for both of parent and child

            List<EvolutionaryEventType> eventsChr1 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[0], cAlleles[0]));
            List<EvolutionaryEventType> eventsChr2 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[1], cAlleles[1]));

            if (eventsChr1 != null)
                events.addAll(eventsChr1);

            if (eventsChr2 != null)
                events.addAll(eventsChr2);
        } else if (pGT.length() > 3) {
            // MuDelIns model with insertions for with parent

            List<EvolutionaryEventType> eventsChr1;
            List<EvolutionaryEventType> eventsChr2;

            if (pGT.charAt(0) == cGT.charAt(0) || (pGT.charAt(0) != cGT.charAt(0) && pGT.charAt(0) != cGT.charAt(2))) {
                eventsChr1 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[0], cAlleles[0]));
                eventsChr2 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[1], cAlleles[1]));
            } else {
                eventsChr1 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[0], cAlleles[1]));
                eventsChr2 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[1], cAlleles[0]));
            }

            if (eventsChr1 != null)
                events.addAll(eventsChr1);

            if (eventsChr2 != null)
                events.addAll(eventsChr2);
        } else if (cGT.length() > 3) {
            // MuDelIns model with insertions for with child

            List<EvolutionaryEventType> eventsChr1;
            List<EvolutionaryEventType> eventsChr2;

            if (pGT.charAt(0) == cGT.charAt(0) || (pGT.charAt(0) != cGT.charAt(0) && cGT.charAt(0) != pGT.charAt(2))) {
                eventsChr1 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[0], cAlleles[0]));
                eventsChr2 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[1], cAlleles[1]));
            } else {
                eventsChr1 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[0], cAlleles[1]));
                eventsChr2 = convertAlphabetEvolutionaryEvents(collectAlphabetEvolutionaryEvents(pAlleles[1], cAlleles[0]));
            }

            if (eventsChr1 != null)
                events.addAll(eventsChr1);

            if (eventsChr2 != null)
                events.addAll(eventsChr2);
        } else {
            // No insertion

            List<EvolutionaryEventType> eventsResult;
            try {
                eventsResult = convertAlphabetEvolutionaryEvents(
                        collectAlphabetEvolutionaryEvents(
                                String.join("", pAlleles),
                                String.join("", cAlleles)
                        )
                );
            } catch (IllegalStateException e) {
                throw new IllegalStateException(e.getMessage() + pGT + " -> " + cGT);
            }

            if (eventsResult != null)
                events.addAll(eventsResult);
        }

        if (events.size() == 0)
            return null;

        events.sort(new EvolutionaryEventType.EventsComparator());
        return events.toArray(new EvolutionaryEventType[0]);
    } // getEvolutionaryEvents

    private List<Character>[] collectAlphabetEvolutionaryEvents(final String parent, final String child) {
        List<Character>[] results = new ArrayList[2];

        // 0/1 -> 0/01
        if (parent.length() < child.length() && child.length() > 1 && child.charAt(0) != child.charAt(1)) {
            results[0] = new ArrayList<>();
            results[1] = new ArrayList<>();

            if (parent.charAt(0) == child.charAt(0)) {
                results[0].add('0');
                results[1].add('1');
            } else {
                results[0].add('1');
                results[1].add('0');
            }
        }

        List<Integer> pLeftIdx = new ArrayList<>();
        List<Integer> cMatchedIdx = new ArrayList<>();
        List<Integer> cLeftIdx = new ArrayList<>();

        for (int pIdx = 0; pIdx < parent.length(); pIdx++) {
            boolean matched = false;

            for (int cIdx = 0; cIdx < child.length(); cIdx++) {
                if (cMatchedIdx.contains(cIdx))
                    continue;

                if (parent.charAt(pIdx) == child.charAt(cIdx)) {
                    matched = true;
                    cMatchedIdx.add(cIdx);
                    break;
                }
            }

            if (!matched) {
                if (!Character.isDigit(parent.charAt(pIdx)))
                    throw new IllegalStateException(INVALID_TRANSITION);

                pLeftIdx.add(pIdx);
            }
        }

        for (int cIdx = 0; cIdx < child.length(); cIdx++) {
            if (!cMatchedIdx.contains(cIdx))
                cLeftIdx.add(cIdx);
        }

        if (pLeftIdx.size() > 0) {
            if (results[0] == null)
                results[0] = new ArrayList<>();

            for (int pIdx : pLeftIdx)
                results[0].add(parent.charAt(pIdx));
        }

        if (cLeftIdx.size() > 0) {
            if (results[1] == null)
                results[1] = new ArrayList<>();

            for (int cIdx : cLeftIdx)
                results[1].add(child.charAt(cIdx));
        }

        return results;
    } // collectAlphabetEvolutionaryEvents

    private List<EvolutionaryEventType> convertAlphabetEvolutionaryEvents(List<Character>[] events) {
        final int size = Math.max(events[0] == null ? 0 : events[0].size(), events[1] == null ? 0 : events[1].size());

        if (size == 0)
            return null;

        List<EvolutionaryEventType> results = new ArrayList<>();

        for (int idx = 0; idx < size; idx++) {
            final char pEvent = events[0] == null ? '?' : (idx < events[0].size() ? events[0].get(idx) : '?');
            final char cEvent = events[1] == null ? '?' : (idx < events[1].size() ? events[1].get(idx) : '?');

            results.add(matchEvolutionaryEvent(pEvent, cEvent));
        }

        return results;
    } // convertAlphabetEvolutionaryEvents

    /**
     * Match an evolutionary event. Only supporting '0', '1', and non-digit alleles.
     *
     * @param from from which allele?
     * @param to   to which allele?
     * @return matched evolutionary event
     */
    protected EvolutionaryEventType matchEvolutionaryEvent(final char from, final char to) {
        final char convertedFrom = from == '0' || from == '1' ? from : '?';
        final char convertedTo = to == '0' || to == '1' ? to : '?';

        if (convertedFrom == convertedTo)
            return null;

        final String query = convertedFrom + String.valueOf(convertedTo);
        switch (query) {
            case "01":
                return EvolutionaryEventType.SINGLE_MUTATION;
            case "10":
                return EvolutionaryEventType.SINGLE_BACK_MUTATION;
            case "0?":
            case "1?":
                return EvolutionaryEventType.DELETION;
            case "?0":
            case "?1":
                return EvolutionaryEventType.INSERTION;
            default:
                throw new IllegalStateException("Unexpected value: " + query);
        }
    } // matchEvolutionaryEvent

} // ScsSubstitutionModelBase
