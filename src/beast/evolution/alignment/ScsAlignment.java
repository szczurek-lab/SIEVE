package beast.evolution.alignment;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.Map;
import beast.core.util.Log;
import beast.evolution.alignment.sequence.ScsSequence;
import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;
import beast.math.util.MathFunctions;
import beast.util.BEASTClassLoader;
import beast.util.PackageManager;

import java.io.PrintStream;
import java.util.*;

@Description("Class representing single cell sequencing alignment data")
public class ScsAlignment extends Map<String> {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    @Override
    protected Class<?> mapType() {
        return String.class;
    }

    /**
     * default data type
     */
    protected final static String FULL_SUPS_COV_DATATYPE = "full supports-coverage";

    /**
     * directory to pick up data types from *
     */
    final static String[] IMPLEMENTATION_DIR = {"beast.evolution.datatype"};

    /**
     * ascertainment bias correction
     * <p>
     * none: no correction applied
     * lewis: conditional likelihood correction method (lewis, 2001), when the number of invariant sites is unavailable
     */
    protected enum ascBiasCorrection {none, lewis, felsenstein}

    protected enum meanAscBiasCorrection {mean, geometric}

    /**
     * list of data type descriptions, obtained from DataType classes
     */
    static List<String> types = new ArrayList<>();

    static {
        findDataTypes();
    }

    @SuppressWarnings("deprecation")
    public static void findDataTypes() {
        List<String> m_sDataTypes = PackageManager.find(beast.evolution.datatype.ReadCounts.class, IMPLEMENTATION_DIR);
        for (String dataTypeName : m_sDataTypes) {
            try {
                ReadCounts dataType = (ReadCounts) BEASTClassLoader.forName(dataTypeName).newInstance();
                String description = dataType.getTypeDescription();
                if (!types.contains(description)) {
                    types.add(description);
                }
            } catch (IllegalAccessException | InstantiationException | ClassNotFoundException e) {
                throw new IllegalArgumentException(e.getMessage());
            }
        }
    }

    final public Input<ScsLociInfo> scsLociInfoInput = new Input<>("lociInfo", "loci information " +
            "of single-cell sequencing read counts data");

    final public Input<List<ScsBackgroundInfo>> scsBackgroundInfoInput = new Input<>("backgroundInfo",
            "background information of single-cell WGS or WES sequencing read counts data", new ArrayList<>(),
            Input.Validate.OPTIONAL);

    final public Input<ScsTaxonSet> scsTaxonSetInput = new Input<>("taxa", "An optional taxon-set used " +
            "only to sort the sequences into the same order as they appear in the taxon-set.", new ScsTaxonSet(),
            Input.Validate.OPTIONAL);

    final public Input<List<ScsSequence>> scsSequenceInput = new Input<>("sequence", "single-cell " +
            "sequencing read counts data", new ArrayList<>(), Input.Validate.REQUIRED);

    final public Input<String> dataTypeInput = new Input<>("dataType", "data type, one of " + types,
            FULL_SUPS_COV_DATATYPE, types.toArray(new String[0]));

    final public Input<ascBiasCorrection> ascBiasCorrectionInput = new Input<>(
            "ascertained",
            "ascertainment bias correction, e.g., conditioning the Felsenstein likelihood on excluding constant " +
                    "sites from the alignment, one of " + Arrays.toString(ascBiasCorrection.values()),
            ascBiasCorrection.none,
            ascBiasCorrection.values()
    );

    final public Input<meanAscBiasCorrection> meanAscBiasCorrectionInput = new Input<>(
            "meanAscBiasCorrection",
            "type of mean log probabilities among constant sites. Only applicable to Felsenstein's " +
                    "bias correction; that is, when \"ascertained = felsenstein\". One of " +
                    Arrays.toString(meanAscBiasCorrection.values()),
            meanAscBiasCorrection.mean,
            meanAscBiasCorrection.values()
    );

    final public Input<Long> manipulatedNrOfBackgroundSitesInput = new Input<>("manipulatedBgSitesNum",
            "the manipulated number of background sites only for acquisition bias correction, not for " +
                    "background likelihood. This number should be set only when the number of background sites is " +
                    "not enough for acquisition bias correction.");


    /**
     * list of sequences in the alignment
     */
    protected List<ScsSequence> sequences = new ArrayList<>();

    /**
     * list of parsed sequences in the alignment
     */
    protected List<List<int[]>> parsedSequences = new ArrayList<>();

    /**
     * list of taxa names defined through the sequences in the alignment
     */
    protected List<String> taxaNames = new ArrayList<>();

    /**
     * store background information names if backgroundExistence == true
     */
    protected List<String> backgroundNames = null;

    /**
     * store list of background information in the alignment if backgroundExistence == true
     * only for WGS and WES data
     * background.size() should equal m_DataType.stateCount + 1
     */
    protected List<List<long[]>> backgroundInfo = null;

    /**
     * number of background sites
     */
    protected long nrOfBackgroundSites;

    /**
     * the manipulated number of background sites
     */
    protected long manipulatedNrOfBackgroundSites;

    protected List<Long> nrOfBackgroundPoints;

    /**
     * a flag to show whether the loci information is provided or not
     */
    protected boolean lociExistence = false;

    /**
     * loci information
     */
    protected List<VariantSiteInfo> loci = null;

    /**
     * number of loci, should be unanimous across all cells
     * if loci information is provided, lociNr is initialized with it and the number of loci for all sequences
     * is compared with it
     * if not, then lociNr is initialized with the length is first sequence, then it is compared for the length
     * of all other sequences
     */
    protected int lociNr = 0;
    protected int lociNrTotal = 0;

    /**
     * e.g., if 600 loci are provided, then bitsTotal is 3.
     */
    protected int bitsTotal = 0;

    /**
     * weight over the columns of a matrix *
     */
    protected int[] patternWeight;

    /**
     * pattern state encodings *
     */
    protected int[][][] sitePatterns; // #patterns x #taxa

    /**
     * maps site nr to pattern nr *
     */
    protected int[] patternIndex;

    /**
     * data type, useful for converting String sequence to Code sequence, and back *
     */
    protected ReadCounts m_dataType;

    protected boolean returnConstSum;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public ScsAlignment() {
    }

    /**
     * constructor for testing purposes
     *
     * @param sequences
     * @param dataType
     */
    public ScsAlignment(List<ScsSequence> sequences, String dataType) {
        for (ScsSequence sequence : sequences) {
            scsSequenceInput.setValue(sequence, this);
        }
        dataTypeInput.setValue(dataType, this);
        initAndValidate();
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        if (scsSequenceInput.get().size() == 0) {
            throw new IllegalArgumentException("No single-cell read counts sequencing data defined!");
        }

        if (scsLociInfoInput.get() != null) {
            if (scsLociInfoInput.get().isLoci()) {
                lociExistence = true;
                lociNr = scsLociInfoInput.get().getLociNr();
                loci = scsLociInfoInput.get().getLociList();
            }
        } else
            throw new RuntimeException("Loci information is mandatory!");

        // initialize the data type
        initDataType();

        // initialize the sequence List
        this.sequences = scsSequenceInput.get();

        // initialize the alignment from the given list of sequences
        initializeWithSequenceList(sequences);

        // initialize background information
        initializeBackgroundInfo();

        Log.info.println(toString(false));
    } // initAndValidate

    protected void initDataType() {
        if (!types.contains(dataTypeInput.get())) {
            throw new IllegalArgumentException("Data type + '" + dataTypeInput.get() + "' cannot be found. " +
                    "Choose one of " + Arrays.toString(types.toArray(new String[0])));
        }

        List<String> dataTypes = PackageManager.find(ReadCounts.class, IMPLEMENTATION_DIR);
        for (String dataTypeName : dataTypes) {
            ReadCounts dataType;
            try {
                dataType = (ReadCounts) BEASTClassLoader.forName(dataTypeName).newInstance();
                if (dataTypeInput.get().equals(dataType.getTypeDescription())) {
                    m_dataType = dataType;
                    break;
                }
            } catch (IllegalAccessException | InstantiationException | ClassNotFoundException e) {
                throw new IllegalArgumentException(e.getMessage());
            }
        }
    } // initDataType

    private void initializeWithSequenceList(List<ScsSequence> sequences) {
        taxaNames.clear();
        parsedSequences.clear();

        try {
            for (ScsSequence seq : sequences) {
                List<int[]> tmp = seq.getReadCounts(m_dataType);

                if (lociNr == 0)
                    lociNr = tmp.size();

                if (lociNr != tmp.size())
                    throw new RuntimeException("Incompatible loci number: expecting " + lociNr + ", observing " +
                            tmp.size() + "(" + this.getClass().getName() + ")");

                parsedSequences.add(tmp);

                if (taxaNames.contains(seq.getTaxon()))
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon() + "(" +
                            this.getClass().getName() + ")");

                taxaNames.add(seq.getTaxon());

                if (parsedSequences.size() == 0)
                    throw new RuntimeException("Sequence data expected, but none found. (" +
                            this.getClass().getName() + ")");

            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        lociNrTotal = lociNr;

        calcPatterns();
    } // initializeWithSequenceList

    private void initializeBackgroundInfo() {
        nrOfBackgroundSites = 0;
        manipulatedNrOfBackgroundSites = 0;

        if (scsBackgroundInfoInput.get() != null && scsBackgroundInfoInput.get().size() > 0) {
            java.util.Map<String, List<long[]>> bgName2Info = new HashMap<>();
            nrOfBackgroundPoints = new ArrayList<>();

            backgroundNames = new ArrayList<>();
            for (ScsBackgroundInfo bg : scsBackgroundInfoInput.get()) {
                if (!backgroundNames.contains(bg.getBackgroundName())) {
                    if (!bg.sanityCheckDataType(m_dataType))
                        throw new IllegalArgumentException("Error! Incompatible background information name given " +
                                "datatype " + m_dataType.getTypeDescription() + "(" + this.getClass().getName() + ")");

                    backgroundNames.add(bg.getBackgroundName());
                    bgName2Info.put(bg.getBackgroundName(), bg.getBackgroundList());
                    nrOfBackgroundPoints.add(bg.getBackgroundInfoPoints());
                } else
                    throw new IllegalArgumentException("Duplicate background information name found. " +
                            bg.getBackgroundName() + " defined at least twice! (" + this.getClass().getName() + ")");
            }

            // sanity check
            sanityCheckBackgroundInfo();

            nrOfBackgroundSites = nrOfBackgroundPoints.get(0) / getTaxonCount();
            manipulatedNrOfBackgroundSites = nrOfBackgroundSites;

            // sort background names in the order we prefer
            sortBackgroundInfo(bgName2Info);
        }

        // if the manipulated number of background sites is manually set
        if (manipulatedNrOfBackgroundSitesInput.get() != null) {
            if (ascBiasCorrectionInput.get() != ascBiasCorrection.felsenstein)
                throw new IllegalArgumentException("Error! 'manipulatedNrOfBackgroundSitesInput' must only be " +
                        "specified when using Felsenstein's bias correction.");

            manipulatedNrOfBackgroundSites = manipulatedNrOfBackgroundSitesInput.get();

            if (manipulatedNrOfBackgroundSites <= nrOfBackgroundSites)
                throw new IllegalArgumentException("Error! The manipulated number of background sites (" +
                        manipulatedNrOfBackgroundSites + ") must be no smaller than the number of background sites (" +
                        + nrOfBackgroundSites + ") inferred from the data; otherwise it makes no sense in the " +
                        "acquisition bias correction. Please use a larger number or remove the definition.");
        }
    } // initializeBackgroundInfo

    private void sanityCheckBackgroundInfo() {
        if ((new HashSet<>(nrOfBackgroundPoints)).size() != 1)
            throw new IllegalArgumentException("Error! The number of background sites is not consistent among " +
                    "different background information string. (" + this.getClass().getName() + ")");

        if (nrOfBackgroundPoints.get(0) % getTaxonCount() != 0)
            throw new RuntimeException("The number of background sites is not an integer. Background points: " +
                    nrOfBackgroundPoints.get(0) + ", number of cells: " + getTaxonCount() + " (" +
                    this.getClass().getName() + ")");

        if ((m_dataType instanceof CovSup && backgroundNames.size() != 3) ||
                (m_dataType instanceof FullSupsCov && backgroundNames.size() != 5))
            throw new IllegalArgumentException("Error! The number of background information is incompatible with " +
                    "the given datatype.");
    } // sanityCheckBackgroundInfo

    private void sortBackgroundInfo(final java.util.Map<String, List<long[]>> bgName2Info) {
        if (this.m_dataType instanceof CovSup)
            this.backgroundNames = ScsBackgroundInfo.getOrderedNamesCovSup();

        if (this.m_dataType instanceof FullSupsCov)
            this.backgroundNames = ScsBackgroundInfo.getOrderedNamesFullSupsCov();

        // save ordered background information list
        this.backgroundInfo = new ArrayList<>();
        for (String name : this.backgroundNames) {
            this.backgroundInfo.add(bgName2Info.get(name));
        }
    } // sortBackgroundInfo

    /**
     * calculate patterns from sequence data
     */
    private void calcPatterns() {
        int taxonNr = parsedSequences.size();

        // convert parsedSequences to transposed data
        int[][][] data = new int[lociNr][taxonNr][];
        for (int i = 0; i < taxonNr; i++) {
            for (int j = 0; j < lociNr; j++) {
                data[j][i] = parsedSequences.get(i).get(j);
            }
        }

        // sort data
        SiteComparator comparator = new SiteComparator();
        Arrays.sort(data, comparator);

        // count patterns in sorted data
        int patterns = 1;
        int[] weights = new int[lociNr];
        weights[0] = 1;
        for (int i = 1; i < lociNr; i++) {
            if (comparator.compare(data[i - 1], data[i]) != 0) {
                patterns++;
                data[patterns - 1] = data[i];
            }
            weights[patterns - 1]++;
        }

        // reserve memory for patterns
        patternWeight = new int[patterns];
        sitePatterns = new int[patterns][taxonNr][];
        for (int i = 0; i < patterns; i++) {
            patternWeight[i] = weights[i];
            System.arraycopy(data[i], 0, sitePatterns[i], 0, taxonNr);
            System.arraycopy(data[i], 0, sitePatterns[i], 0, taxonNr);
        }

        // find patterns for the loci
        patternIndex = new int[lociNr];
        int[][] tmp = new int[taxonNr][];
        for (int i = 0; i < lociNr; i++) {
            for (int j = 0; j < taxonNr; j++) {
                tmp[j] = parsedSequences.get(j).get(i);
            }

            patternIndex[i] = Arrays.binarySearch(sitePatterns, tmp, comparator);
        }
    } // calcPatterns

    /**
     * Pretty printing of vital statistics of an alignment including id, #taxa, #sites, #patterns
     *
     * @param singleLine true if the string should fit on one line
     * @return string representing this alignment
     */
    public String toString(boolean singleLine) {
        StringBuilder builder = new StringBuilder();
        builder.append(getClass().getSimpleName() + "(" + getID() + ")");

        if (singleLine) {
            builder.append(": [taxa, loci, unique patterns] = [" + getTaxonCount() + ", " + getLociNr() + ", " +
                    getPatternCount() + "]");
        } else {
            builder.append("\n");
            builder.append("  " + getTaxonCount() + " taxa");
            builder.append("\n");
            builder.append("  " + getLociNr() + (getLociNr() == 1 ? " locus" : " loci"));
            builder.append("\n");
            if (getLociNr() > 1) {
                builder.append("  " + getPatternCount() + " unique patterns");
                builder.append("\n");
            }
        }
        return builder.toString();
    } // toString


    //************************************************
    //*                Nested classes                *
    //************************************************

    /**
     * ListImmuListComparator is used for ordering the loci,
     * which makes it easy to identify patterns.
     */
    static class SiteComparator implements Comparator<int[][]> {
        @Override
        public int compare(int[][] o1, int[][] o2) {
            final int len = o1[0].length;

            for (int i = 0; i < o1.length; i++) {
                for (int j = 0; j < len; j++) {
                    if (o1[i][j] > o2[i][j])
                        return 1;

                    if (o1[i][j] < o2[i][j])
                        return -1;
                }
            }
            return 0;
        }
    } // class SiteComparator


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    /*
     * assorted getters and setters *
     */
    public List<String> getTaxaNames() {
        if (taxaNames.size() == 0) {
            try {
                initAndValidate();
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }
        return taxaNames;
    } // getTaxaNames

    /**
     * Returns a List of immutable Integer Lists where each Integer List represents
     * the sequence corresponding to a taxon.  The taxon is identified by
     * the position of the immutable Integer List in the outer List, which corresponds
     * to the nodeNr of the corresponding leaf node and the position of the
     * taxon name in the taxaNames list.
     *
     * @return integer representation of sequence alignment
     */
    public List<List<int[]>> getParsedSequences() {
        return parsedSequences;
    } // getParsedSequences

    public ReadCounts getDataType() {
        return m_dataType;
    } // getReadCount

    /**
     * @return number of taxa in Alignment.
     */
    public int getTaxonCount() {
        return taxaNames.size();
    } // getTaxonCount

    public int getTaxonIndex(String id) {
        return taxaNames.indexOf(id);
    } // getTaxonIndex

    /**
     * @return Number of unique character patterns in alignment.
     */
    public int getPatternCount() {
        return sitePatterns.length;
    } // getPatternCount

    public int[][][] getPattern() {
        return sitePatterns;
    } // getPattern

    public int[][] getPattern(int patternIndex_) {
        return sitePatterns[patternIndex_];
    } // getPattern

    public int[] getPattern(int taxonIndex, int patternIndex_) {
        return sitePatterns[patternIndex_][taxonIndex];
    } // getPattern

    public int getSequencingCoverageOfPattern(int taxonIndex, int patternIndex_) {
        return m_dataType.getSequencingCoverage(sitePatterns[patternIndex_][taxonIndex]);
    } // getSequencingCoverage

    public int getSequencingCoverage(int[] counts) {
        return m_dataType.getSequencingCoverage(counts);
    } // getSequencingCoverage

    /**
     * Retrieve the "weight" of a particular pattern: the number of sites
     * having that pattern.
     *
     * @param patternIndex_ Index into pattern array.
     * @return pattern weight
     */
    public int getPatternWeight(int patternIndex_) {
        return patternWeight[patternIndex_];
    } // getPatternWeight

    /**
     * Retrieve index of pattern corresponding to a particular site.
     *
     * @param locus Index of site.
     * @return Index of pattern.
     */
    public int getPatternIndex(int locus) {
        return patternIndex[locus];
    } // getPatternIndex

    public int[] getPatternIndexArray() {
        return patternIndex;
    } // getPatternIndex

    /**
     * @return Total number of sites in alignment.
     */
    public int getSiteCount() {
        return patternIndex.length;
    }

    /**
     * @return Total number of loci in alignment.
     */
    public int getLociNr() {
        return lociNr;
    } // getLociNr

    public VariantSiteInfo getLociInfo(int locus) {
        return loci.get(locus);
    } // getLociInfo

    /**
     * Get the coverage of {@param taxonIndex} at {@param index}.
     *
     * @param taxonIndex which taxon?
     * @param locusIndex which locus (NOT pattern!)?
     * @return coverage
     */
    public int getSequencingCoverageOfLocus(int taxonIndex, int locusIndex) {
        return this.sequences.get(taxonIndex).getCoverage(locusIndex);
    } // getCoverage

    /**
     * Get alternative nucleotides of {@param taxonIndex} at {@param locusIndex}.
     *
     * @param taxonIndex which taxon?
     * @param locusIndex which locus (NOT pattern!)?
     * @return alternative nucleotides
     */
    public char[] getAltNucs(int taxonIndex, int locusIndex) {
        return this.sequences.get(taxonIndex).getAltNucs(locusIndex);
    } // getAltNucs

    /**
     * For writing to VCF file.
     *
     * @param taxonIndex which taxon?
     * @param locusIndex which locus (NOT pattern)?
     * @param altNucs    locus-wise significant alt nucs
     * @return the first is read counts of ref alt, the rest are read counts of alt nucs specified in {@param altNucs}
     */
    public int[] getAlleleDepth(int taxonIndex, int locusIndex, char[] altNucs) {
        int[] results = new int[(altNucs != null && altNucs.length > 0) ? altNucs.length + 1 : 1];

        results[0] = this.sequences.get(taxonIndex).getRefReads(locusIndex);

        if (altNucs != null && altNucs.length > 0) {
            final char[] altNucsCell = this.sequences.get(taxonIndex).getAltNucs(locusIndex);

            for (int i = 0; i < altNucs.length; i++) {
                for (int j = 0; j < altNucsCell.length; j++) {
                    if (Character.toUpperCase(altNucs[i]) == Character.toUpperCase(altNucsCell[j])) {
                        results[i + 1] = this.sequences.get(taxonIndex).getAltReads(locusIndex, j);
                        break;
                    }
                }
            }
        }

        return results;
    } // getAlleleDepth

    /**
     * Retrieve an array containing the number of times each character pattern
     * occurs in the alignment.
     *
     * @return Pattern weight array.
     */
    public int[] getWeights() {
        return patternWeight;
    } // getWeights

    public List<String> getBackgroundNames() {
        return backgroundNames;
    }

    public List<List<long[]>> getBackgroundInfo() {
        return backgroundInfo;
    }

    public long getNrOfBackgroundSites() {
        return nrOfBackgroundSites;
    }

    public long getManipulatedNrOfBackgroundSites() {
        return manipulatedNrOfBackgroundSites;
    }

    public boolean isAscBiasCorrection() {
        return !(ascBiasCorrectionInput.get().equals(ascBiasCorrection.none));
    } // isAscBiasCorrection

    /**
     * When branch rate model is not strict molecular clock model (Relaxed Molecular Clock),
     * turning on Felsenstein's bias correction will force the entire tree to be dirty to avoid numeric issues.
     *
     * @return whether using Felsenstein's bias correction
     */
    public boolean isForcingTreeDirtyRMC() {
        return ascBiasCorrectionInput.get().equals(ascBiasCorrection.felsenstein);
    } // isForcingTreeDirtyRMC

    /**
     * get ascertainment bias correction value
     *
     * @param logConstRoot   log-likelihood for constant site
     * @param returnConstSum return only the sum of the likelihoods of constant sites
     * @return ascertainment bias correction value
     */
    public double getAscBiasCorrection(
            final double[] logConstRoot,
            boolean returnConstSum
    ) {
        if (isAscBiasCorrection()) {
            switch (ascBiasCorrectionInput.get()) {
                case lewis:
                    this.returnConstSum = false;
                    return getAscBiasCorrectionLewis(logConstRoot);
                case felsenstein:
                    this.returnConstSum = returnConstSum;
                    return getAscBiasCorrectionFelsenstein(logConstRoot, returnConstSum);
                default:
                    return 0.0;
            }
        } else
            return 0.0;
    } // getAscBiasCorrection

    /**
     * conditional likelihood correction method (Paul Lewis, 2001)
     * <p>
     * -n * log(1.0 - c)
     * n: number of candidate mutated sites
     *
     * @param logConstRoot log-likelihood for constant site
     * @return ascertainment bias correction value
     */
    protected double getAscBiasCorrectionLewis(double[] logConstRoot) {
        double corr = 0.0;
        for (int i = 0; i < getPatternCount(); i++)
            corr -= Math.log(1 - Math.exp(logConstRoot[i])) * getPatternWeight(i);

        return corr;
    }

    /**
     * conditional likelihood correction method (Adam D. Leaché, 2015)
     * <p>
     * w * log(c)
     * w: number of background sites
     *
     * @param logConstRoot   log-likelihood for constant site
     * @param returnConstSum return only the sum of the likelihoods of constant sites
     * @return ascertainment bias correction value
     */
    protected double getAscBiasCorrectionFelsenstein(double[] logConstRoot, boolean returnConstSum) {
        double corr;
        double corrected;

        switch (meanAscBiasCorrectionInput.get()) {
            case mean:
                corr = MathFunctions.logSumExp(logConstRoot);
                corrected = (corr - Math.log(logConstRoot.length)) * getManipulatedNrOfBackgroundSites();
                break;
            case geometric:
                corr = MathFunctions.sum(logConstRoot);
                corrected = corr * getManipulatedNrOfBackgroundSites() / logConstRoot.length;
                break;
            default:
                throw new IllegalArgumentException("Error! Illegal value for meanAscBiasCorrection.");
        }

        return returnConstSum ? corr : corrected;
    }

    /**
     * get ascertainment bias correction value
     *
     * @param logConstRoot bias correction terms
     * @return ascertainment bias correction value
     */
    public double getAscBiasCorrection(final double[] logConstRoot) {
        if (isAscBiasCorrection()) {
            switch (ascBiasCorrectionInput.get()) {
                case felsenstein:
                    return getAscBiasCorrectionFelsenstein(logConstRoot);
                default:
                    return 0.0;
            }
        } else
            return 0.0;
    } // getAscBiasCorrection

    /**
     * conditional likelihood correction method (Adam D. Leaché, 2015)
     * <p>
     * w * log(c)
     * w: number of background sites
     *
     * @param logConstRoot bias correction terms
     * @return ascertainment bias correction value
     */
    protected double getAscBiasCorrectionFelsenstein(final double[] logConstRoot) {
        switch (meanAscBiasCorrectionInput.get()) {
            case mean:
                return (MathFunctions.logSumExp(logConstRoot) - Math.log(lociNr)) * getManipulatedNrOfBackgroundSites();
            case geometric:
                return MathFunctions.sum(logConstRoot) * getManipulatedNrOfBackgroundSites() / logConstRoot.length;
            default:
                throw new IllegalArgumentException("Error! Illegal value for meanAscBiasCorrection.");
        }
    }

    public boolean isReturnConstSum() {
        return returnConstSum;
    }


    //*******************************************
    //*                 Logger                  *
    //*******************************************

    /**
     * write site header
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    public void logSiteHeader(int start, PrintStream out) {
        for (int i = 0; i != lociNr; i++) {
            out.print("site" + getSiteIndexString(i + start));

            if (i < lociNr - 1)
                out.print("\t");
        }
    } // logSiteHeader

    /**
     * write loci information to out in a form of:
     * site code->chromosome,site number,reference nucleotide,alternative nucleotides
     *
     * @param start the start number of site code
     * @param out   apparently
     */
    public void logSiteMap(int start, PrintStream out) {
        for (int i = 0; i != lociNr; i++) {
            out.print("#site" + getSiteIndexString(i + start));
            out.print("->");
            out.println(loci.get(i).toString(",", true));
        }
    } // logSiteInfo

    public void setBitsTotal() {
        this.bitsTotal = getBits(this.lociNrTotal);
    } // setBitsTotal

    public int getBits(int num) {
        int i, divider = 10, result;
        if (num >= 10)
            i = 1;
        else if (num >= 0)
            i = 0;
        else
            return 0;

        do {
            i++;
            result = num / divider;
            divider *= 10;
        } while (result >= 10);

        return i;
    } // getBits

    public String getSiteIndexString(int index) {
        if (this.bitsTotal == 0)
            setBitsTotal();

        int bits = getBits(index);
        String result = String.join("", Collections.nCopies(this.bitsTotal - bits, "0"));
        result += index;

        return result;
    } // getSiteIndexString

} // class ScsAlignment
