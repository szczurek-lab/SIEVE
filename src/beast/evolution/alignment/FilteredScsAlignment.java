package beast.evolution.alignment;

import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;

import java.util.ArrayList;
import java.util.Arrays;

@Description("Alignment based on a filter operation on another alignment")
public class FilteredScsAlignment extends ScsAlignment {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    final public Input<String> filterInput = new Input<>("filter", "specifies which of the patterns in the input alignment should be selected " +
            "First site is 1." +
            "Filter specs are comma separated, either a singleton, a range [from]-[to] or iteration [from]:[to]:[step]; " +
            "1-100 defines a range, " +
            "1-100\3 or 1:100:3 defines every third in range 1-100, " +
            "1::3,2::3 removes every third site. " +
            "Default for range [1]-[last pattern], default for iterator [1]:[last pattern]:[1]", Input.Validate.REQUIRED);

    final public Input<ScsAlignment> alignmentInput = new Input<>("scsData", "alignment to be filtered", Input.Validate.REQUIRED);

    final public Input<Long> nrOfBackgroundSitesInput = new Input<>("nrOfBackgroundSites", "number of background sites", Input.Validate.REQUIRED);

    // these triples specify a range for(i=From; i <= To; i += Step) for patterns
    int[] from;
    int[] to;
    int[] step;

    /**
     * list of indices filtered from input alignment *
     */
    int[] filter;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public FilteredScsAlignment() {
        scsSequenceInput.setRule(Input.Validate.OPTIONAL);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        ScsAlignment data = alignmentInput.get();
        m_dataType = data.m_dataType;
        taxaNames = data.taxaNames;
        nrOfBackgroundSites = nrOfBackgroundSitesInput.get();
        lociExistence = data.lociExistence;
        sequences = data.sequences;

        parseFilterSpec();
        calcFilter();
        calcPatterns();
    } // initAndValidate

    private void parseFilterSpec() {
        // parse filter specification
        String filterString = filterInput.get();
        String[] filters = filterString.split(",");
        from = new int[filters.length];
        to = new int[filters.length];
        step = new int[filters.length];
        for (int i = 0; i < filters.length; i++) {
            filterString = " " + filters[i] + " ";
            if (filterString.matches(".*-.*")) {
                // range, e.g. 1-100/3
                if (filterString.indexOf('\\') >= 0) {
                    String str2 = filterString.substring(filterString.indexOf('\\') + 1);
                    step[i] = parseInt(str2, 1);
                    filterString = filterString.substring(0, filterString.indexOf('\\'));
                } else {
                    step[i] = 1;
                }
                String[] strs = filterString.split("-");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignmentInput.get().getSiteCount()) - 1;
            } else if (filterString.matches(".*:.*:.+")) {
                // iterator, e.g. 1:100:3
                String[] strs = filterString.split(":");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignmentInput.get().getSiteCount()) - 1;
                step[i] = parseInt(strs[2], 1);
            } else if (filterString.trim().matches("[0-9]*")) {
                from[i] = parseInt(filterString.trim(), 1) - 1;
                to[i] = from[i];
                step[i] = 1;
            } else {
                throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
            }
        }
    } // parseFilterSpec

    int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    } // parseInt

    private void calcFilter() {
        boolean[] isUsed = new boolean[alignmentInput.get().getLociNr()];
        for (int i = 0; i < to.length; i++) {
            for (int k = from[i]; k <= to[i]; k += step[i]) {
                isUsed[k] = true;
            }
        }

        // count
        int k = 0;
        for (boolean b : isUsed) {
            if (b) {
                k++;
            }
        }

        // set up index set
        filter = new int[k];
        k = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (isUsed[i]) {
                filter[k++] = i;
            }
        }
    } // calcFilter

    /**
     * calculate patterns from sequence data
     */
    private void calcPatterns() {
        parsedSequences = alignmentInput.get().getParsedSequences();
        int taxonNr = parsedSequences.size();
        lociNr = filter.length;
        lociNrTotal = alignmentInput.get().lociNrTotal;

        // convert parsedSequences to transposed data
        int[][][] data = new int[lociNr][taxonNr][];
        for (int i = 0; i < taxonNr; i++) {
            for (int j = 0; j < lociNr; j++) {
                data[j][i] = parsedSequences.get(i).get(filter[j]);
            }
        }

        // save loci information
        loci = new ArrayList<>();
        for (int i = 0; i < lociNr; i++) {
            loci.add(alignmentInput.get().getLociInfo(filter[i]));
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
        }

        // find patterns for the loci
        patternIndex = new int[lociNr];
        int[][] tmp = new int[taxonNr][];
        for (int i = 0; i < lociNr; i++) {
            for (int j = 0; j < taxonNr; j++) {
                tmp[j] = parsedSequences.get(j).get(filter[i]);
            }
            patternIndex[i] = Arrays.binarySearch(sitePatterns, tmp, comparator);
        }

        // report some statistics
        Log.info.println("Filter " + filterInput.get());
        Log.info.println(getTaxonCount() + " taxa");
        Log.info.println(getSiteCount() + " sites");
        Log.info.println(getPatternCount() + " patterns");
    } // calcPatterns

    /**
     * return indices of the sites that the filter uses
     */
    public int[] indices() {
        return filter.clone();
    } // indices

    /**
     * Get the coverage of {@param taxonIndex} at {@param locusIndex}.
     *
     * @param taxonIndex which taxon?
     * @param locusIndex which locus (NOT pattern!)?
     * @return coverage
     */
    @Override
    public int getSequencingCoverageOfLocus(int taxonIndex, int locusIndex) {
        return super.getSequencingCoverageOfLocus(taxonIndex, this.filter[locusIndex]);
    }

    /**
     * Get alternative nucleotides of {@param taxonIndex} at {@param locusIndex}.
     *
     * @param taxonIndex which taxon?
     * @param locusIndex which locus (NOT pattern!)?
     * @return alternative nucleotides
     */
    @Override
    public char[] getAltNucs(int taxonIndex, int locusIndex) {
        return super.getAltNucs(taxonIndex, this.filter[locusIndex]);
    }

    /**
     * For writing to VCF file.
     *
     * @param taxonIndex which taxon?
     * @param locusIndex which locus (NOT pattern)?
     * @param altNucs    locus-wise significant alt nucs
     * @return the first is read counts of ref alt, the rest are read counts of alt nucs specified in {@param altNucs}
     */
    @Override
    public int[] getAlleleDepth(int taxonIndex, int locusIndex, char[] altNucs) {
        return super.getAlleleDepth(taxonIndex, this.filter[locusIndex], altNucs);
    }

}
