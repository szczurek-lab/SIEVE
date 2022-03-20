package beast.app.variantcaller;

import beast.evolution.alignment.VariantSiteInfo;
import org.jetbrains.annotations.NotNull;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.*;

public class GenotypeAdoStateProcessor {

    private final String inputFileName;
    private final double burnin;

    private final int threadNum;
    private final ExecutorService threadPool;
    private GenotypeAdoState[] genotypeAdoStates = null;

    private final String genotypeFileName;
    private final String adoStateFileName;

    List<SampleProcessor> sampleProcessors;

    private static PrintStream out;

    public GenotypeAdoStateProcessor(
            int threadNum,
            String inputFileName,
            String genotypeFileName,
            String adoStateFileName,
            int burninPercentage,
            final PrintStream out
    ) {
        this.threadNum = threadNum;
        this.threadPool = Executors.newFixedThreadPool(threadNum);

        this.inputFileName = inputFileName;
        this.genotypeFileName = genotypeFileName;
        this.adoStateFileName = adoStateFileName;
        this.burnin = burninPercentage / 100.0;

        this.sampleProcessors = new ArrayList<>();

        GenotypeAdoStateProcessor.out = out;
        SampleProcessor.setLogStream(out);
    }

    public void processGenotypeAdoStateSamples() throws IOException, InterruptedException, ExecutionException {
        int cellNum = 0, siteNum = 0, sampleIndex = 0;
        String[] cellNames, pseudoSiteNames = null;
        Map<String, VariantSiteInfo> siteNamesMap = new HashMap<>();
        String line;

        Future<GenotypeAdoState[]>[] threadedGenotypeAdoStates = new Future[this.threadNum];

        /*
         * 0: processing the first line which contains pseudo site names.
         * 1: reaching samples but the number of cells is unknown.
         * 2: processing samples.
         * 3: processing real cell names.
         * 4: processing real site names.
         */
        int mode = 0;

        BufferedReader fin = new BufferedReader(new FileReader(this.inputFileName));
        if (!fin.ready())
            throw new IOException(this.inputFileName + " appears empty.");
        while (fin.ready()) {
            line = fin.readLine().trim();

            if (line.isEmpty()) continue;

            if (line.startsWith("#Cells:"))
                mode = 3;
            else if (line.startsWith("#Sites map")) {
                mode = 4;
                continue;
            } else if (mode != 4 && line.startsWith("#"))
                continue;

            switch (mode) {
                case 0:
                    pseudoSiteNames = getPseudoSiteNames(line);
                    siteNum = pseudoSiteNames.length;
                    mode++;
                    break;
                case 1:
                    cellNum = getCellNum(line);
                    genotypeAdoStates = new GenotypeAdoState[cellNum * siteNum];
                    for (int i = 0; i < this.threadNum; i++)
//                        sampleProcessors.add(
//                                new SampleProcessor(
//                                        i,
//                                        cellNum,
//                                        siteNum,
//                                        pseudoSiteNames
//                                )
//                        );
                    mode++;
                case 2:
                    sampleProcessors.get(sampleIndex % this.threadNum).addLine(line);
                    sampleIndex++;
                    break;
                case 3:
                    for (int i = 0; i < this.threadNum; i++)
                        threadedGenotypeAdoStates[i] = threadPool.submit(sampleProcessors.get(i));
                    cellNames = getCellNames(line);
                    if (cellNum != cellNames.length)
                        throw new RuntimeException("Error! The number of cells in samples does not match the number of real cell names.");
                    break;
                case 4:
                    processSiteMap(line, siteNamesMap);
                    break;
                default:
                    throw new RuntimeException("Error! Unrecognized mode. The file containing samples of genotypes and ado states is illegally formed.");
            }
        }
        fin.close();

        GenotypeAdoState.setValidSampleIndex((int) Math.ceil(this.burnin * sampleIndex));

        for (int i = 0; i < siteNum; i++) {
            for (int j = 0; j < cellNum; j++) {
                final int index = i * cellNum + j;

                for (int k = 0; k < threadNum; k++) {
                    if (k == 0)
                        genotypeAdoStates[index] = threadedGenotypeAdoStates[k].get()[index];
                    else
                        genotypeAdoStates[index].addSamples(threadedGenotypeAdoStates[k].get()[index].getSamples());
                }

                assert genotypeAdoStates[index].checkSampleSize(sampleIndex);
                genotypeAdoStates[index].setSiteName(siteNamesMap.get(genotypeAdoStates[index].getPseudoSiteName()));
                genotypeAdoStates[index].sortAndCount(true);
            }
        }
    }

    private String[] getPseudoSiteNames(final String line) {
        String[] parsedLine = line.split("\t");
        return Arrays.copyOfRange(parsedLine, 1, parsedLine.length);
    }

    private int getCellNum(final String line) {
        return line.split("\t")[1].split(";").length;
    }

    private String[] getCellNames(final String line) {
        return line.replaceFirst("^#Cells:", "").trim().replaceAll("\\s+", "").split(",");
    }

    private void processSiteMap(final String line, @NotNull Map<String, VariantSiteInfo> siteNamesMap) {
        final String[] parsedLine = line.replaceFirst("#", "").trim().replaceAll("\\s+", "").split("->");
        final String[] parsedSite = parsedLine[1].split(",");

        assert parsedSite.length == 4;

        siteNamesMap.put(
                parsedLine[0],
                new VariantSiteInfo(
                        parsedSite[0],
                        Long.parseLong(parsedSite[1]),
                        parsedSite[2],
                        parsedSite[3].split("")
                )
                );
    }

}

class SampleProcessor implements Callable<GenotypeAdoState[]> {

    private final int id;
    private final List<String> lines;
    private final int cellNum;
    private final int siteNum;
    private static PrintStream out;
    private final String[] pseudoSiteNames;

    public SampleProcessor(
            int id,
            int cellNum,
            int fromSite,
            int toSite,
            String[] pseudoSiteNames
    ) {
        this.id = id;
        this.lines = new ArrayList<>();
        this.cellNum = cellNum;
        this.siteNum = toSite - fromSite;
        this.pseudoSiteNames = new String[this.siteNum];
        System.arraycopy(pseudoSiteNames, fromSite, this.pseudoSiteNames, 0, this.siteNum);
    }

    public static void setLogStream(final PrintStream s) {
        out = s;
    }

    public void addLine(final String line) {
        this.lines.add(line);
    }

    @Override
    public GenotypeAdoState[] call() {
        synchronized (SampleProcessor.class) {
            out.println("[Thread " + this.id + "] Starting to process samples of genotypes and ado states...");
        }

        GenotypeAdoState[] genotypeAdoStates = new GenotypeAdoState[this.siteNum * this.cellNum];

        for (String line : lines) {
            final String[] parsedLine = line.trim().split("\t");
            final long sample = Long.parseLong(parsedLine[0]);

            assert parsedLine.length == siteNum + 1;

            for (int i = 0; i < siteNum; i++) {
                final String[] parsedSite = parsedLine[i + 1].trim().split(";");

                assert parsedSite.length == cellNum;

                for (int j = 0; j < cellNum; j++) {
                    final int index = i * cellNum + j;

                    if (genotypeAdoStates[index] == null)
                        genotypeAdoStates[index] = new GenotypeAdoState(pseudoSiteNames[i]);

                    final String[] parsedCell = parsedSite[j].trim().split(",");
                    genotypeAdoStates[index].addSample(
                            new Sample(
                                    sample,
                                    parsedCell[0],
                                    Integer.parseInt(parsedCell[1])
                            )
                    );
                }
            }
        }

        synchronized (SampleProcessor.class) {
            out.println("[Thread " + this.id + "] Done.");
        }

        return genotypeAdoStates;
    }

}

class GenotypeAdoState {

    static private int validSampleIndex = 0;

    private final String pseudoSiteName;
    private VariantSiteInfo siteName;

    private final List<Sample> samples;

    private final Map<String, Integer> genotypeCount;
    private final Map<Integer, Integer> adoStateCount;

    private final Map<String, Double> genotypeFreq;
    private final Map<Integer, Double> adoStateFreq;

    public static void setValidSampleIndex(int val) {
        validSampleIndex = val;
    }

    public GenotypeAdoState(String pseudoSiteName) {
        this.pseudoSiteName = pseudoSiteName;

        this.samples = new ArrayList<>();

        this.genotypeCount = new HashMap<>();
        this.adoStateCount = new HashMap<>();

        this.genotypeFreq = new HashMap<>();
        this.adoStateFreq = new HashMap<>();
    }

    public boolean checkSampleSize(int val) {
        return val == samples.size();
    }

    public void setSiteName(VariantSiteInfo siteName) {
        this.siteName = siteName;
    }

    public String getPseudoSiteName() {
        return this.pseudoSiteName;
    }

    public void addSample(Sample val) {
        if (this.samples.contains(val))
            throw new RuntimeException("Error! Duplicate samples found for sample " + val.getSample() + " and site " + pseudoSiteName);

        this.samples.add(val);
    }

    public void addSamples(List<Sample> vals) {
        this.samples.addAll(vals);
    }

    public List<Sample> getSamples() {
        return this.samples;
    }

    public void sortAndCount(boolean doSort) {
        if (doSort)
            Collections.sort(this.samples);

        for (int i = validSampleIndex; i < samples.size(); i++) {
            final String genotype = samples.get(i).getGenotype();
            final int adoState = samples.get(i).getAdoState();

            genotypeCount.put(genotype, genotypeCount.getOrDefault(genotype, 0) + 1);
            adoStateCount.put(adoState, adoStateCount.getOrDefault(adoState, 0) + 1);
        }

        getFreqFromCount(genotypeCount, genotypeFreq);
        getFreqFromCount(adoStateCount, adoStateFreq);
    }

    private static <K> void getFreqFromCount(
            @NotNull Map<K, Integer> count,
            @NotNull Map<K, Double> freq
    ) {
        final int valSum = Arrays.stream(count.values().toArray(new Integer[0])).mapToInt(Integer::intValue).sum();

        for (K key : count.keySet())
            freq.put(key, ((double) count.get(key)) / valSum);
    }

}

class Sample implements Comparable<Sample> {

    private final long sample;
    private final String genotype;
    private final int adoState;

    public Sample() {
        throw new RuntimeException("Error! Unsupported constructor.");
    }

    public Sample(
            long sample,
            String genotype,
            int adoState
    ) {
        this.sample = sample;
        this.genotype = genotype;
        this.adoState = adoState;
    }

    public long getSample() {
        return sample;
    }

    public String getGenotype() {
        return genotype;
    }

    public int getAdoState() {
        return adoState;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Sample sample1 = (Sample) o;

        if (sample != sample1.sample) return false;
        if (adoState != sample1.adoState) return false;
        return genotype.equals(sample1.genotype);
    }

    @Override
    public int hashCode() {
        int result = (int) (sample ^ (sample >>> 32));
        result = 31 * result + genotype.hashCode();
        result = 31 * result + adoState;
        return result;
    }

    @Override
    public int compareTo(@NotNull Sample o) {
        if (this == o) return 0;

        if (this.sample > o.sample)
            return 1;
        else if (this.sample < o.sample)
            return -1;

        final int val = this.genotype.compareTo(o.genotype);
        if (val != 0) return val;

        if (this.adoState > o.adoState)
            return 1;
        else if (this.adoState < o.adoState)
            return -1;
        return 0;
    }

}
