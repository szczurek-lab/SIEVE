package beast.app.variantcaller;

import beast.core.Description;
import beast.math.util.MathFunctions;
import beast.util.NoFileHeaderFoundException;
import beast.util.NoMatchLengthException;
import beast.util.NotSingleException;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

@Description("Load allelic sequencing log in an independent thread for variant calling.")
public class AllelicSeqLogProcessor implements Runnable {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private int numOfThreads = 1;
    private ExecutorService pool = null;

    private final String inputFileName;

    private List<AllelicSeqLogEstimatesCaller> allelicSeqLogEstimatesCaller;

    private int numOfMatrices = 0;

    private Map<String, String> sitesMap;
    private String sitesMapStr;
    private Map<String, Map<String, double[][]>> allelicSeqEstimates;

    private boolean isSuccessful;

    static PrintStream out;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public AllelicSeqLogProcessor(int numOfThreads,
                                  final String inputFileName,
                                  final PrintStream out,
                                  EstimatesTypeCollection.EstimatesType estimatesType,
                                  EstimatesTypeCollection.ModeKDEType modeKDEType,
                                  int burninPercentage) {
        this.numOfThreads = numOfThreads;
        this.inputFileName = inputFileName;
        AllelicSeqLogProcessor.out = out;

        if (numOfThreads > 1) {
            pool = Executors.newFixedThreadPool(numOfThreads);
        }

        this.sitesMap = new HashMap<>();

        allelicSeqLogEstimatesCaller = new ArrayList<>();
        for (int i = 0; i < numOfThreads; i++) {
            allelicSeqLogEstimatesCaller.add(new AllelicSeqLogEstimatesCaller(i,
                    estimatesType, modeKDEType, burninPercentage));
        }

        this.allelicSeqEstimates = new HashMap<>();

        this.isSuccessful = false;
    }

    public AllelicSeqLogProcessor(
            final String inputFileName,
            final PrintStream out,
            Map<String, Map<String, Double>> mcmcSamplesEstimates,
            EstimatesTypeCollection estimatesTypeCollection
    ) throws NoMatchLengthException, NotSingleException, IOException {
        if (estimatesTypeCollection == null) {
            throw new NullPointerException("Empty variable.");
        }

        this.inputFileName = inputFileName;
        AllelicSeqLogProcessor.out = out;

        this.sitesMap = new HashMap<>();
        this.allelicSeqEstimates = new HashMap<>();

        getCachedEstimates(mcmcSamplesEstimates, estimatesTypeCollection);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public void getAllelicSeqInfo() throws IOException, NoFileHeaderFoundException, NoMatchLengthException {
        BufferedReader fin = new BufferedReader(new FileReader(this.inputFileName));
        String str;
        String[] parsedLine;
        String[] sites = new String[0];
        int numOfSites = 0;

        if (!fin.ready()) {
            throw new IOException(this.inputFileName + " appears to be empty.");
        }

        while (fin.ready()) {
            str = fin.readLine().trim();

            if (!str.startsWith("#")) {
                parsedLine = str.split("\t");

                if (parsedLine[0].trim().equals("Sample")) {
                    numOfSites = parsedLine.length - 1;
                    sites = new String[numOfSites];
                    System.arraycopy(parsedLine, 1, sites, 0, numOfSites);
                } else {
                    if (numOfSites == 0) {
                        throw new NoFileHeaderFoundException(this.inputFileName + " contains no headers. Cannot " +
                                "match the IDs.");
                    } else if (numOfSites != parsedLine.length - 1) {
                        throw new NoMatchLengthException(this.inputFileName + " contains a sample (" + parsedLine[0] +
                                ") which has different values (" + (parsedLine.length - 1) + ") than expected (" +
                                numOfSites + ").");
                    } else {
                        for (int i = 0; i < numOfSites; i++) {
                            String[] tmpParsed = parsedLine[i + 1].trim().split(";");

                            if (this.numOfMatrices == 0) {
                                this.numOfMatrices = tmpParsed.length;
                            } else {
                                if (this.numOfMatrices != tmpParsed.length) {
                                    throw new NoMatchLengthException("Expected number of matrices is " + this.numOfMatrices +
                                            ", but " + tmpParsed.length + " matrices are found for site " + (i + 1) +
                                            " at sample " + parsedLine[0]);
                                }
                            }

                            for (int j = 0; j < this.numOfMatrices; j++) {
                                put(i % this.numOfThreads, sites[i], this.numOfMatrices, j, tmpParsed[j].split(","));
                            }
                        }
                    }
                }
            } else {
                // comments
                if (str.startsWith("#site")) {
                    this.sitesMapStr = appendToStr(this.sitesMapStr, str);

                    str = str.substring(1);
                    parsedLine = str.split("->");
                    sitesMap.put(parsedLine[1], parsedLine[0]);
                }
            }
        }

        fin.close();
    } // getAllelicSeqInfo

    public void getCachedEstimates(Map<String, Map<String, Double>> mcmcSamplesEstimates,
                                   EstimatesTypeCollection estimatesTypeCollection) throws IOException, NullPointerException, NoMatchLengthException, NotSingleException {
        String str;
        String[] parsedLine;
        Set<String> sites = new HashSet<>();

        boolean mcmcCommentFlag = false;
        boolean mcmcHeaderFlag = false;
        boolean allelicCommentFlag = false;
        boolean allelicHeaderFlag = false;

        BufferedReader fin = new BufferedReader(new FileReader(this.inputFileName));

        if (mcmcSamplesEstimates == null)
            throw new NullPointerException();

        if (!fin.ready())
            throw new IOException(this.inputFileName + " appears empty.");

        while (fin.ready()) {
            str = fin.readLine().trim();

            if (str.startsWith("#")) {
                if (str.startsWith("#site")) {
                    // Enter site map

                    str = str.substring(1);
                    parsedLine = str.split("->");

                    if (sites.contains(parsedLine[0])) {
                        sites.remove(parsedLine[0]);
                        this.sitesMap.put(parsedLine[1], parsedLine[0]);
                    } else
                        throw new NoMatchLengthException("Missing estimates for " + parsedLine[0]);

                } else if (str.startsWith("#MCMC")) {
                    // Enter estimates for MCMC samples

                    mcmcCommentFlag = true;
                    allelicCommentFlag = false;
                } else if (str.startsWith("#Allelic")) {
                    // Enter estimates for allelic sequencing information

                    allelicCommentFlag = true;
                    mcmcCommentFlag = false;
                } else if (str.startsWith("#Type of estimates")) {
                    // Enter estimates type
                    boolean found = false;

                    parsedLine = str.split("->");
                    for (EstimatesTypeCollection.EstimatesType i : EstimatesTypeCollection.EstimatesType.values()) {
                        if (i.toString().toLowerCase().contains(parsedLine[1].toLowerCase())) {
                            if (found)
                                throw new NotSingleException("Multiple matches of estimates type. Be more specific.");
                            else {
                                found = true;

                                estimatesTypeCollection.setEstimateType(i);
                            }
                        }
                    }
                } else if (str.startsWith("#Distribution of KDE")) {
                    // Enter distribution of KDE

                    parsedLine = str.split("->");
                    if (parsedLine[1].toLowerCase().contains("null"))
                        estimatesTypeCollection.setModeKDEType(null);
                    else {
                        boolean found = false;

                        for (EstimatesTypeCollection.ModeKDEType i : EstimatesTypeCollection.ModeKDEType.values()) {
                            if (i.toString().toLowerCase().contains(parsedLine[1].toLowerCase())) {
                                if (found)
                                    throw new NotSingleException("Multiple matches of the distribution of KDE. Be more specific.");
                                else {
                                    found = true;

                                    estimatesTypeCollection.setModeKDEType(i);
                                }
                            }
                        }
                    }
                }
            } else {

                if (str.startsWith("id")) {
                    mcmcHeaderFlag = true;
                    allelicHeaderFlag = false;
                    continue;
                } else if (str.startsWith("sites")) {
                    allelicHeaderFlag = true;
                    mcmcHeaderFlag = false;
                    continue;
                }

                parsedLine = str.split("\t");

                // mcmc variables estimates
                if (mcmcCommentFlag && mcmcHeaderFlag) {
                    if (!mcmcSamplesEstimates.containsKey(parsedLine[0]))
                        mcmcSamplesEstimates.put(parsedLine[0], new HashMap<>());

                    mcmcSamplesEstimates.get(parsedLine[0])
                            .put(parsedLine[1],
                                    Double.valueOf(parsedLine[2]));
                }

                // allelic sequencing information estimates
                if (allelicCommentFlag && allelicHeaderFlag) {
                    sites.add(parsedLine[0]);

                    if (!this.allelicSeqEstimates.containsKey(parsedLine[0]))
                        this.allelicSeqEstimates.put(parsedLine[0], new HashMap<>());

                    String[] parsedMatrices = parsedLine[2].split(";");

                    if (this.numOfMatrices == 0)
                        this.numOfMatrices = parsedMatrices.length;
                    else if (this.numOfMatrices != parsedMatrices.length)
                        throw new NoMatchLengthException("Expected number of matrices is " + this.numOfMatrices +
                                ", but " + parsedMatrices.length + " matrices are found for " + parsedLine[0] +
                                " at estimates type " + parsedLine[1]);

                    double[][] estimates = new double[this.numOfMatrices][2];
                    for (int i = 0; i < this.numOfMatrices; i++) {
                        String[] pair = parsedMatrices[i].split(",");

                        estimates[i][0] = Double.parseDouble(pair[0]);
                        estimates[i][1] = Double.parseDouble(pair[1]);
                    }

                    this.allelicSeqEstimates.get(parsedLine[0]).put(parsedLine[1], estimates);
                }
            }
        }

        // Sanity check

        if (sites.size() != 0)
            throw new NoMatchLengthException("Missing map to real site information for: " + sites.toString());

        // The number of estimates type should be consistent
        int numOfEstimates = 0;

        for (String i : mcmcSamplesEstimates.keySet()) {
            int len = mcmcSamplesEstimates.get(i).size();

            if (numOfEstimates == 0)
                numOfEstimates = len;
            else if (len != numOfEstimates)
                throw new NoMatchLengthException("The number of estimates type for entry " + i + " is " + len +
                        ", different from that of other entries (" + numOfEstimates + ").");
        }

        for (String i : this.allelicSeqEstimates.keySet()) {
            int len = this.allelicSeqEstimates.get(i).size();

            if (numOfEstimates == 0)
                numOfEstimates = len;
            else if (len != numOfEstimates)
                throw new NoMatchLengthException("The number of estimates type for entry " + i + " is " + len +
                        ", different from that of other entries (" + numOfEstimates + ").");
        }

        fin.close();
    } // getCachedEstimates

    public void put(int id, String key, int numOfMatrices) {
        if (this.allelicSeqLogEstimatesCaller.get(id).getLog(key) == null) {
            this.allelicSeqLogEstimatesCaller.get(id).putLog(key, new ArrayList[numOfMatrices][2]);
        }
    } // put

    public void put(int id, String key, int numOfMatrices, int matrixIndex, String[] values) {
        put(id, key, numOfMatrices);
        this.allelicSeqLogEstimatesCaller.get(id).putLog(key, matrixIndex, values);
    } // put

    @Override
    public void run() {
        try {
            getAllelicSeqInfo();
        } catch (Exception e) {
            e.printStackTrace();
            out.println("\nError parsing allelic sequencing log: " + e.getMessage());
            return;
        }

        this.isSuccessful = true;
        out.println("\n>>> Allelic sequencing information ready. [Thread " + Thread.currentThread().getId() + "]");
    } // run

    public boolean isSuccessful() {
        return isSuccessful;
    } // isSuccessful

    private String appendToStr(String targetStr, String str) {
        if (targetStr == null || targetStr.length() == 0) {
            targetStr = "";
        }

        targetStr += str + "\n";
        return targetStr;
    } // appendToStr

    public void collectAllelicSeqEstimates() {
        try {
            if (this.numOfThreads > 1) {
                pool.invokeAll(allelicSeqLogEstimatesCaller);
            } else {
                allelicSeqLogEstimatesCaller.get(0).call();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
            System.exit(1);
        }

        for (AllelicSeqLogEstimatesCaller i : allelicSeqLogEstimatesCaller) {
            allelicSeqEstimates.putAll(i.getAllelicSeqEstimates());
        }
    } // collectAllelicSeqEstimates

    public Map<String, Map<String, double[][]>> getAllelicSeqEstimates() {
        return this.allelicSeqEstimates;
    } // getAllelicSeqEstimates

    public Map<String, String> getSitesMap() {
        return this.sitesMap;
    } // getSitesMap

    public String getSitesMapStr() {
        return this.sitesMapStr;
    } // getSitesMapStr

    public int getNumOfMatrices() {
        return this.numOfMatrices;
    } // getNumOfMatrices


    //**********************************************
    //*               Nested classes               *
    //**********************************************

    class AllelicSeqLogEstimatesCaller implements Callable<Double> {

        // 1st dimension: number of matrices
        // 2nd dimension: 0: cov; 1: var
        private Map<String, List<Double>[][]> allelicSeqLog;
        private Map<String, Map<String, double[][]>> allelicSeqEstimates;

        private final int threadNr;

        private EstimatesTypeCollection.EstimatesType estimatesType;
        private EstimatesTypeCollection.ModeKDEType modeKDEType;
        private int burninPercentage;

        public AllelicSeqLogEstimatesCaller(int threadNr,
                                            EstimatesTypeCollection.EstimatesType estimatesType,
                                            EstimatesTypeCollection.ModeKDEType modeKDEType,
                                            int burninPercentage) {
            this.threadNr = threadNr;
            this.estimatesType = estimatesType;
            this.modeKDEType = modeKDEType;
            this.burninPercentage = burninPercentage;

            allelicSeqLog = new HashMap<>();
            allelicSeqEstimates = new HashMap<>();
        }

        public void putLog(String key, List<Double>[][] value) {
            allelicSeqLog.put(key, value);
        } // putLog

        public void putLog(String key, int matrixIndex, String[] values) {
            List<Double>[][] tmp = getLog(key);

            if (tmp[matrixIndex][0] == null) {
                tmp[matrixIndex][0] = new ArrayList<>();
            }
            tmp[matrixIndex][0].add(Double.valueOf(values[0]));

            if (tmp[matrixIndex][1] == null) {
                tmp[matrixIndex][1] = new ArrayList<>();
            }
            tmp[matrixIndex][1].add(Double.valueOf(values[1]));
        } // putLog

        public List<Double>[][] getLog(String key) {
            return allelicSeqLog.get(key);
        } // getLog

        public Map<String, Map<String, double[][]>> getAllelicSeqEstimates() {
            return this.allelicSeqEstimates;
        } // collectAllelicSeqEstimates

        @Override
        public Double call() {
            AllelicSeqLogProcessor.out.println(">>> Thread " + this.threadNr + " is starting...");

            try {
                MathFunctions.getArrEstimates(
                        this.estimatesType,
                        this.burninPercentage,
                        this.allelicSeqLog,
                        this.allelicSeqEstimates
                );
            } catch (Exception e) {
                AllelicSeqLogProcessor.out.println("Something went wrong in thread " + this.threadNr);
                e.printStackTrace();
                System.exit(1);
            }

            AllelicSeqLogProcessor.out.println(">>> Thread " + this.threadNr + " done.");

            return null;
        } // call

    } // class AllelicSeqLogEstimatesCaller

}
