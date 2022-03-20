package beast.app.variantcaller;

import beast.app.BeastMCMC;
import beast.app.tools.LogCombiner;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.app.variantcallingmodeladaptor.VariantCallingModelAdaptor;
import beast.core.Runnable;
import beast.core.VariantCall;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.math.util.MathFunctions;
import beast.util.*;
import jam.console.ConsoleApplication;
import org.jetbrains.annotations.NotNull;
import org.json.JSONException;
import org.xml.sax.SAXException;

import javax.swing.*;
import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class VariantCaller {

    final static private String FLOAT_FORMAT = "%.10f";


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public static PrintStream progressStream = Log.err;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public VariantCaller(
            int numOfThreads,
            int burninPercentage,
            EstimatesTypeCollection.EstimatesType estimatesType,
            EstimatesTypeCollection.ModeKDEType modeKDEType,
            String mcmcSamplesFileName,
            String allelicInfoFileName,
            String gtAdoSamplesFileName,
            String modelConfigFileName,
            String inputTreeFileName,
            String outputVCFileName,
            boolean saveDetails,
            final int cellThreshold,
            final boolean useMeanRate,
            final boolean useOnlyBranchLength
    ) throws IOException, InterruptedException {
        Map<String, List<Double>> mcmcSamples = new HashMap<>();
        Map<String, Map<String, Double>> mcmcSamplesEstimates = new HashMap<>();

        AllelicSeqLogProcessor allelicSeqLog = null;
        if (allelicInfoFileName != null)
            allelicSeqLog = new AllelicSeqLogProcessor(
                    1,
                    allelicInfoFileName,
                    Log.err,
                    estimatesType,
                    modeKDEType,
                    burninPercentage
            );

        GenotypeAdoStateProcessor gtAdoLog = null;
        String genotypeFileName, adoStateFileName;
        if (gtAdoSamplesFileName != null) {
            genotypeFileName = getVariantCallingFileName(gtAdoSamplesFileName, outputVCFileName, "genotype");
            adoStateFileName = getVariantCallingFileName(gtAdoSamplesFileName, outputVCFileName, "ado");

            gtAdoLog = new GenotypeAdoStateProcessor(
                    numOfThreads,
                    gtAdoSamplesFileName,
                    genotypeFileName,
                    adoStateFileName,
                    burninPercentage,
                    Log.err
            );
        }

        VariantCallingModelAdaptor variantCallingModelAdaptor;

        int numOfModels;
        String[] modelFileNames;
        EstimatesTypeCollection[] modelEstimateTypes;
        Runnable[] models;
        double[] logPosteriors;
        double[] logTreeLikelihoods;

        final long startTime = System.currentTimeMillis();

        // 1. Get allelic sequencing log and estimates (possibly in an independent thread)
        progressStream.println();
        progressStream.print(">>> Reading the allelic sequencing log...");
        if (allelicSeqLog != null) {
            try {
                allelicSeqLog.getAllelicSeqInfo();
            } catch (Exception e) {
                e.printStackTrace();
                Log.err.println("Error parsing allelic sequencing log: " + e.getMessage());
            }

            progressStream.println(">>> Allelic sequencing information ready. [Thread " +
                    Thread.currentThread().getId() + "]");
        }

        // 2. Get MCMC samples
        progressStream.println();
        progressStream.println(">>> Reading the MCMC samples log... [Thread " + Thread.currentThread().getId() + "]");
        try {
            getMCMCSamples(mcmcSamplesFileName, mcmcSamples);
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error parsing MCMC samples log: " + e.getMessage());
            return;
        }
        progressStream.println(">>> MCMC samples ready. [Thread " + Thread.currentThread().getId() + "]");

        // 3. Get the estimates of MCMC samples
        progressStream.println();
        progressStream.println(">>> Estimating model variables... [Thread " + Thread.currentThread().getId() + "]");
        try {
            MathFunctions.getEstimates(estimatesType, modeKDEType, burninPercentage, mcmcSamples, mcmcSamplesEstimates);
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error getting estimates for allelic sequencing information: " + e.getMessage());
            return;
        }
        progressStream.println(">>> Estimates of model variables ready. [Thread " + Thread.currentThread().getId() + "]");

        // 4. Get the estimates of allelic sequencing log
        if (allelicSeqLog != null) {
            progressStream.println();
            progressStream.println(">>> Estimating allelic sequencing coverage and raw variance...");
            allelicSeqLog.collectAllelicSeqEstimates();
            progressStream.println(">>> Estimates of allelic sequencing information ready.");
        }

        // 5. Save the computed estimates
        progressStream.println();
        progressStream.print(">>> Saving estimates of mcmc samples and allelic sequencing information to ");
        String estimatesOutputFileName = getVariantCallingFileName(mcmcSamplesFileName, outputVCFileName, "estimates");
        progressStream.println("'" + estimatesOutputFileName + "'");
        saveEstimatesTo(
                estimatesOutputFileName,
                estimatesType,
                modeKDEType,
                mcmcSamplesEstimates,
                allelicSeqLog
        );
        progressStream.println(">>> Done.");

        // 6. Process samples of genotype and ado states collected during MCMC.
        if (gtAdoLog != null) {
            progressStream.println();
            progressStream.println(">>> Processing samples of genotype and ado states collected during MCMC...");
            try {
                gtAdoLog.processGenotypeAdoStateSamples();
            } catch (ExecutionException e) {
                e.printStackTrace();
                Log.err.println("Error! " + e.getMessage());
                return;
            }
            progressStream.println(">>> Done.");
        }

        // 7. Modify the provided model configuration file used in MCMC for variant calling
        progressStream.println();
        progressStream.println(">>> Modifying template configuration file for variant calling...");
        variantCallingModelAdaptor = new VariantCallingModelAdaptor(
                modelConfigFileName,
                inputTreeFileName,
                outputVCFileName,
                new EstimatesTypeCollection(estimatesType, modeKDEType),
                mcmcSamplesEstimates,
                allelicSeqLog,
                saveDetails,
                cellThreshold,
                useMeanRate,
                useOnlyBranchLength
        );
        progressStream.println(">>> Done.");
        progressStream.flush();

        // 8. Call variants
        BeastMCMC.m_nThreads = numOfThreads;
        numOfModels = variantCallingModelAdaptor.getNumOfAdaptedConfigFiles();
        modelFileNames = variantCallingModelAdaptor.getAdaptedConfigFileNames().toArray(new String[0]);
        modelEstimateTypes = variantCallingModelAdaptor.getAdaptedConfigEstimateTypes().toArray(new EstimatesTypeCollection[0]);
        models = new Runnable[numOfModels];
        logPosteriors = new double[numOfModels];
        logTreeLikelihoods = new double[numOfModels];
        callVariants(
                numOfModels,
                modelFileNames,
                modelEstimateTypes,
                models,
                logPosteriors,
                logTreeLikelihoods,
                progressStream,
                outputVCFileName
        );
        progressStream.println();
        progressStream.println(">>> Finished. [Thread " + Thread.currentThread().getId() + "]");

        final long endTime = System.currentTimeMillis();
        progressStream.println();
        progressStream.println(">>> Variant calling completed in " + (endTime - startTime) / 1000 + "s.");

    }

    public VariantCaller(
            int numOfThreads,
            String cachedEstimatesFileName,
            String modelConfigFileName,
            String inputTreeFileName,
            String outputVCFileName,
            boolean saveDetails,
            final int cellThreshold,
            final boolean useMeanRate,
            final boolean useOnlyBranchLength
    ) {
        Map<String, Map<String, Double>> mcmcSamplesEstimates = new HashMap<>();

        AllelicSeqLogProcessor allelicSeqLog = null;

        EstimatesTypeCollection estimatesTypeCollection = new EstimatesTypeCollection();

        VariantCallingModelAdaptor variantCallingModelAdaptor;

        int numOfModels;
        String[] modelFileNames;
        EstimatesTypeCollection[] modelEstimateTypes;
        Runnable[] models;
        double[] logPosteriors;
        double[] logTreeLikelihoods;

        final long startTime = System.currentTimeMillis();

        // 1. Get cached estimates
        try {
            allelicSeqLog = new AllelicSeqLogProcessor(
                    cachedEstimatesFileName,
                    Log.err,
                    mcmcSamplesEstimates,
                    estimatesTypeCollection
            );
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error parsing cached estimates file: " + e.getMessage());
            return;
        }

        // 2. Modify the provided model configuration file used in MCMC for variant calling
        progressStream.println();
        progressStream.println(">>> Modifying template configuration file for variant calling...");
        variantCallingModelAdaptor = new VariantCallingModelAdaptor(
                modelConfigFileName,
                inputTreeFileName,
                null,
                estimatesTypeCollection,
                mcmcSamplesEstimates,
                allelicSeqLog,
                saveDetails,
                cellThreshold,
                useMeanRate,
                useOnlyBranchLength
        );
        progressStream.println(">>> Done.");

        // 3. Call variants
        BeastMCMC.m_nThreads = numOfThreads;
        numOfModels = variantCallingModelAdaptor.getNumOfAdaptedConfigFiles();
        modelFileNames = variantCallingModelAdaptor.getAdaptedConfigFileNames().toArray(new String[0]);
        modelEstimateTypes = variantCallingModelAdaptor.getAdaptedConfigEstimateTypes().toArray(new EstimatesTypeCollection[0]);
        models = new Runnable[numOfModels];
        logPosteriors = new double[numOfModels];
        logTreeLikelihoods = new double[numOfModels];
        callVariants(
                numOfModels,
                modelFileNames,
                modelEstimateTypes,
                models,
                logPosteriors,
                logTreeLikelihoods,
                progressStream,
                outputVCFileName
        );
        progressStream.println();
        progressStream.println(">>> Finished.");

        final long endTime = System.currentTimeMillis();
        progressStream.println();
        progressStream.println(">>> Variant calling completed in " + (endTime - startTime) / 1000 + "s.");
    }

    public VariantCaller(
            int numOfThreads,
            @NotNull String runConfigFileName
    ) {
        Runnable model;

        progressStream.println();
        progressStream.println(">>> Loading the model...");
        try {
            model = getModel(runConfigFileName);
        } catch (JSONException | XMLParserException | IOException e) {
            e.printStackTrace();
            progressStream.println("Error parsing model configuration file: " + e.getMessage());
            return;
        }
        progressStream.println(">>> Done.");

        progressStream.println();
        progressStream.println(">>> Run the model...");
        try {
            model.run();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }
        progressStream.println(">>> Log posterior: " + String.format(FLOAT_FORMAT, ((VariantCall) model).getLogPosterior()));
        progressStream.println(">>> Log tree likelihood: " + String.format(FLOAT_FORMAT, ((VariantCall) model).getLogTreeLikelihood()));

        progressStream.println();
        progressStream.println(">>> Calling variants... ");
        ((VariantCall) model).callVariantsAndGenerateVCF();
        progressStream.println(">>> Done.");
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public Tree getTree(final String inputTreeFileName) throws IOException, NotSingleException {
        List<Tree> parsedTrees;
        boolean isNexus = true;
        BufferedReader fin;
        String str;

        // Nexus tree or Newick tree?
        fin = new BufferedReader(new FileReader(new File(inputTreeFileName)));
        if (!fin.ready()) {
            throw new IOException(inputTreeFileName + " appears empty.");
        }
        str = fin.readLine();
        if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
            // Newick tree found
            isNexus = false;
        }

        // Parse tree
        if (isNexus) {
            NexusParser nexusParser = new NexusParser();
            nexusParser.parseFile(new File(inputTreeFileName));
            parsedTrees = nexusParser.trees;
        } else {
            parsedTrees = new ArrayList<>();
            while (fin.ready()) {
                str = fin.readLine().trim();

                Tree thisTree;
                try {
                    thisTree = new TreeParser(null, str, 0, false);
                } catch (ArrayIndexOutOfBoundsException e) {
                    thisTree = new TreeParser(null, str, 1, false);
                }

                parsedTrees.add(thisTree);
            }
        }
        fin.close();

        // Verify; only one input tree is allowed
        if (parsedTrees.size() != 1) {
            throw new NotSingleException("Only one input tree is expected, but " + parsedTrees.size() + " provided.");
        }

        return parsedTrees.get(0);
    } // getTree

    public void getMCMCSamples(final String mcmcSamplesFileName, Map<String, List<Double>> mcmcSamples) throws IOException, NoFileHeaderFoundException, NoMatchLengthException {
        BufferedReader fin = new BufferedReader(new FileReader(mcmcSamplesFileName));
        String str;
        String[] parsedLine;
        String[] variables = new String[0];
        int numOfVariables = 0;

        if (!fin.ready()) {
            throw new IOException(mcmcSamplesFileName + " appears empty.");
        }

        while (fin.ready()) {
            str = fin.readLine().trim();

            if (!str.startsWith("#")) {
                parsedLine = str.split("\t");

                if (str.startsWith("Sample")) {
                    numOfVariables = parsedLine.length - 1;

                    variables = new String[numOfVariables];
                    System.arraycopy(parsedLine, 1, variables, 0, numOfVariables);

                    for (int i = 1; i != numOfVariables + 1; i++) {
                        mcmcSamples.put(parsedLine[i], new ArrayList<>());
                    }
                } else {
                    if (numOfVariables == 0) {
                        throw new NoFileHeaderFoundException(mcmcSamplesFileName + " contains no headers. Cannot " +
                                "match the IDs.");
                    } else if (numOfVariables != parsedLine.length - 1) {
                        throw new NoMatchLengthException(mcmcSamplesFileName + " contains a sample (" + parsedLine[0] +
                                ") which has more values (" + (parsedLine.length - 1) + ") than expected (" +
                                numOfVariables + ").");
                    } else {
                        for (int i = 0; i != numOfVariables; i++) {
                            mcmcSamples.get(variables[i]).add(Double.valueOf(parsedLine[i + 1]));
                        }
                    }
                }
            }
        }

        fin.close();
    } // getMCMCSamples

    public String getVariantCallingFileName(final String baseFileName, final String rootFileName, final String suffix) {
        String results;
        results = FileNameProcessor.getRootPath(rootFileName) +
                FileNameProcessor.getBaseName(baseFileName) +
                "." + suffix;

        return results;
    } // getVariantCallingFileName

    public void saveEstimatesTo(
            String outputFileName,
            final EstimatesTypeCollection.EstimatesType estimatesType,
            final EstimatesTypeCollection.ModeKDEType modeKDEType,
            final Map<String, Map<String, Double>> mcmcSamplesEstimates,
            final AllelicSeqLogProcessor allelicSeqInfo
    ) throws NullPointerException, FileNotFoundException {
        if (outputFileName == null)
            throw new NullPointerException("No output file name provided.");

        if (System.getProperty("variant.calling.file.prefix") != null)
            outputFileName = System.getProperty("variant.calling.file.prefix") + outputFileName;

        File output = new File(outputFileName);
        if (outputFileName.contains("/"))
            output.getParentFile().mkdirs();

        final PrintStream stream = new PrintStream(output);

        // Print header
        stream.println("#Estimates for MCMC samples and allelic sequencing information.");
        stream.println("#Format: variables(id, for mcmc)/sites(for allelic sequencing information)\t" +
                "estimate type(mean, median, mode, or all)\t" +
                "estimates('allelic sequencing coverage,raw variance' for allelic sequencing information, which is " +
                "also semicolon separated for site-wise categories)");
        stream.println("#");

        stream.print("#Type of estimates->");
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MEDIAN || estimatesType == EstimatesTypeCollection.EstimatesType.MEAN)
            stream.println(estimatesType.toString().toLowerCase());
        else if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE)
            stream.println("mode");
        else
            stream.println("all");

        stream.print("#Distribution of KDE->");
        if (estimatesType == EstimatesTypeCollection.EstimatesType.MEDIAN ||
                estimatesType == EstimatesTypeCollection.EstimatesType.MEAN ||
                estimatesType == EstimatesTypeCollection.EstimatesType.ALL)
            stream.println("null");
        else
            stream.println(modeKDEType.toString().toLowerCase());
        stream.println("#");

        // Print estimates for MCMC samples
        stream.println("#MCMC samples");
        stream.println("id\testimate type\tvalue");
        for (String var : mcmcSamplesEstimates.keySet())
            for (String type : mcmcSamplesEstimates.get(var).keySet())
                stream.println(var + "\t" + type.toLowerCase() + "\t" + mcmcSamplesEstimates.get(var).get(type));
        stream.println("#");

        // Print estimates for allelic sequencing information
        if (allelicSeqInfo != null) {
            Map<String, Map<String, double[][]>> allelicSeqEstimates = allelicSeqInfo.getAllelicSeqEstimates();

            stream.println("#Allelic sequencing information");
            stream.println("sites\testimate type\tallelic sequencing coverage\traw variance");
            List<String> keys = new ArrayList<>(allelicSeqEstimates.keySet());
            Collections.sort(keys);
            for (String var : keys) {
                for (String type : allelicSeqEstimates.get(var).keySet()) {
                    double[][] tmp = allelicSeqEstimates.get(var).get(type);

                    stream.print(var + "\t" + type.toLowerCase() + "\t");
                    for (int i = 0; i < tmp.length; i++) {
                        stream.print(tmp[i][0] + "," + tmp[i][1]);

                        if (i < tmp.length - 1)
                            stream.print(";");
                    }

                    stream.println();
                }
            }
            stream.println("#");

            // Print sites map
            stream.println("#Map (chromosome number,locus number,reference nucleotide,alternative nucleotide)");
            stream.print(allelicSeqInfo.getSitesMapStr());
            stream.print("#");
        }

        stream.close();
    } // saveEstimatesTo

    public Runnable getModel(final String modelConfigFileName) throws JSONException, XMLParserException, IOException {
        Runnable model;

        if (modelConfigFileName.toLowerCase().endsWith(".json")) {
            model = new JSONParser(new HashMap<>()).parseFile(new File(modelConfigFileName));
        } else {
            try {
                model = new XMLParser(new HashMap<>(), null, false).parseFile(new File(modelConfigFileName));
            } catch (ParserConfigurationException | SAXException e) {
                throw new IllegalArgumentException(e);
            }
        }

        return model;
    } // getModel

    public void callVariants(
            int numOfModels,
            final String[] modelFileNames,
            final EstimatesTypeCollection[] modelEstimateTypes,
            Runnable[] models,
            double[] logPosteriors,
            double[] logTreeLikelihoods,
            PrintStream out,
            String outputVCFileName
    ) {
        assert numOfModels >= 1;
        assert modelFileNames != null;
        assert models != null;
        assert logPosteriors != null;

        List<Integer> maxIndex = new ArrayList<>();

        out.println();
        out.println(">>> Loading model configuration file(s)... [Thread " + Thread.currentThread().getId() + "]");
        out.flush();

        // Load configuration file and compute posterior
        for (int i = 0; i < numOfModels; i++) {
            out.println();

            // Load configuration file
            out.println(">>> Loading '" + modelFileNames[i] + "'... [Thread " + Thread.currentThread().getId() + "]");
            out.flush();
            try {
                models[i] = getModel(modelFileNames[i]);
            } catch (JSONException | XMLParserException | IOException e) {
                e.printStackTrace();
                out.println("Error parsing model configuration file: " + e.getMessage());
                return;
            }
            out.println(">>> Done.");
            out.flush();

            // Compute posterior
            try {
                models[i].run();
            } catch (Exception e) {
                e.printStackTrace();
                out.println("Error! Something is wrong when computing posterior for '" + modelFileNames[i] + "'");
                return;
            }
            logPosteriors[i] = ((VariantCall) models[i]).getLogPosterior();
            logTreeLikelihoods[i] = ((VariantCall) models[i]).getLogTreeLikelihood();
            out.println(">>> Log posterior: " + String.format(FLOAT_FORMAT, logPosteriors[i]));
            out.println(">>> Log tree likelihood: " + String.format(FLOAT_FORMAT, logTreeLikelihoods[i]));
            out.flush();
        }

        // If more than one configuration file provided, find the one with the largest posterior
        if (numOfModels > 1) {

            double max = -Double.MAX_VALUE;
            for (int i = 0; i < numOfModels; i++) {
                if (logPosteriors[i] > max) {
                    max = logPosteriors[i];

                    maxIndex.clear();
                    maxIndex.add(i);
                } else if (logPosteriors[i] == max) {
                    maxIndex.add(i);
                }
            }

            out.println();
            out.print(">>> The largest log posterior is " + String.format(FLOAT_FORMAT, max) + ", when using ");
            int j = 0;
            for (int i : maxIndex) {
                if (modelEstimateTypes[i].getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEDIAN) {
                    out.print("median estimates (log tree likelihood = " +
                            String.format(FLOAT_FORMAT, logTreeLikelihoods[i]) + ")");
                } else if (modelEstimateTypes[i].getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEAN) {
                    out.print("mean estimates (log tree likelihood = " +
                            String.format(FLOAT_FORMAT, logTreeLikelihoods[i]) + ")");
                } else if (modelEstimateTypes[i].getEstimatesType() == EstimatesTypeCollection.EstimatesType.MODE) {
                    out.print("mode estimates with " + modelEstimateTypes[i].getModeKDEType().toString().toLowerCase() +
                            " distribution as KDE kernel (log tree likelihood = " +
                            String.format(FLOAT_FORMAT, logTreeLikelihoods[i]) + ")");
                }

                if (j < maxIndex.size() - 1) {
                    out.print(" and ");
                }
            }
            out.println(".");

        } else
            maxIndex.add(0);

        // Set or update the output path
        System.setProperty("variant.calling.file.prefix",
                (System.getProperty("variant.calling.file.prefix") == null ?
                        "" :
                        System.getProperty("variant.calling.file.prefix")) +
                        FileNameProcessor.getRootPath(outputVCFileName));

        // Call variants and output to local file(s)
        out.println();
        for (int i = 0; i < maxIndex.size(); i++) {
            out.println(">>> Calling and writing variants [" + (i + 1) + "/" + maxIndex.size() + "]...");
            ((VariantCall) models[maxIndex.get(i)]).callVariantsAndGenerateVCF();
            out.println(">>> Done.");
        }

    } // callVariants


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("variantcaller", "");
    } // printUsage


    //***********************************************
    //*                 Main method                 *
    //***********************************************

    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // Define variables
        int numOfThreads;
        int burninPercentage;

        boolean isLoadCachedEstimates;
        String cachedEstimatesFileName = null;
        EstimatesTypeCollection.EstimatesType estimatesType = null;
        EstimatesTypeCollection.ModeKDEType modeKDEType = null;
        String mcmcSamplesFileName = null;
        String allelicInfoFileName = null;
        String gtAdoSamplesFileName = null;

        String modelConfigFileName = null;
        String inputTreeFileName = null;
        String outputVCFileName = null;

        boolean saveDetails;
        int cellThreshold = 1;
        boolean useMeanRate = true;
        boolean useOnlyBranchLength = true;

        String runConfigFileName = null;

        // No arguments provided, launch GUI
        if (args.length == 0) {
            Utils.loadUIManager();
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");
            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;
            if (url != null) {
                icon = new javax.swing.ImageIcon(url);
            }

            // Construct a new console
            new ConsoleApplication(null, null, icon, true);
            Log.info = System.out;
            Log.err = System.err;
            progressStream = System.out;

            // TODO: print some information here
            System.out.println("VariantCaller");

            VariantCallerDialog dialog = new VariantCallerDialog(false, new JFrame(), null);
            if (!dialog.showDialog("VariantCaller")) {
                return;
            }

            // Get parameters
            numOfThreads = dialog.getNumOfThreads();
            if (numOfThreads < 0) {
                Log.warning.println("number of threads is " + numOfThreads + " but should be positive. Setting it to one.");
                numOfThreads = 1;
            }
            if (numOfThreads > 1000) {
                Log.err.println("number of threads is " + numOfThreads + " but should be no larger than 1000.");
                return;
            }

            burninPercentage = dialog.getBurninText();
            if (burninPercentage < 0) {
                Log.warning.println("burnin percentage is " + burninPercentage + " but should be non-negative. Setting it to zero");
                burninPercentage = 0;
            }
            if (burninPercentage >= 100) {
                Log.err.println("burnin percentage is " + burninPercentage + " but should be less than 100.");
                return;
            }

            isLoadCachedEstimates = dialog.isLoadCachedEstimates();
            if (isLoadCachedEstimates) {
                // Load estimates from cached file

                cachedEstimatesFileName = dialog.getCachedEstimatesFileName();
                if (cachedEstimatesFileName == null) {
                    Log.err.println("No cached estimates file specified!");
                    return;
                }

            } else {
                // Compute estimates from log files

                estimatesType = dialog.getEstimatesType();
                progressStream.println("Use " + estimatesType.toString().toLowerCase() + " value of parameters computed " +
                        "from MCMC samples to perform variant calling.");

                if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE) {
                    modeKDEType = dialog.getModeKDEType();
                    progressStream.println("Use " + modeKDEType.toString().toLowerCase() + " distribution in the KDE of " +
                            "mode estimates.");
                }

                mcmcSamplesFileName = dialog.getMCMCSamplesFileName();
                if (mcmcSamplesFileName == null) {
                    Log.err.println("No MCMC samples log file specified!");
                    return;
                }

                allelicInfoFileName = dialog.getAllelicInfoFileName();

                gtAdoSamplesFileName = dialog.getGtAdoSamplesFileName();

            }

            modelConfigFileName = dialog.getModelConfigFileName();
            if (modelConfigFileName == null) {
                Log.err.println("No model configuration file specified!");
                return;
            }

            inputTreeFileName = dialog.getInputTreeFileName();
            if (inputTreeFileName == null) {
                Log.err.println("No best inferred tree file specified!");
                return;
            }

            outputVCFileName = dialog.getOutputVCFileName();

            saveDetails = dialog.getSaveDetails();

            cellThreshold = dialog.getCellThreshold();

            useMeanRate = dialog.isUseMedianRate();

            useOnlyBranchLength = dialog.isUseBranchLengthAndRate();

            try {
                if (isLoadCachedEstimates) {
                    new VariantCaller(
                            numOfThreads,
                            cachedEstimatesFileName,
                            modelConfigFileName,
                            inputTreeFileName,
                            outputVCFileName,
                            saveDetails,
                            cellThreshold,
                            useMeanRate,
                            useOnlyBranchLength
                    );
                } else {
                    new VariantCaller(
                            numOfThreads,
                            burninPercentage,
                            estimatesType,
                            modeKDEType,
                            mcmcSamplesFileName,
                            allelicInfoFileName,
                            gtAdoSamplesFileName,
                            modelConfigFileName,
                            inputTreeFileName,
                            outputVCFileName,
                            saveDetails,
                            cellThreshold,
                            useMeanRate,
                            useOnlyBranchLength
                    );
                }
            } catch (Exception e) {
                progressStream.println("Exception: " + e.getMessage());
            }

            progressStream.println("Finished - Quit program to exit.");
            while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

        } else {
            // TODO: print some information here
            Arguments arguments = new Arguments(
                    new Arguments.Option[]{
                            new Arguments.Option("help", "option to print this message -> OPTIONAL"),
                            new Arguments.IntegerOption("threads", 1, 1000, " specifies the number of threads (default 1); recommending to use more -> OPTIONAL"),
                            new Arguments.StringOption("prefix", "output_file_prefix", "specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL"),

                            // Run configuration file directly
                            new Arguments.StringOption("runconfig", "run_config_file", "specifies a modified configuration file and calls variants -> OPTIONAL, having priority over other flags except for -threads and -prefix"),

                            // Generate configuration file and run
                            new Arguments.IntegerOption("burnin", 0, 99, "specifies the percentage of samples to be considered as 'burn-in' -> MANDATORY (or -b)"),
                            new Arguments.IntegerOption("b", 0, 99, "the same as 'burn-in' -> MANDATORY (or -burnin)"),
                            new Arguments.StringOption("cached", "cached_estimates_file", "specifies the cached estimates file -> OPTIONAL, conflicting with -estimates, -kdedist, -mcmclog, -allelic"),
                            new Arguments.StringOption("estimates", new String[]{"mean", "median", "mode", "all"}, false, "specifies which kind of estimates of samples will be used to perform variant calling (one of 'mean', 'median' (default), 'mode', and 'all', where the last option compares likelihoods between all situations and chooses the largest one) -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("kdedist", new String[]{"gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"}, false, "specifies which KDE distribution will be used if mode estimates is selected (one of 'gaussian' (default), 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine', 'optcosine') -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("mcmclog", "mcmc_log_file", "specifies the log file containing MCMC samples -> OPTIONAL, conflicting with -cached"),
//                            new Arguments.StringOption("allelic", "allelic_seq_info_file", "specifies the allelic sequencing coverage and raw variance sampled during MCMC -> OPTIONAL, conflicting with -cached"),
//                            new Arguments.StringOption("gtado", "genotype_ado_state_samples_file", "specifies the genotypes and ado states sampled during MCMC -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("config", "config_file", "specifies the configuration file used to performing phylogenetic analisys (either xml or json) -> MANDATORY"),
                            new Arguments.StringOption("tree", "tree_file", "specifies the tree summarized from ScsTreeAnnotator -> MANDATORY"),
                            new Arguments.Option("details", "marks whether to save details, including inferred genotypes and ternary matrix, if specified -> OPTIONAL"),
                            new Arguments.IntegerOption("cells", 1, Integer.MAX_VALUE, "specifies the number of mutated cells used to filter invariant sites, default 1 -> OPTIONAL"),
                            new Arguments.Option("medianrate", "use median rate rather than mean rate from the input tree -> OPTIONAL"),
                            new Arguments.Option("userate", "use both branch length and rate from the input tree in variant calling -> OPTIONAL")
                    }
            );

            // Parse command line arguments
            try {
                arguments.parseArguments(args);
            } catch (Arguments.ArgumentException e) {
                progressStream.println(e);
                printUsage(arguments);
                System.exit(1);
            }

            // Print help message
            if (arguments.hasOption("help") || arguments.hasOption("h")) {
                printUsage(arguments);
                System.exit(0);
            }

            // Set numOfThreads
            numOfThreads = 1;
            if (arguments.hasOption("threads")) {
                numOfThreads = arguments.getIntegerOption("threads");
            }
            if (numOfThreads > 1000) {
                Log.err.println("number of threads is " + numOfThreads + " but should be no larger than 1000.");
                System.exit(1);
            }

            if (arguments.hasOption("prefix"))
                System.setProperty("variant.calling.file.prefix", arguments.getStringOption("prefix").trim());

            if (arguments.hasOption("runconfig")) {

                runConfigFileName = arguments.getStringOption("runconfig");

                Log.warning.println("Run provided configuration file: '" + runConfigFileName + "'. " +
                        "Any other flags will be ignored except for '-threads' and '-prefix'.");

                try {
                    new VariantCaller(
                            numOfThreads,
                            runConfigFileName
                    );
                } catch (Exception e) {
                    e.printStackTrace();
                }

            } else {

                // Set burninPercentage
                burninPercentage = -1;
                if (arguments.hasOption("burnin")) {
                    burninPercentage = arguments.getIntegerOption("burnin");
                } else if (arguments.hasOption("b")) {
                    burninPercentage = arguments.getIntegerOption("b");
                } else {
                    printUsage(arguments);
                    Log.err.println("Option -burnin or -b is missing.");
                    System.exit(1);
                }
                if (burninPercentage >= 100) {
                    Log.err.println("burnin percentage is " + burninPercentage + " but should be less than 100.");
                    System.exit(1);
                }

                // Set cachedEstimatesFileName
                if (arguments.hasOption("cached")) {
                    // Load estimates from cached file

                    isLoadCachedEstimates = true;

                    // Sanity check
                    if (arguments.hasOption("estimates") ||
                            arguments.hasOption("kdedist") ||
                            arguments.hasOption("mcmclog") ||
                            arguments.hasOption("allelic") ||
                            arguments.hasOption("gtado")) {
                        Log.err.println("-cached flag conflicts with -estimates, -kdedist, -mcmclog, and -allelic");
                        System.exit(1);
                    }

                    cachedEstimatesFileName = arguments.getStringOption("cached");
                } else {
                    // Compute estimates from log files

                    isLoadCachedEstimates = false;

                    // Set estimatesType
                    estimatesType = EstimatesTypeCollection.EstimatesType.MEDIAN;
                    if (arguments.hasOption("estimates")) {
                        String value = arguments.getStringOption("estimates");
                        if (value.equalsIgnoreCase("mean")) {
                            estimatesType = EstimatesTypeCollection.EstimatesType.MEAN;
                        } else if (value.equalsIgnoreCase("median")) {
                            estimatesType = EstimatesTypeCollection.EstimatesType.MEDIAN;
                        } else if (value.equalsIgnoreCase("mode")) {
                            estimatesType = EstimatesTypeCollection.EstimatesType.MODE;
                        } else if (value.equalsIgnoreCase("all")) {
                            estimatesType = EstimatesTypeCollection.EstimatesType.ALL;
                        }
                    }

                    // Set modeKDEType
                    if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE) {
                        modeKDEType = EstimatesTypeCollection.ModeKDEType.GAUSSIAN;
                        if (arguments.hasOption("kdedist")) {
                            String value = arguments.getStringOption("kdedist");
                            if (value.equalsIgnoreCase("gaussian")) {
                                modeKDEType = EstimatesTypeCollection.ModeKDEType.GAUSSIAN;
                            }
                        }
                    }

                    // Set mcmcSamplesFileName
                    if (arguments.hasOption("mcmclog")) {
                        mcmcSamplesFileName = arguments.getStringOption("mcmclog");
                    } else {
                        printUsage(arguments);
                        Log.err.println("Option -mcmclog is missing.");
                        System.exit(1);
                    }

                    // Set allelicInfoFileName
                    if (arguments.hasOption("allelic")) {
                        allelicInfoFileName = arguments.getStringOption("allelic");
                    }

                    // Set gtAdoSamplesFileName
                    if (arguments.hasOption("gtado")) {
                        gtAdoSamplesFileName = arguments.getStringOption("gtado");
                    }

                }

                // Set modelConfigFileName
                if (arguments.hasOption("config")) {
                    modelConfigFileName = arguments.getStringOption("config");
                } else {
                    printUsage(arguments);
                    Log.err.println("Option -estimates is missing.");
                    System.exit(1);
                }

                // Set inputTreeFileName
                if (arguments.hasOption("tree")) {
                    inputTreeFileName = arguments.getStringOption("tree");
                } else {
                    printUsage(arguments);
                    Log.err.println("Option -tree is missing.");
                    System.exit(1);
                }

                // Set saveToTernary
                saveDetails = arguments.hasOption("details");

                // Set cellThreshold
                if (arguments.hasOption("cells")) {
                    cellThreshold = arguments.getIntegerOption("cells");
                }

                // Set useMeanRate
                if (arguments.hasOption("medianrate"))
                    useMeanRate = false;

                // Set useOnlyBranchLength
                if (arguments.hasOption("userate"))
                    useOnlyBranchLength = false;

                try {
                    if (isLoadCachedEstimates) {
                        new VariantCaller(
                                numOfThreads,
                                cachedEstimatesFileName,
                                modelConfigFileName,
                                inputTreeFileName,
                                outputVCFileName,
                                saveDetails,
                                cellThreshold,
                                useMeanRate,
                                useOnlyBranchLength
                        );
                    } else {
                        new VariantCaller(
                                numOfThreads,
                                burninPercentage,
                                estimatesType,
                                modeKDEType,
                                mcmcSamplesFileName,
                                allelicInfoFileName,
                                gtAdoSamplesFileName,
                                modelConfigFileName,
                                inputTreeFileName,
                                outputVCFileName,
                                saveDetails,
                                cellThreshold,
                                useMeanRate,
                                useOnlyBranchLength
                        );
                    }
                } catch (IOException e) {
                    throw e;
                } catch (Exception e) {
                    e.printStackTrace();
                }

            }

        }

        System.exit(0);

    } // main

}
