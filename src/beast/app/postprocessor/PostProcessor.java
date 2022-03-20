package beast.app.postprocessor;

import beast.app.tools.LogCombiner;
import beast.app.treeannotator.ScsTreeAnnotator;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.app.variantcaller.EstimatesTypeCollection;
import beast.app.variantcaller.VariantCaller;
import beast.core.util.Log;
import beast.util.FileNameProcessor;
import jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class PostProcessor {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public static ExecutorService g_exec = Executors.newFixedThreadPool(1);

    public static PrintStream progressStream = Log.err;


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("postprocessor", "");
    } // printUsage


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // Define variables
        // Overall
        int burninPercentage;

        // For ScsTreeAnnotator
        String targetTreeFileName = null;
        String inputTreesFileName = null;
        String outputTreeFileName = null;
        String outputSimpleTreeFileName = null;
        double posteriorLimit;
        double hpd2D;
        boolean lowMem;
        boolean saveSimpleTree = false;
        boolean hasTrunk = true;
        boolean forceIntegerToDiscrete = false;
        boolean processSA = true;
        boolean SAmode = false;
        boolean processBivariateAttributes = true;

        // For VariantCaller
        int numOfThreads;
        boolean isLoadCachedEstimates;
        String cachedEstimatesFileName = null;
        EstimatesTypeCollection.EstimatesType estimatesType = null;
        EstimatesTypeCollection.ModeKDEType modeKDEType = null;
        String mcmcSamplesFileName = null;
        String allelicInfoFileName = null;
        String gtAdoSamplesFileName = null;
        String modelConfigFileName = null;
        String outputVCFileName = null;
        boolean saveDetails;
        int cellThreshold = 1;
        boolean useMeanRate = true;
        boolean useOnlyBranchLength = true;

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

            // Console application
            new ConsoleApplication(null, null, icon, true);
            Log.info = System.out;
            Log.err = System.err;

            // TODO: print some information here
            progressStream.println("PostProcessor");

            PostProcessorDialog dialog = new PostProcessorDialog(new JFrame());
            if (!dialog.showDialog("PostProcessor")) {
                return;
            }

            // Get parameters for ScsTreeAnnotator
            burninPercentage = dialog.getBurninText();
            if (burninPercentage < 0) {
                progressStream.println("burnin percentage is " + burninPercentage + " but should be non-negative. Setting it to zero");
                burninPercentage = 0;
            }
            if (burninPercentage >= 100) {
                progressStream.println("burnin percentage is " + burninPercentage + " but should be less than 100.");
                return;
            }

            posteriorLimit = dialog.getTADialog().getPosteriorLimit();
            hpd2D = 0.80;
            ScsTreeAnnotator.Target targetOption = dialog.getTADialog().getTargetOption();
            ScsTreeAnnotator.HeightsSummary heightsOption = dialog.getTADialog().getHeightsOption();

            targetTreeFileName = dialog.getTADialog().getTargetFileName();
            if (targetOption == ScsTreeAnnotator.Target.USER_TARGET_TREE && targetTreeFileName == null) {
                progressStream.println("No target file specified");
                return;
            }

            inputTreesFileName = dialog.getTADialog().getInputFileName();
            if (inputTreesFileName == null) {
                progressStream.println("No input tree file specified");
                return;
            }

            outputTreeFileName = dialog.getTADialog().getOutputFileName();
            if (outputTreeFileName == null) {
                progressStream.println("No output tree file specified");
                return;
            }

            lowMem = dialog.getTADialog().useLowMem();
            saveSimpleTree = dialog.getTADialog().saveSimpleTree();
            hasTrunk = dialog.getTADialog().hasTrunk();

            if (saveSimpleTree) {
                String delimiter = String.valueOf(File.separatorChar);
                String[] tmp1 = outputTreeFileName.split(delimiter);
                String[] tmp2 = tmp1[tmp1.length - 1].split("\\.");
                tmp1[tmp1.length - 1] = tmp2[0] + "_simple." + tmp2[1];
                outputSimpleTreeFileName = String.join(delimiter, tmp1);
            }

            // Get parameters for VariantCaller
            numOfThreads = dialog.getVCDialog().getNumOfThreads();
            if (numOfThreads < 0) {
                progressStream.println("number of threads is " + numOfThreads + " but should be positive. Setting it to one.");
                numOfThreads = 1;
            }
            if (numOfThreads > 1000) {
                progressStream.println("number of threads is " + numOfThreads + " but should be no larger than 1000.");
                return;
            }

            isLoadCachedEstimates = dialog.getVCDialog().isLoadCachedEstimates();
            if (isLoadCachedEstimates) {
                // Load estimates from cached file

                cachedEstimatesFileName = dialog.getVCDialog().getCachedEstimatesFileName();
                if (cachedEstimatesFileName == null) {
                    progressStream.println("No cached estimates file specified!");
                    return;
                }

            } else {
                // Compute estimates from log files

                estimatesType = dialog.getVCDialog().getEstimatesType();
                progressStream.println("Use " + estimatesType.toString().toLowerCase() + " value of parameters computed " +
                        "from MCMC samples to perform variant calling.");

                if (estimatesType == EstimatesTypeCollection.EstimatesType.MODE) {
                    modeKDEType = dialog.getVCDialog().getModeKDEType();
                    progressStream.println("Use " + modeKDEType.toString().toLowerCase() + " distribution in the KDE of " +
                            "mode estimates.");
                }

                mcmcSamplesFileName = dialog.getVCDialog().getMCMCSamplesFileName();
                if (mcmcSamplesFileName == null) {
                    progressStream.println("No MCMC samples log file specified!");
                    return;
                }

                allelicInfoFileName = dialog.getVCDialog().getAllelicInfoFileName();

                gtAdoSamplesFileName = dialog.getVCDialog().getGtAdoSamplesFileName();

            }

            modelConfigFileName = dialog.getVCDialog().getModelConfigFileName();
            if (modelConfigFileName == null) {
                progressStream.println("No model configuration file specified!");
                return;
            }

            outputVCFileName = dialog.getVCDialog().getOutputVCFileName();
            saveDetails = dialog.getVCDialog().getSaveDetails();
            cellThreshold = dialog.getVCDialog().getCellThreshold();

            useMeanRate = dialog.getVCDialog().isUseMedianRate();

            useOnlyBranchLength = dialog.getVCDialog().isUseBranchLengthAndRate();

            try {
                new ScsTreeAnnotator(
                        processBivariateAttributes,
                        forceIntegerToDiscrete,
                        processSA,
                        SAmode,
                        saveSimpleTree,
                        hasTrunk,
                        burninPercentage,
                        lowMem,
                        heightsOption,
                        posteriorLimit,
                        hpd2D,
                        targetOption,
                        targetTreeFileName,
                        inputTreesFileName,
                        outputTreeFileName,
                        outputSimpleTreeFileName
                );
            } catch (Exception ex) {
                progressStream.println("Exception: " + ex.getMessage());
            }

            try {
                if (isLoadCachedEstimates) {
                    new VariantCaller(
                            numOfThreads,
                            cachedEstimatesFileName,
                            modelConfigFileName,
                            outputTreeFileName,
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
                            outputTreeFileName,
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

            Arguments arguments = new Arguments(
                    new Arguments.Option[]{
                            new Arguments.Option("help", "option to print this message -> OPTIONAL"),
                            new Arguments.IntegerOption("burnin", 0, 99, "the percentage of states to be considered as 'burn-in' -> MANDATORY, the same as -b"),
                            new Arguments.StringOption("prefix", "output_file_prefix", "specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL"),

                            // ScsTreeAnnotator
                            //new Arguments.StringOption("target", new String[] { "maxclade", "maxtree" }, false, "an option of 'maxclade' or 'maxtree'"),
                            new Arguments.StringOption("heights", new String[]{"keep", "median", "mean", "ca"}, false,
                                    "an option of 'keep' (default), 'median', 'mean' or 'ca' -> OPTIONAL"),
                            // allow -b as burnin option, just like other apps
                            new Arguments.IntegerOption("b", 0, 99, "the percentage of states to be considered as 'burn-in' -> MANDATORY, the same as -burnin"),
                            new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annotated -> OPTIONAL"),
                            new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated -> OPTIONAL"),
                            new Arguments.Option("forceDiscrete", "forces integer traits to be treated as discrete traits. -> OPTIONAL"),
                            new Arguments.Option("lowMem", "use less memory, which is a bit slower -> OPTIONAL"),
                            new Arguments.RealOption("hpd2D", "the HPD interval to be used for the bivariate traits -> OPTIONAL"),
                            new Arguments.Option("nohpd2D", "suppress calculation of HPD intervals for the bivariate traits -> OPTIONAL"),
                            new Arguments.Option("noSA", "interpret the tree set as begin from a not being from a sampled ancestor analysis, even if there are zero branch lengths in the tree set -> OPTIONAL"),
                            new Arguments.Option("simpleTree", "simple output tree only containing tree heights (the output file will be labeled with 'simple'), default disabled -> OPTIONAL"),
                            new Arguments.Option("noTrunk", "trees processed do not have a trunk connecting the root of samples and the root of entire tree, default enabled -> OPTIONAL"),
                            new Arguments.StringOption("trees", "input_trees_file", "specifies the file containing sampled trees -> MANDATORY"),
                            new Arguments.StringOption("outTree", "summarized_best_tree", "specifies the file containing the summarized best tree -> MANDATORY"),

                            // VariantCaller
                            new Arguments.IntegerOption("threads", 1, 1000, " specifies the number of threads (default 1); recommending to use more -> OPTIONAL"),
                            new Arguments.StringOption("cached", "cached_estimates_file", "specifies the cached estimates file -> OPTIONAL, conflicting with -estimates, -kdedist, -mcmclog, -allelic"),
                            new Arguments.StringOption("estimates", new String[]{"mean", "median", "mode", "all"}, false, "specifies which kind of estimates of samples will be used to perform variant calling (one of 'mean', 'median' (default), 'mode', and 'all', where the last option compares likelihoods between all situations and chooses the largest one) -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("kdedist", new String[]{"gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"}, false, "specifies which KDE distribution will be used if mode estimates is selected (one of 'gaussian' (default), 'epanechnikov', 'rectangular', 'triangular', 'biweight', 'cosine', 'optcosine') -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("mcmclog", "mcmc_log_file", "specifies the log file containing MCMC samples -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("allelic", "allelic_seq_info_file", "specifies the allelic sequencing coverage and raw variance sampled during MCMC -> OPTIONAL, conflicting with -cached"),
                            new Arguments.StringOption("config", "config_file", "specifies the configuration file used to performing phylogenetic analisys (either xml or json) -> MANDATORY"),
                            new Arguments.StringOption("out", "vcf_file", "specifies a vcf file containing the output variant calling results -> OPTIONAL"),
                            new Arguments.Option("details", "marks whether to save details, including inferred genotypes and ternary matrix, if specified -> OPTIONAL"),
                            new Arguments.IntegerOption("cells", 1, Integer.MAX_VALUE, "specifies the number of mutated cells used to filter invariant sites, default 1 -> OPTIONAL"),
                            new Arguments.Option("medianrate", "use median rate rather than mean rate from the input tree -> OPTIONAL"),
                            new Arguments.Option("userate", "use both branch length and rate from the input tree in variant calling -> OPTIONAL")
                    });

            try {
                arguments.parseArguments(args);
            } catch (Arguments.ArgumentException ae) {
                progressStream.println(ae);
                printUsage(arguments);
                System.exit(1);
            }

            if (arguments.hasOption("help")) {
                printUsage(arguments);
                System.exit(0);
            }

            if (arguments.hasOption("prefix")) {
                System.setProperty("tree.annotation.file.prefix", arguments.getStringOption("prefix").trim());
                System.setProperty("variant.calling.file.prefix", arguments.getStringOption("prefix").trim());
            }

            // ScsTreeAnnotator
            if (arguments.hasOption("nohpd2D")) {
                processBivariateAttributes = false;
            }

            if (arguments.hasOption("noSA")) {
                processSA = false;
            }

            if (arguments.hasOption("forceDiscrete")) {
                Log.info.println("  Forcing integer traits to be treated as discrete traits.");
                forceIntegerToDiscrete = true;
            }

            lowMem = false;
            if (arguments.hasOption("lowMem")) {
                lowMem = true;
            }

            ScsTreeAnnotator.HeightsSummary heights = ScsTreeAnnotator.HeightsSummary.CA_HEIGHTS;
            if (arguments.hasOption("heights")) {
                String value = arguments.getStringOption("heights");
                if (value.equalsIgnoreCase("mean")) {
                    heights = ScsTreeAnnotator.HeightsSummary.MEAN_HEIGHTS;
                } else if (value.equalsIgnoreCase("median")) {
                    heights = ScsTreeAnnotator.HeightsSummary.MEDIAN_HEIGHTS;
                } else if (value.equalsIgnoreCase("ca")) {
                    heights = ScsTreeAnnotator.HeightsSummary.CA_HEIGHTS;
                    Log.info.println("Please cite: Heled and Bouckaert: Looking for trees in the forest:\n" +
                            "summary tree from posterior samples. BMC Evolutionary Biology 2013 13:221.");
                }
            }

            burninPercentage = -1;
            if (arguments.hasOption("burnin")) {
                burninPercentage = arguments.getIntegerOption("burnin");
            } else if (arguments.hasOption("b")) {
                burninPercentage = arguments.getIntegerOption("b");
            }
            if (burninPercentage >= 100) {
                Log.err.println("burnin percentage is " + burninPercentage + " but should be less than 100.");
                System.exit(1);
            }

            posteriorLimit = 0.0;
            if (arguments.hasOption("limit")) {
                posteriorLimit = arguments.getRealOption("limit");
            }

            hpd2D = 0.80;
            if (arguments.hasOption("hpd2D")) {
                hpd2D = arguments.getRealOption("hpd2D");
                if (hpd2D <= 0 || hpd2D >= 1) {
                    Log.err.println("hpd2D is a fraction and should be in between 0.0 and 1.0.");
                    System.exit(1);
                }
                processBivariateAttributes = true;
            }

            ScsTreeAnnotator.Target target = ScsTreeAnnotator.Target.MAX_CLADE_CREDIBILITY;
            if (arguments.hasOption("target")) {
                target = ScsTreeAnnotator.Target.USER_TARGET_TREE;
                targetTreeFileName = arguments.getStringOption("target");
            }

            if (arguments.hasOption("simpleTree")) {
                saveSimpleTree = true;
            }

            if (arguments.hasOption("noTrunk")) {
                Log.warning.println("are you sure your sampled trees do not contain an extra trunk?");
                hasTrunk = false;
            }

            if (arguments.hasOption("trees")) {
                inputTreesFileName = arguments.getStringOption("trees");
            } else {
                Log.err.println("-trees must be specified.");
                printUsage(arguments);
                System.exit(1);
            }

            if (arguments.hasOption("outTree")) {
                outputTreeFileName = arguments.getStringOption("outTree");
            } else {
                Log.err.println("-outTree must be specified.");
                printUsage(arguments);
                System.exit(1);
            }

            if (saveSimpleTree) {
                assert outputTreeFileName != null;
                outputSimpleTreeFileName = FileNameProcessor.getRootPath(outputTreeFileName) +
                        FileNameProcessor.getBaseName(outputTreeFileName) +
                        "_simple" +
                        FileNameProcessor.getSuffix(outputTreeFileName);
            }

            // VariantCalling
            // Set numOfThreads
            numOfThreads = 1;
            if (arguments.hasOption("threads")) {
                numOfThreads = arguments.getIntegerOption("threads");
            }
            if (numOfThreads > 1000) {
                Log.err.println("number of threads is " + numOfThreads + " but should be no larger than 1000.");
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
                        arguments.hasOption("allelic")) {
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
                        } /*else if (value.equalsIgnoreCase("epanechnikov")) {
                            modeKDEType = ModeKDEType.EPANECHNIKOV;
                        } else if (value.equalsIgnoreCase("rectangular")) {
                            modeKDEType = ModeKDEType.RECTANGULAR;
                        } else if (value.equalsIgnoreCase("triangular")) {
                            modeKDEType = ModeKDEType.TRIANGULAR;
                        } else if (value.equalsIgnoreCase("biweight")) {
                            modeKDEType = ModeKDEType.BIWEIGHT;
                        } else if (value.equalsIgnoreCase("cosine")) {
                            modeKDEType = ModeKDEType.COSINE;
                        } else if (value.equalsIgnoreCase("optcosine")) {
                            modeKDEType = ModeKDEType.OPTCOSINE;
                        }*/
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

            // Set outputVCFileName
            if (arguments.hasOption("out")) {
                outputVCFileName = arguments.getStringOption("out");
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

            // ScsTreeAnnotator
            try {
                new ScsTreeAnnotator(
                        processBivariateAttributes,
                        forceIntegerToDiscrete,
                        processSA,
                        SAmode,
                        saveSimpleTree,
                        hasTrunk,
                        burninPercentage,
                        lowMem,
                        heights,
                        posteriorLimit,
                        hpd2D,
                        target,
                        targetTreeFileName,
                        inputTreesFileName,
                        outputTreeFileName,
                        outputSimpleTreeFileName
                );
            } catch (IOException e) {
                throw e;
            } catch (Exception e) {
                e.printStackTrace();
            }

            if (System.getProperty("tree.annotation.file.prefix") != null) {
                outputTreeFileName = System.getProperty("tree.annotation.file.prefix") + outputTreeFileName;
            }

            // VariantCaller
            try {
                if (isLoadCachedEstimates) {
                    new VariantCaller(
                            numOfThreads,
                            cachedEstimatesFileName,
                            modelConfigFileName,
                            outputTreeFileName,
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
                            outputTreeFileName,
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

        g_exec.shutdown();
        g_exec.shutdownNow();
        System.exit(0);

    } // main

}
