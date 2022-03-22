package beast.app.geneannotator;

import beast.app.tools.LogCombiner;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.app.utils.ChromosomeLabel;
import beast.core.util.Log;
import beast.evolution.alignment.VariantSiteInfo;
import beast.evolution.substitutionmodel.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.NexusParser;
import beast.util.NotSingleException;
import beast.util.TreeParser;
import jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static beast.evolution.variantsinfo.GenericVariantsInfoVCF.META_DATA_GENOTYPES;
import static beast.util.FileNameProcessor.getBaseName;
import static beast.util.FileNameProcessor.getRootPath;
import static beast.util.TreeUtils.processMetaData;

public class GeneAnnotator {

    private final static String GENOTYPES_FORMAT_ERROR = "The genotypes parsed from the input tree can only be numbers.";

    private final static String METADATA_INDEX = "index";
    private final static String METADATA_CHR = "chr";
    private final static String METADATA_POS = "pos";
    private final static String METADATA_REF_NUC = "ref_nuc";
    private final static String METADATA_ALT_NUC = "alt_nuc";
    private final static String METADATA_EVENT_TYPE = "event_type";
    private final static String METADATA_GENE = "gene";
    private final static String METADATA_ISA_SUFFIX = "_isa";
    private final static String METADATA_FSA_SUFFIX = "_fsa";

    private final static String OUTPUT_SUFFIX = ".gene_tree";

    private final static String GENE_NAME_HEADER = "Gene.refGene";

    static PrintStream progressStream = Log.err;

    private final NodeGenesInfo[] nodeGenesInfos;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public GeneAnnotator(
            int substModelLabel,
            String treeFileName,
            String snvFileName,
            String mutationMapFileName,
            boolean mapFromAnnovar,
            String mapSeparator,
            String filteringGenesFileName,
            String filteringGenesSeparator,
            int filteringGenesColIndex,
            String outputFileName
    ) throws IOException, NotSingleException {
        ScsSubstitutionModelBase substModel;
        Tree tree;
        List<VariantSiteInfo> snvSites;
        List<String> filteringGenes = null;
        Map<ChromosomeLabel, List<MutationMap>> mutationMap;

        // 1. Get substitution model
        substModel = getSubstModel(substModelLabel);

        // 2. Get input tree
        tree = getTree(treeFileName);

        // 3. Get variant sites info
        snvSites = getSNVSitesInfo(snvFileName);

        // 4. Get filtering genes
        if (filteringGenesFileName != null)
            filteringGenes = getFilteringGenes(filteringGenesFileName, filteringGenesSeparator, filteringGenesColIndex);

        // 5. Get mutation map
        if (mapFromAnnovar)
            mutationMap = getMutationMapFromAnnovar(mutationMapFileName, mapSeparator, filteringGenes);
        else
            mutationMap = getMutationMap(mutationMapFileName, filteringGenes);

        // 6. Initialize some variables
        this.nodeGenesInfos = new NodeGenesInfo[tree.getNodeCount()];

        // 7. Collect genes for all nodes
        collectGenes4Tree(
                tree.getRoot(),
                META_DATA_GENOTYPES,
                null,
                substModel,
                snvSites,
                mutationMap
        );

        // 8. Process genes for all nodes
        processGenes();

        // 9. Annotate the tree
        for (Node node : tree.getNodesAsArray()) {
            node.removeMetaData(META_DATA_GENOTYPES);

            final NodeGenesInfo info = nodeGenesInfos[node.getNr()];

            if (info == null) continue;

            node.setMetaData(METADATA_INDEX, node.getNr());

            if (info.getFullGenes() != null && info.getFullGenes().size() > 0) {
                node.setMetaData(METADATA_CHR, info.getFullChr().toArray(new String[0]));
                node.setMetaData(METADATA_POS, info.getFullPos().toArray(new Long[0]));
                node.setMetaData(METADATA_REF_NUC, info.getFullRefNuc().toArray(new String[0]));
                node.setMetaData(METADATA_ALT_NUC, info.getFullAltNuc().toArray(new String[0]));
                node.setMetaData(METADATA_EVENT_TYPE, info.getFullEventType().toArray(new String[0]));
                node.setMetaData(METADATA_GENE, info.getFullGeneName().toArray(new String[0]));
            }

            if (info.getISAGenes() != null && info.getISAGenes().size() > 0) {
                node.setMetaData(METADATA_CHR + METADATA_ISA_SUFFIX, info.getISAChr().toArray(new String[0]));
                node.setMetaData(METADATA_POS + METADATA_ISA_SUFFIX, info.getISAPos().toArray(new Long[0]));
                node.setMetaData(METADATA_REF_NUC + METADATA_ISA_SUFFIX, info.getISARefNuc().toArray(new String[0]));
                node.setMetaData(METADATA_ALT_NUC + METADATA_ISA_SUFFIX, info.getISAAltNuc().toArray(new String[0]));
                node.setMetaData(METADATA_EVENT_TYPE + METADATA_ISA_SUFFIX, info.getISAEventType().toArray(new String[0]));
                node.setMetaData(METADATA_GENE + METADATA_ISA_SUFFIX, info.getISAGeneName().toArray(new String[0]));
            }

            if (info.getFSAGenes() != null && info.getFSAGenes().size() > 0) {
                node.setMetaData(METADATA_CHR + METADATA_FSA_SUFFIX, info.getFSAChr().toArray(new String[0]));
                node.setMetaData(METADATA_POS + METADATA_FSA_SUFFIX, info.getFSAPos().toArray(new Long[0]));
                node.setMetaData(METADATA_REF_NUC + METADATA_FSA_SUFFIX, info.getFSARefNuc().toArray(new String[0]));
                node.setMetaData(METADATA_ALT_NUC + METADATA_FSA_SUFFIX, info.getFSAAltNuc().toArray(new String[0]));
                node.setMetaData(METADATA_EVENT_TYPE + METADATA_FSA_SUFFIX, info.getFSAEventType().toArray(new String[0]));
                node.setMetaData(METADATA_GENE + METADATA_FSA_SUFFIX, info.getFSAGeneName().toArray(new String[0]));
            }
        }

        // 10. Save the annotated tree to disk
        processMetaData(tree.getRoot());
        final PrintStream treeOut = new PrintStream(outputFileName);
        tree.init(treeOut);
        treeOut.println();
        treeOut.print("tree TREE1 = ");
        treeOut.print(tree.getRoot().toSortedNewick(new int[1], true));
        treeOut.println(";");
        tree.close(treeOut);
        treeOut.println();
        treeOut.close();

        progressStream.println("Successful!");
        progressStream.println("Annotated tree has been save to " + outputFileName);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    private ScsSubstitutionModelBase getSubstModel(final int substModelLabel) {
        switch (substModelLabel) {
            case 0:
                return new ScsFiniteMuExtendedModel();
            case 1:
                return new ScsFiniteMuModel();
            case 2:
                return new ScsFiniteMuDelModel();
            case 3:
                return new ScsFiniteMuDelInsModel();
            default:
                throw new RuntimeException("Error! Invalid label for substitution model: " + substModelLabel);
        }
    } // getSubstModel

    private Tree getTree(final String treeFileName) throws IOException, NotSingleException {
        List<Tree> parsedTrees;
        boolean isNexus = true;
        BufferedReader fin;
        String str;

        // Nexus tree or Newick tree?
        fin = new BufferedReader(new FileReader(treeFileName));
        if (!fin.ready()) {
            throw new IOException(treeFileName + " appears empty.");
        }
        str = fin.readLine();
        if (!str.toUpperCase().trim().startsWith("#NEXUS")) {
            // Newick tree found
            isNexus = false;
        }

        // Parse tree
        if (isNexus) {
            NexusParser nexusParser = new NexusParser();
            nexusParser.parseFile(new File(treeFileName));
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

    private List<VariantSiteInfo> getSNVSitesInfo(final String snvFileName) throws IOException {
        List<VariantSiteInfo> results = new ArrayList<>();
        BufferedReader fin = new BufferedReader(new FileReader(snvFileName));

        if (!fin.ready()) {
            throw new IOException(snvFileName + " appears empty.");
        }

        while (fin.ready()) {
            String[] comp = fin.readLine().trim().split("\\s+");

            try {
                results.add(
                        new VariantSiteInfo(
                                comp[0],
                                Long.parseLong(comp[1]),
                                comp[2],
                                comp[3].split(",")
                        )
                );
            } catch (ArithmeticException e) {
                throw new IllegalArgumentException("Error! Make sure the position is a number: " +
                        String.join(";", comp));
            }
        }

        fin.close();
        return results;
    } // getSNVSitesInfo

    private List<String> getFilteringGenes(
            String filteringGenesFileName,
            String filteringGenesSeparator,
            int filteringGenesColIndex
    ) throws IOException {
        if (filteringGenesFileName == null)
            return null;

        List<String> results = new ArrayList<>();
        BufferedReader fin = new BufferedReader(new FileReader(filteringGenesFileName));

        if (!fin.ready())
            throw new IOException(filteringGenesFileName + " appears empty.");

        int lineNumber = 0;
        while (fin.ready()) {
            final String[] line = fin.readLine().trim().replaceAll("\"", "").split(filteringGenesSeparator);

            if (lineNumber > 0) {
                final String gene = line[filteringGenesColIndex];

                if (!results.contains(gene))
                    results.add(gene);
            }

            lineNumber++;
        }

        fin.close();

        return results;
    }

    private Map<ChromosomeLabel, List<MutationMap>> getMutationMapFromAnnovar(
            String mutationMapFileName,
            String mapSeparator,
            List<String> filteringGenes
    ) throws IOException {
        Map<ChromosomeLabel, List<MutationMap>> results = new HashMap<>();
        BufferedReader fin = new BufferedReader(new FileReader(mutationMapFileName));

        if (!fin.ready())
            throw new IOException(mutationMapFileName + " appears empty.");

        int lineNumber = 0;
        int geneNameIndex = -1;
        while (fin.ready()) {
            final String[] line = fin.readLine().trim().replaceAll("\"", "").split(mapSeparator);

            if (lineNumber == 0) {
                // process the header

                for (int i = 0; i < line.length; i++) {
                    if (line[i].equalsIgnoreCase(GENE_NAME_HEADER))
                        geneNameIndex = i;
                }

                if (geneNameIndex == -1)
                    throw new IOException("Error! The mutation map does not contain a column named " + GENE_NAME_HEADER);
            } else {
                // process each line

                MutationMap mm = new MutationMap(new String[]{line[0].trim(), line[1].trim(), line[2].trim(), line[geneNameIndex].trim()});
                final ChromosomeLabel key = mm.getChromosome();

                if (filteringGenes == null || filteringGenes.contains(mm.getGeneName())) {
                    if (results.containsKey(key))
                        results.get(key).add(mm);
                    else
                        results.put(key, Stream.of(mm).collect(Collectors.toList()));
                }
            }

            lineNumber++;
        }

        final MutationMap.ListComparator comparator = new MutationMap.ListComparator();
        for (ChromosomeLabel chromosomeLabel : results.keySet()) {
            results.get(chromosomeLabel).sort(comparator);
        }

        fin.close();
        return results;
    } // getMutationMapFromAnnovar

    private Map<ChromosomeLabel, List<MutationMap>> getMutationMap(
            String mutationMapFileName,
            List<String> filteringGenes
    ) throws IOException {
        Map<ChromosomeLabel, List<MutationMap>> results = new HashMap<>();
        BufferedReader fin = new BufferedReader(new FileReader(mutationMapFileName));

        if (!fin.ready())
            throw new IOException(mutationMapFileName + " appears empty.");

        while (fin.ready()) {
            String line = fin.readLine().trim();
            if (!line.startsWith("#")) {
                MutationMap mm = new MutationMap(line.split("\\s+"));
                final ChromosomeLabel key = mm.getChromosome();

                if (filteringGenes == null || filteringGenes.contains(mm.getGeneName())) {
                    if (results.containsKey(key))
                        results.get(key).add(mm);
                    else
                        results.put(key, Stream.of(mm).collect(Collectors.toList()));
                }
            }
        }

        final MutationMap.ListComparator comparator = new MutationMap.ListComparator();
        for (ChromosomeLabel chromosomeLabel : results.keySet()) {
            results.get(chromosomeLabel).sort(comparator);
        }

        fin.close();
        return results;
    } // getMutationMap

    private void collectGenes4Tree(
            Node node,
            String genotypeKey,
            Object parentGenotypes,
            final ScsSubstitutionModelBase substModel,
            final List<VariantSiteInfo> snvSites,
            final Map<ChromosomeLabel, List<MutationMap>> mutationMap
    ) throws IllegalArgumentException {
        if (!node.isRoot() && !(parentGenotypes instanceof Integer[]) && !(parentGenotypes instanceof Double[]))
            throw new IllegalArgumentException(GENOTYPES_FORMAT_ERROR);

        Object genotypes = node.getMetaData(genotypeKey);
        if (node.isLeaf() && !(genotypes instanceof Integer[]) && !(genotypes instanceof Double[]))
            throw new IllegalArgumentException(GENOTYPES_FORMAT_ERROR);

        // Walk the tree in preorder
        for (Node child : node.getChildren())
            collectGenes4Tree(child, genotypeKey, genotypes, substModel, snvSites, mutationMap);

        if (!node.isRoot()) {
            for (int snv = 0; snv < snvSites.size(); snv++) {
                String geneName = getGeneName(snvSites.get(snv), mutationMap);

                // This snv site is located in the range of an important gene.
                if (geneName != null) {
                    final int childGenotype = (((Double[]) genotypes)[snv]).intValue();
                    final int parentGenotype = (((Double[]) parentGenotypes)[snv]).intValue();

                    if (parentGenotype != childGenotype) {
                        if (nodeGenesInfos[node.getNr()] == null)
                            nodeGenesInfos[node.getNr()] = new NodeGenesInfo();

                        ScsSubstitutionModelBase.EvolutionaryEventType[] events = substModel.getEvolutionaryEvents(
                                parentGenotype,
                                childGenotype
                        );

                        if (events != null && events.length > 0) {
                            // replace ';' with '/' in geneName to avoid subsequent parsing problems
                            geneName = geneName.replaceAll(";", "/");

                            nodeGenesInfos[node.getNr()].addGeneInfo(
                                    new GeneInfo(
                                            snvSites.get(snv),
                                            events,
                                            geneName
                                    )
                            );
                        }
                    }
                }
            }
        }
    } // annotateGenes2Tree

    private String getGeneName(final VariantSiteInfo info, final Map<ChromosomeLabel, List<MutationMap>> mutationMap) {
        final ChromosomeLabel cl = info.getChromosomeLabel();
        if (mutationMap.containsKey(cl)) {
            for (MutationMap map : mutationMap.get(cl)) {
                final long pos = info.getPosition();
                if (pos >= map.getStartPos() && pos <= map.getEndPos())
                    return map.getGeneName();
            }
        }

        return null;
    } // getGeneName

    private void processGenes() {
        List<GeneInfo> genes = new ArrayList<>();
        List<GeneInfo> genesISA = new ArrayList<>();
        List<GeneInfo> genesFSA = new ArrayList<>();

        for (NodeGenesInfo i : this.nodeGenesInfos) {
            if (i != null && i.getFullGenes() != null)
                genes.addAll(i.getFullGenes());
        }

        genes.sort(new GeneInfo.GeneInfoComparator());

        for (GeneInfo i : new HashSet<>(genes)) {
            final int j = Collections.frequency(genes, i);
            if (j > 1 || (j == 1 && i.violateISA()))
                genesFSA.add(i);
            else
                genesISA.add(i);
        }

        for (NodeGenesInfo i : this.nodeGenesInfos) {
            if (i != null && i.getFullGenes() != null) {
                for (GeneInfo j : i.getFullGenes()) {
                    if (genesISA.contains(j)) {
                        i.addISAGeneInfo(j);
                    } else if (genesFSA.contains(j)) {
                        i.addFSAGeneInfo(j);
                    }
                }

                i.convertGeneInfo();
            }
        }
    } // processGenes


    //**********************************************
    //*               Static methods               *
    //**********************************************

    private static String getOutputFileName(final String fileName) {
        return getRootPath(fileName) + getBaseName(fileName) + OUTPUT_SUFFIX;
    } // getOutputFileName

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("geneannotator", "");
    } // printUsage


    //***********************************************
    //*                 Main method                 *
    //***********************************************

    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // Variables
        int substModelLabel;
        String treeFileName;
        String snvFileName;
        String mutationMapFileName;
        String outputFileName;
        boolean mapFromAnnovar = false;
        int mapSeparator = -1;
        String filteringGenesFileName = null;
        int filteringGenesSeparator = 0;
        int filteringGenesColIndex = 0;

        if (args.length == 0) {

            Utils.loadUIManager();
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");
            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;
            if (url != null)
                icon = new javax.swing.ImageIcon(url);

            // Construct a new console
            new ConsoleApplication(null, null, icon, true);
            Log.info = System.out;
            Log.err = System.err;
            progressStream = System.out;

            // TODO: print some information here
            System.out.println("GeneAnnotator");

            GeneAnnotatorDialog dialog = new GeneAnnotatorDialog(new JFrame());
            if (!dialog.showDialog("GeneAnnotator"))
                return;

            // Get parameters
            treeFileName = dialog.getTreeFileName();
            if (treeFileName == null) {
                Log.err.println("No tree file specified!");
                return;
            }

            snvFileName = dialog.getSNVFileName();
            if (snvFileName == null) {
                Log.err.println("No snv file specified!");
                return;
            }

            substModelLabel = dialog.getSubstModelLabel();

            mutationMapFileName = dialog.getMutationMapFileName();
            if (mutationMapFileName == null) {
                Log.err.println("No mutation map file specified!");
                return;
            }

            mapFromAnnovar = dialog.isMapFromAnnovar();
            if (mapFromAnnovar)
                mapSeparator = dialog.getMapSeparator();

            filteringGenesFileName = dialog.getFilteringFileName();
            filteringGenesSeparator = dialog.getFilteringFileSeparator();
            filteringGenesColIndex = dialog.getColIndexFilteringGenes();

            outputFileName = dialog.getOutputFileName();
            if (outputFileName == null)
                outputFileName = getOutputFileName(treeFileName);

            try {
                new GeneAnnotator(
                        substModelLabel,
                        treeFileName,
                        snvFileName,
                        mutationMapFileName,
                        mapFromAnnovar,
                        mapSeparator != 0 && mapSeparator != 1 ? null : (mapSeparator == 0 ? "\t" : ","),
                        filteringGenesFileName,
                        filteringGenesSeparator != 0 && filteringGenesSeparator != 1 ? null : (filteringGenesSeparator == 0 ? "\t" : ","),
                        filteringGenesColIndex,
                        outputFileName
                );
            } catch (Exception e) {
                e.printStackTrace();
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
                            new Arguments.StringOption("tree", "input_tree_file", "specifies a tree file annotated with genotypes. This is usually generated by VariantCaller and ended with '.intermediate_tree'. -> MANDATORY"),
                            new Arguments.StringOption("snv", "snv_file", "specifies a file containing variant sites information. This is usually generated by VariantCaller and ended with '.loci_info'. -> MANDATORY"),
                            new Arguments.IntegerOption("subst", 0, 3, "specifies the substitution model used to infer phylogeny and call variants (0: ScsFiniteMuExtendedModel, 1: ScsFiniteMuModel, 2: ScsFiniteMuDelModel, 3: ScsFiniteMuDelInsModel). -> MANDATORY"),
                            new Arguments.StringOption("map", "mutation_map_file", "specifies the mutation map. -> MANDATORY"),
                            new Arguments.IntegerOption("anv", 0, 1, "specifies when the value of --map option is an output of Annovar (0: tab-separated, 1: comma-separated). -> OPTIONAL"),
                            new Arguments.StringOption("filter", "filtering_genes_file", "specifies filtering genes file (headers required). -> OPTIONAL"),
                            new Arguments.IntegerOption("sep", 0, 1, "specifies the separator of filtering genes file (0: tab-separated; 1: comma-separated; default: 0). -> OPTIONAL"),
                            new Arguments.IntegerOption("col", 0, 1000, "specifies the column index of gene names in the filtering genes file (starting from 0; default: 0). -> OPTIONAL"),
                            new Arguments.StringOption("out", "out_tree_file", "specifies the output tree file with genes annotated. -> OPTIONAL")
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

            // Set treeFileName
            if (arguments.hasOption("tree"))
                treeFileName = arguments.getStringOption("tree");
            else {
                Log.err.println("No tree file specified!");
                return;
            }

            // Set snvFileName
            if (arguments.hasOption("snv"))
                snvFileName = arguments.getStringOption("snv");
            else {
                Log.err.println("No snv file specified!");
                return;
            }

            // Set substModelLabel
            if (arguments.hasOption("subst"))
                substModelLabel = arguments.getIntegerOption("subst");
            else {
                Log.err.println("No substitution model specified! Check the usage.");
                return;
            }

            // Set mutationMapFileName
            if (arguments.hasOption("map"))
                mutationMapFileName = arguments.getStringOption("map");
            else {
                Log.err.println("No mutation map file specified!");
                return;
            }

            // Set mapFromAnnovar and mapSeparator
            if (arguments.hasOption("anv")) {
                mapFromAnnovar = true;
                mapSeparator = arguments.getIntegerOption("anv");
            }

            // Set filteringGenesFileName
            if (arguments.hasOption("filter"))
                filteringGenesFileName = arguments.getStringOption("filter");

            // Set filteringGenesSeparator
            if (arguments.hasOption("sep"))
                filteringGenesSeparator = arguments.getIntegerOption("sep");

            // Set filteringGenesColIndex
            if (arguments.hasOption("col"))
                filteringGenesColIndex = arguments.getIntegerOption("col");

            // Set outputFileName
            if (arguments.hasOption("out"))
                outputFileName = arguments.getStringOption("out");
            else
                outputFileName = getOutputFileName(treeFileName);

            try {
                new GeneAnnotator(
                        substModelLabel,
                        treeFileName,
                        snvFileName,
                        mutationMapFileName,
                        mapFromAnnovar,
                        mapSeparator != 0 && mapSeparator != 1 ? null : (mapSeparator == 0 ? "\t" : ","),
                        filteringGenesFileName,
                        filteringGenesSeparator != 0 && filteringGenesSeparator != 1 ? null : (filteringGenesSeparator == 0 ? "\t" : ","),
                        filteringGenesColIndex,
                        outputFileName
                );
            } catch (IOException e) {
                throw e;
            } catch (Exception e) {
                e.printStackTrace();
            }

            progressStream.println();
            progressStream.println("Successful!");
        }

    } // main

}
