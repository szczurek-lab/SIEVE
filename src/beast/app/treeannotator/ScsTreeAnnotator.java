package beast.app.treeannotator;

import beast.app.BEASTVersion;
import beast.app.beauti.BeautiDoc;
import beast.app.tools.LogCombiner;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeUtils;
import beast.math.statistic.DiscreteStatistics;
import beast.util.*;
import jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class ScsTreeAnnotator extends TreeAnnotator {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private final static BEASTVersion version = new BEASTVersion();

    private final static boolean USE_R = false;

    private static boolean forceIntegerToDiscrete = false;

    static boolean processSA = true;

    private boolean SAmode = false;

    private static boolean processBivariateAttributes = true;

    private static boolean saveSimpleTree = false;

    private static boolean hasTrunk = true;

    private final List<TreeAnnotationPlugin> beastObjects = new ArrayList<>();

    public enum Target {
        MAX_CLADE_CREDIBILITY("Maximum clade credibility tree"),
        MAX_SUM_CLADE_CREDIBILITY("Maximum sum of clade credibilities"),
        USER_TARGET_TREE("User target tree");

        String desc;

        Target(String s) {
            desc = s;
        }

        @Override
        public String toString() {
            return desc;
        }
    }

    public enum HeightsSummary {
        CA_HEIGHTS("Common Ancestor heights"),
        MEDIAN_HEIGHTS("Median heights"),
        MEAN_HEIGHTS("Mean heights"),
        KEEP_HEIGHTS("Keep target heights");

        String desc;

        HeightsSummary(String s) {
            desc = s;
        }

        @Override
        public String toString() {
            return desc;
        }
    }


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public ScsTreeAnnotator(
            final int burninPercentage,
            boolean lowMemory, // allowSingleChild was defunct (always set to false), now replaced by flag to say how much
            HeightsSummary heightsOption,
            double posteriorLimit,
            double hpd2D,
            Target targetOption,
            String targetTreeFileName,
            String inputFileName,
            String outputFileName,
            String simpleOutputFileName
    ) throws IOException {

        this.posteriorLimit = posteriorLimit;
        this.hpd2D = hpd2D;

        attributeNames.add("height");
        attributeNames.add("length");

        final long startTime = System.currentTimeMillis();

        CladeSystem cladeSystem = new ScsCladeSystem(hasTrunk);

        totalTrees = 10000;
        totalTreesUsed = 0;

        try {
            if (lowMemory) {
                treeSet = new MemoryFriendlyTreeSet(inputFileName, burninPercentage);
            } else {
                treeSet = new FastTreeSet(inputFileName, burninPercentage);
            }
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }

        if (targetOption != Target.USER_TARGET_TREE) {
            try {
                treeSet.reset();
                cladeSystem.setProcessSA(false);
                while (treeSet.hasNext()) {
                    Tree tree = treeSet.next();
                    tree.getLeafNodeCount();
                    if (tree.getDirectAncestorNodeCount() > 0 && !SAmode && processSA) {
                        SAmode = true;
                        Log.err.println("A tree with a sampled ancestor is found. Turning on\n the sampled ancestor " +
                                "summary analysis.");
                        if (heightsOption == HeightsSummary.CA_HEIGHTS) {
                            throw new RuntimeException("The common ancestor height is not \n available for trees with sampled " +
                                    "ancestors. Please choose \n another height summary option");
                        }
                        cladeSystem.setProcessSA(true);
                    }
                    cladeSystem.add(tree, false);
                    totalTreesUsed++;
                }
                totalTrees = totalTreesUsed * 100 / (100 - Math.max(burninPercentage, 0));
            } catch (Exception e) {
                Log.err.println(e.getMessage());
                return;
            }

            progressStream.println();
            progressStream.println();

            if (totalTrees < 1) {
                Log.err.println("No trees");
                return;
            }
            if (totalTreesUsed <= 1) {
                if (burninPercentage > 0) {
                    Log.err.println("No trees to use: burnin too high");
                    return;
                }
            }
            cladeSystem.calculateCladeCredibilities(totalTreesUsed);

            progressStream.println("Total number of trees " + totalTrees + ", where " + totalTreesUsed + " are used.");

            progressStream.println("Total unique clades: " + cladeSystem.getCladeMap().keySet().size());
            progressStream.println();
        } else {
            // even when a user specified target tree is provided we still need to count the totalTreesUsed for subsequent steps.
            treeSet.reset();
            while (treeSet.hasNext()) {
                Tree tree = treeSet.next();
                tree.getLeafNodeCount();
                if (tree.getDirectAncestorNodeCount() > 0 && !SAmode && processSA) {
                    SAmode = true;
                    Log.err.println("A tree with a sampled ancestor is found. Turning on\n the sampled ancestor " +
                            "summary analysis.");
                    if (heightsOption == HeightsSummary.CA_HEIGHTS) {
                        throw new RuntimeException("The common ancestor height is not \n available for trees with sampled " +
                                "ancestors. Please choose \n another height summary option");
                    }
                }
                totalTreesUsed++;
            }
        }

        Tree targetTree = null;

        switch (targetOption) {
            case USER_TARGET_TREE: {
                if (targetTreeFileName != null) {
                    progressStream.println("Reading user specified target tree, " + targetTreeFileName);

                    String tree = BeautiDoc.load(targetTreeFileName);

                    if (tree.startsWith("#NEXUS")) {
                        NexusParser parser2 = new NexusParser();
                        parser2.parseFile(new File(targetTreeFileName));
                        targetTree = parser2.trees.get(0);
                    } else {
                        try {
                            TreeParser parser2 = new TreeParser();
                            parser2.initByName("IsLabelledNewick", true, "newick", tree);
                            targetTree = parser2;
                        } catch (Exception e) {
                            Log.err.println("Error Parsing Target Tree: " + e.getMessage());
                            return;
                        }
                    }
                } else {
                    Log.err.println("No user target tree specified.");
                    return;
                }
                break;
            }
            case MAX_CLADE_CREDIBILITY: {
                progressStream.println("Finding maximum credibility tree...");
                targetTree = summarizeTrees(cladeSystem, false).copy();
                break;
            }
            case MAX_SUM_CLADE_CREDIBILITY: {
                progressStream.println("Finding maximum sum clade credibility tree...");
                targetTree = summarizeTrees(cladeSystem, true).copy();
                break;
            }
        }

        // statistics
        progressStream.println("Collecting node information...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = Math.max(totalTreesUsed / 60, 1);
        int reported = 0;

        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        cladeSystem = new ScsCladeSystem(hasTrunk);
        cladeSystem.setProcessSA(processSA);
        cladeSystem.add(targetTree, true);
        int totalTreesUsedNew = 0;
        try {
            int counter = 0;
            treeSet.reset();
            while (treeSet.hasNext()) {
                Tree tree = treeSet.next();
                if (counter == 0) {
                    setupAttributes(tree);
                }
                cladeSystem.collectAttributes(tree, attributeNames);
                if (counter > 0 && counter % stepSize == 0 && reported < 61) {
                    while (1000 * reported < 61000 * (counter + 1) / this.totalTreesUsed) {
                        progressStream.print("*");
                        reported++;
                    }
                    progressStream.flush();
                }
                totalTreesUsedNew++;
                counter++;
            }

            cladeSystem.removeClades(targetTree.getRoot(), true);
            this.totalTreesUsed = totalTreesUsedNew;
            cladeSystem.calculateCladeCredibilities(totalTreesUsedNew);
        } catch (Exception e) {
            Log.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        progressStream.println();
        progressStream.println();

        progressStream.println("Annotating target tree...");

        try {
            annotateTree(cladeSystem, targetTree.getRoot(), null, heightsOption);

            if (heightsOption == HeightsSummary.CA_HEIGHTS) {
                setTreeHeightsByCA(targetTree, targetOption);
            }
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error to annotate tree: " + e.getMessage() + "\nPlease check the tree log file format.");
            return;
        }

        progressStream.println("Writing annotated tree....");

        processMetaData(targetTree.getRoot());
        try {

            File outputTree = null;

            if (outputFileName != null) {
                if (System.getProperty("tree.annotation.file.prefix") != null) {
                    outputFileName = System.getProperty("tree.annotation.file.prefix") + outputFileName;
                }

                outputTree = new File(outputFileName);
                if (outputFileName.contains("/")) {
                    outputTree.getParentFile().mkdirs();
                }
            }

            final PrintStream stream = outputTree != null ?
                    new PrintStream(outputTree) :
                    System.out;

            targetTree.init(stream);
            stream.println();

            stream.print("tree TREE1 = ");
            int[] dummy = new int[1];
            String newick = targetTree.getRoot().toSortedNewick(dummy, true);
            stream.print(newick);
            stream.println(";");
//            stream.println(targetTree.getRoot().toShortNewick(false));
//            stream.println();
            targetTree.close(stream);
            stream.println();
        } catch (Exception e) {
            Log.err.println("Error to write annotated tree file: " + e.getMessage());
        }

        // save a simpler, no meta data tree if necessary
        if (saveSimpleTree) {
            try {

                File outputTree = null;

                if (simpleOutputFileName != null) {
                    if (System.getProperty("tree.annotation.file.prefix") != null) {
                        simpleOutputFileName = System.getProperty("tree.annotation.file.prefix") + simpleOutputFileName;
                    }

                    outputTree = new File(simpleOutputFileName);
                    if (simpleOutputFileName.contains("/")) {
                        outputTree.getParentFile().mkdirs();
                    }
                }

                final PrintStream stream = outputTree != null ?
                        new PrintStream(outputTree) :
                        System.out;

                targetTree.init(stream);
                stream.println();

                stream.print("tree TREE1 = ");
                int[] dummy = new int[1];
                String newick = targetTree.getRoot().toSortedNewick(dummy, false);
                stream.print(newick);
                stream.println(";");
//            stream.println(targetTree.getRoot().toShortNewick(false));
//            stream.println();
                targetTree.close(stream);
                stream.println();
            } catch (Exception e) {
                Log.err.println("Error to write annotated tree file: " + e.getMessage());
            }
        }

        final long endTime = System.currentTimeMillis();
        progressStream.println();
        progressStream.println(">>> Tree annotation completed in " + (endTime - startTime) / 1000 + "s.");

    }

    public ScsTreeAnnotator(
            boolean processBivariateAttributes,
            boolean forceIntegerToDiscrete,
            boolean processSA,
            boolean SAmode,
            boolean saveSimpleTree,
            boolean hasTrunk,
            final int burninPercentage,
            boolean lowMemory,
            HeightsSummary heightsOption,
            double posteriorLimit,
            double hpd2D,
            Target targetOption,
            String targetTreeFileName,
            String inputFileName,
            String outputFileName,
            String simpleOutputFileName
    ) throws IOException {
        ScsTreeAnnotator.forceIntegerToDiscrete = forceIntegerToDiscrete;
        ScsTreeAnnotator.processSA = processSA;
        this.SAmode = SAmode;
        ScsTreeAnnotator.saveSimpleTree = saveSimpleTree;
        ScsTreeAnnotator.hasTrunk = hasTrunk;

        this.posteriorLimit = posteriorLimit;
        this.hpd2D = hpd2D;

        attributeNames.add("height");
        attributeNames.add("length");

        final long startTime = System.currentTimeMillis();

        CladeSystem cladeSystem = new ScsCladeSystem(hasTrunk);

        totalTrees = 10000;
        totalTreesUsed = 0;

        try {
            if (lowMemory) {
                treeSet = new MemoryFriendlyTreeSet(inputFileName, burninPercentage);
            } else {
                treeSet = new FastTreeSet(inputFileName, burninPercentage);
            }
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }

        if (targetOption != Target.USER_TARGET_TREE) {
            try {
                treeSet.reset();
                cladeSystem.setProcessSA(false);
                while (treeSet.hasNext()) {
                    Tree tree = treeSet.next();
                    tree.getLeafNodeCount();
                    if (tree.getDirectAncestorNodeCount() > 0 && !SAmode && processSA) {
                        SAmode = true;
                        Log.err.println("A tree with a sampled ancestor is found. Turning on\n the sampled ancestor " +
                                "summary analysis.");
                        if (heightsOption == HeightsSummary.CA_HEIGHTS) {
                            throw new RuntimeException("The common ancestor height is not \n available for trees with sampled " +
                                    "ancestors. Please choose \n another height summary option");
                        }
                        cladeSystem.setProcessSA(true);
                    }
                    cladeSystem.add(tree, false);
                    totalTreesUsed++;
                }
                totalTrees = totalTreesUsed * 100 / (100 - Math.max(burninPercentage, 0));
            } catch (Exception e) {
                Log.err.println(e.getMessage());
                return;
            }

            progressStream.println();
            progressStream.println();

            if (totalTrees < 1) {
                Log.err.println("No trees");
                return;
            }
            if (totalTreesUsed <= 1) {
                if (burninPercentage > 0) {
                    Log.err.println("No trees to use: burnin too high");
                    return;
                }
            }
            cladeSystem.calculateCladeCredibilities(totalTreesUsed);

            progressStream.println("Total number of trees " + totalTrees + ", where " + totalTreesUsed + " are used.");

            progressStream.println("Total unique clades: " + cladeSystem.getCladeMap().keySet().size());
            progressStream.println();
        } else {
            // even when a user specified target tree is provided we still need to count the totalTreesUsed for subsequent steps.
            treeSet.reset();
            while (treeSet.hasNext()) {
                Tree tree = treeSet.next();
                tree.getLeafNodeCount();
                if (tree.getDirectAncestorNodeCount() > 0 && !SAmode && processSA) {
                    SAmode = true;
                    Log.err.println("A tree with a sampled ancestor is found. Turning on\n the sampled ancestor " +
                            "summary analysis.");
                    if (heightsOption == HeightsSummary.CA_HEIGHTS) {
                        throw new RuntimeException("The common ancestor height is not \n available for trees with sampled " +
                                "ancestors. Please choose \n another height summary option");
                    }
                }
                totalTreesUsed++;
            }
        }

        Tree targetTree = null;

        switch (targetOption) {
            case USER_TARGET_TREE: {
                if (targetTreeFileName != null) {
                    progressStream.println("Reading user specified target tree, " + targetTreeFileName);

                    String tree = BeautiDoc.load(targetTreeFileName);

                    if (tree.startsWith("#NEXUS")) {
                        NexusParser parser2 = new NexusParser();
                        parser2.parseFile(new File(targetTreeFileName));
                        targetTree = parser2.trees.get(0);
                    } else {
                        try {
                            TreeParser parser2 = new TreeParser();
                            parser2.initByName("IsLabelledNewick", true, "newick", tree);
                            targetTree = parser2;
                        } catch (Exception e) {
                            Log.err.println("Error Parsing Target Tree: " + e.getMessage());
                            return;
                        }
                    }
                } else {
                    Log.err.println("No user target tree specified.");
                    return;
                }
                break;
            }
            case MAX_CLADE_CREDIBILITY: {
                progressStream.println("Finding maximum credibility tree...");
                targetTree = summarizeTrees(cladeSystem, false).copy();
                break;
            }
            case MAX_SUM_CLADE_CREDIBILITY: {
                progressStream.println("Finding maximum sum clade credibility tree...");
                targetTree = summarizeTrees(cladeSystem, true).copy();
                break;
            }
        }

        // statistics
        progressStream.println("Collecting node information...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = Math.max(totalTreesUsed / 60, 1);
        int reported = 0;

        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        cladeSystem = new ScsCladeSystem(hasTrunk);
        cladeSystem.setProcessSA(processSA);
        cladeSystem.add(targetTree, true);
        int totalTreesUsedNew = 0;
        try {
            int counter = 0;
            treeSet.reset();
            while (treeSet.hasNext()) {
                Tree tree = treeSet.next();
                if (counter == 0) {
                    setupAttributes(tree);
                }
                cladeSystem.collectAttributes(tree, attributeNames);
                if (counter > 0 && counter % stepSize == 0 && reported < 61) {
                    while (1000 * reported < 61000 * (counter + 1) / this.totalTreesUsed) {
                        progressStream.print("*");
                        reported++;
                    }
                    progressStream.flush();
                }
                totalTreesUsedNew++;
                counter++;
            }

            cladeSystem.removeClades(targetTree.getRoot(), true);
            this.totalTreesUsed = totalTreesUsedNew;
            cladeSystem.calculateCladeCredibilities(totalTreesUsedNew);
        } catch (Exception e) {
            Log.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        progressStream.println();
        progressStream.println();

        progressStream.println("Annotating target tree...");

        try {
            annotateTree(cladeSystem, targetTree.getRoot(), null, heightsOption);

            if (heightsOption == HeightsSummary.CA_HEIGHTS) {
                setTreeHeightsByCA(targetTree, targetOption);
            }
        } catch (Exception e) {
            e.printStackTrace();
            Log.err.println("Error to annotate tree: " + e.getMessage() + "\nPlease check the tree log file format.");
            return;
        }

        progressStream.println("Writing annotated tree....");

        processMetaData(targetTree.getRoot());
        try {

            File outputTree = null;

            if (outputFileName != null) {
                if (System.getProperty("tree.annotation.file.prefix") != null) {
                    outputFileName = System.getProperty("tree.annotation.file.prefix") + outputFileName;
                }

                outputTree = new File(outputFileName);
                if (outputFileName.contains("/")) {
                    outputTree.getParentFile().mkdirs();
                }
            }

            final PrintStream stream = outputTree != null ?
                    new PrintStream(outputTree) :
                    System.out;

            targetTree.init(stream);
            stream.println();

            stream.print("tree TREE1 = ");
            int[] dummy = new int[1];
            String newick = targetTree.getRoot().toSortedNewick(dummy, true);
            stream.print(newick);
            stream.println(";");
//            stream.println(targetTree.getRoot().toShortNewick(false));
//            stream.println();
            targetTree.close(stream);
            stream.println();
        } catch (Exception e) {
            Log.err.println("Error to write annotated tree file: " + e.getMessage());
        }

        // save a simpler, no meta data tree if necessary
        if (saveSimpleTree) {
            try {

                File outputTree = null;

                if (simpleOutputFileName != null) {
                    if (System.getProperty("tree.annotation.file.prefix") != null) {
                        simpleOutputFileName = System.getProperty("tree.annotation.file.prefix") + simpleOutputFileName;
                    }

                    outputTree = new File(simpleOutputFileName);
                    if (simpleOutputFileName.contains("/")) {
                        outputTree.getParentFile().mkdirs();
                    }
                }

                final PrintStream stream = outputTree != null ?
                        new PrintStream(outputTree) :
                        System.out;

                targetTree.init(stream);
                stream.println();

                stream.print("tree TREE1 = ");
                int[] dummy = new int[1];
                String newick = targetTree.getRoot().toSortedNewick(dummy, false);
                stream.print(newick);
                stream.println(";");
//            stream.println(targetTree.getRoot().toShortNewick(false));
//            stream.println();
                targetTree.close(stream);
                stream.println();
            } catch (Exception e) {
                Log.err.println("Error to write annotated tree file: " + e.getMessage());
            }
        }

        final long endTime = System.currentTimeMillis();
        progressStream.println();
        progressStream.println(">>> Tree annotation completed in " + (endTime - startTime) / 1000 + "s.");

    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    private Tree summarizeTrees(CladeSystem cladeSystem, boolean useSumCladeCredibility) throws IOException {

        Tree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        progressStream.println("Analyzing " + totalTreesUsed + " trees...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = Math.max(totalTreesUsed / 60, 1);
        int reported = 0;

        int counter = 0;
        treeSet.reset();
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            double score = scoreTree(tree, cladeSystem, useSumCladeCredibility);
            if (score > bestScore) {
                bestTree = tree;
                bestScore = score;
            }
            if (counter % stepSize == 0 && reported < 61) {
                while (1000 * reported < 61000 * (counter + 1) / totalTreesUsed) {
                    progressStream.print("*");
                    reported++;
                }
                progressStream.flush();
            }
            counter++;
        }
        progressStream.println();
        progressStream.println();
        if (useSumCladeCredibility) {
            progressStream.println("Highest Sum Clade Credibility: " + bestScore);
        } else {
            progressStream.println("Highest Log Clade Credibility: " + bestScore);
        }

        return bestTree;
    } // summarizeTrees

    public double scoreTree(Tree tree, CladeSystem cladeSystem, boolean useSumCladeCredibility) {
        if (useSumCladeCredibility) {
            return cladeSystem.getSumCladeCredibility(tree.getRoot(), null);
        } else {
            return cladeSystem.getLogCladeCredibility(tree.getRoot(), null);
        }
    } // scoreTree

    private void setupAttributes(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
            Set<String> iter = node.getMetaDataNames();
            if (iter != null) {
                attributeNames.addAll(iter);
            }
        }

        for (TreeAnnotationPlugin beastObject : beastObjects) {
            Set<String> claimed = beastObject.setAttributeNames(attributeNames);
            attributeNames.removeAll(claimed);
        }
    } // setupAttributes

    private void annotateTree(CladeSystem cladeSystem, Node node, BitSet bits, HeightsSummary heightsOption) {

        BitSet bits2 = new BitSet();

        if (node.isLeaf()) {

            int index = cladeSystem.getTaxonIndex(node);
            bits2.set(2 * index);

            annotateNode(cladeSystem, node, bits2, true, heightsOption);
        } else {

            for (int i = 0; i < node.getChildCount(); i++) {

                Node node1 = node.getChild(i);

                annotateTree(cladeSystem, node1, bits2, heightsOption);
            }

            for (int i = 1; i < bits2.length(); i = i + 2) {
                bits2.set(i, false);
            }
            if (node.isFake() && processSA) {
                int index = cladeSystem.getTaxonIndex(node.getDirectAncestorChild());
                bits2.set(2 * index + 1);
            }

            annotateNode(cladeSystem, node, bits2, false, heightsOption);
        }

        if (bits != null) {
            bits.or(bits2);
        }
    } // annotateTree

    private void annotateNode(CladeSystem cladeSystem, Node node, BitSet bits, boolean isTip, HeightsSummary heightsOption) {
        final boolean isTrunkRoot = node.isRoot() && hasTrunk;

        CladeSystem.Clade clade = cladeSystem.cladeMap.get(bits);
        assert clade != null : "Clade missing?";

        final List<Object[]> attributeValues = isTrunkRoot ? ((ScsCladeSystem) cladeSystem).getRootAttributeValues() : clade.attributeValues;

        boolean filter = false;
        if (!isTip) {
            if (isTrunkRoot) {
                node.setMetaData("posterior", 1.0);
            } else {
                final double posterior = clade.getCredibility();
                node.setMetaData("posterior", posterior);
                if (posterior < posteriorLimit) {
                    filter = true;
                }
            }
        }

        int i = 0;
        for (String attributeName : attributeNames) {

            if (attributeValues != null && attributeValues.size() > 0) {
                double[] values = new double[attributeValues.size()];

                HashMap<Object, Integer> hashMap = new HashMap<>();

                Object[] v = attributeValues.get(0);
                if (v[i] != null) {

                    final boolean isHeight = attributeName.equals("height");
                    boolean isBoolean = v[i] instanceof Boolean;

                    boolean isDiscrete = v[i] instanceof String;

                    if (forceIntegerToDiscrete && v[i] instanceof Integer) isDiscrete = true;

                    double minValue = Double.MAX_VALUE;
                    double maxValue = -Double.MAX_VALUE;

                    final boolean isArray = v[i] instanceof Object[];
                    boolean isDoubleArray = isArray && ((Object[]) v[i])[0] instanceof Double;
                    // This is Java, friends - first value type does not imply all.
                    if (isDoubleArray) {
                        for (Object n : (Object[]) v[i]) {
                            if (!(n instanceof Double)) {
                                isDoubleArray = false;
                                break;
                            }
                        }
                    }
                    // todo Handle other types of arrays

                    double[][] valuesArray = null;
                    double[] minValueArray = null;
                    double[] maxValueArray = null;
                    int lenArray = 0;

                    if (isDoubleArray) {
                        lenArray = ((Object[]) v[i]).length;

                        valuesArray = new double[lenArray][attributeValues.size()];
                        minValueArray = new double[lenArray];
                        maxValueArray = new double[lenArray];

                        for (int k = 0; k < lenArray; k++) {
                            minValueArray[k] = Double.MAX_VALUE;
                            maxValueArray[k] = -Double.MAX_VALUE;
                        }
                    }

                    for (int j = 0; j < attributeValues.size(); j++) {
                        Object value = attributeValues.get(j)[i];
                        if (isDiscrete) {
                            final Object s = value;
                            if (hashMap.containsKey(s)) {
                                hashMap.put(s, hashMap.get(s) + 1);
                            } else {
                                hashMap.put(s, 1);
                            }
                        } else if (isBoolean) {
                            values[j] = (((Boolean) value) ? 1.0 : 0.0);
                        } else if (isDoubleArray) {
                            // Forcing to Double[] causes a cast exception. MAS
                            try {
                                Object[] array = (Object[]) value;
                                for (int k = 0; k < lenArray; k++) {
                                    valuesArray[k][j] = ((Double) array[k]);
                                    if (valuesArray[k][j] < minValueArray[k]) minValueArray[k] = valuesArray[k][j];
                                    if (valuesArray[k][j] > maxValueArray[k]) maxValueArray[k] = valuesArray[k][j];
                                }
                            } catch (Exception e) {
                                // ignore
                            }
                        } else {
                            // Ignore other (unknown) types
                            if (value instanceof Number) {
                                values[j] = ((Number) value).doubleValue();
                                if (values[j] < minValue) minValue = values[j];
                                if (values[j] > maxValue) maxValue = values[j];
                            }
                        }
                    }
                    if (isHeight) {
                        if (heightsOption == HeightsSummary.MEAN_HEIGHTS) {
                            final double mean = DiscreteStatistics.mean(values);
                            if (node.isDirectAncestor()) {
                                node.getParent().setHeight(mean);
                            }
                            if (node.isFake() && processSA) {
                                node.getDirectAncestorChild().setHeight(mean);
                            }
                            node.setHeight(mean);
                        } else if (heightsOption == HeightsSummary.MEDIAN_HEIGHTS) {
                            final double median = DiscreteStatistics.median(values);
                            if (node.isDirectAncestor()) {
                                node.getParent().setHeight(median);
                            }
                            if (node.isFake() && processSA) {
                                node.getDirectAncestorChild().setHeight(median);
                            }
                            node.setHeight(median);
                        } else {
                            // keep the existing height
                        }
                    }

                    if (!filter) {
                        boolean processed = false;
                        for (TreeAnnotationPlugin beastObject : beastObjects) {
                            if (beastObject.handleAttribute(node, attributeName, values)) {
                                processed = true;
                            }
                        }

                        if (!processed) {
                            if (!isDiscrete) {
                                if (!isDoubleArray)
                                    annotateMeanAttribute(node, attributeName, values);
                                else {
                                    for (int k = 0; k < lenArray; k++) {
                                        annotateMeanAttribute(node, attributeName + (k + 1), valuesArray[k]);
                                    }
                                }
                            } else {
                                annotateModeAttribute(node, attributeName, hashMap);
                                annotateFrequencyAttribute(node, attributeName, hashMap);
                            }
                            if (!isBoolean && minValue < maxValue && !isDiscrete && !isDoubleArray) {
                                // Basically, if it is a boolean (0, 1) then we don't need the distribution information
                                // Likewise if it doesn't vary.
                                annotateMedianAttribute(node, attributeName + "_median", values);
                                annotateHPDAttribute(node, attributeName + "_95%_HPD", 0.95, values);
                                annotateRangeAttribute(node, attributeName + "_range", values);
                            }

                            if (isDoubleArray) {
                                String name = attributeName;
                                // todo
//                                    if (name.equals(location1Attribute)) {
//                                        name = locationOutputAttribute;
//                                    }
                                boolean want2d = processBivariateAttributes && lenArray == 2;
                                if (name.equals("dmv")) {  // terrible hack
                                    want2d = false;
                                }
                                for (int k = 0; k < lenArray; k++) {
                                    if (minValueArray[k] < maxValueArray[k]) {
                                        annotateMedianAttribute(node, name + (k + 1) + "_median", valuesArray[k]);
                                        annotateRangeAttribute(node, name + (k + 1) + "_range", valuesArray[k]);
                                        if (!want2d)
                                            annotateHPDAttribute(node, name + (k + 1) + "_95%_HPD", 0.95, valuesArray[k]);
                                    }
                                }
                                // 2D contours
                                if (want2d) {

                                    boolean variationInFirst = (minValueArray[0] < maxValueArray[0]);
                                    boolean variationInSecond = (minValueArray[1] < maxValueArray[1]);

                                    if (variationInFirst && !variationInSecond)
                                        annotateHPDAttribute(node, name + "1" + "_95%_HPD", 0.95, valuesArray[0]);

                                    if (variationInSecond && !variationInFirst)
                                        annotateHPDAttribute(node, name + "2" + "_95%_HPD", 0.95, valuesArray[1]);

                                    if (variationInFirst && variationInSecond)
                                        annotate2DHPDAttribute(node, name, "_" + (int) (100 * hpd2D) + "%HPD", hpd2D, valuesArray);
                                }
                            }
                        }
                    }
                }
            }
            i++;
        }
    } // annotateNode

    private void annotateMeanAttribute(Node node, String label, double[] values) {
        double mean = DiscreteStatistics.mean(values);
        node.setMetaData(label, mean);
    } // annotateMeanAttribute

    private void annotateMedianAttribute(Node node, String label, double[] values) {
        double median = DiscreteStatistics.median(values);
        node.setMetaData(label, median);
    } // annotateMedianAttribute

    private void annotateModeAttribute(Node node, String label, HashMap<Object, Integer> values) {
        Object mode = null;
        int maxCount = 0;
        int totalCount = 0;
        int countInMode = 1;

        for (Object key : values.keySet()) {
            int thisCount = values.get(key);
            if (thisCount == maxCount) {
                // I hope this is the intention
                mode = mode.toString().concat("+" + key);
                countInMode++;
            } else if (thisCount > maxCount) {
                mode = key;
                maxCount = thisCount;
                countInMode = 1;
            }
            totalCount += thisCount;
        }
        double freq = (double) maxCount / (double) totalCount * countInMode;
        node.setMetaData(label, mode);
        node.setMetaData(label + ".prob", freq);
    } // annotateModeAttribute

    private void annotateFrequencyAttribute(Node node, String label, HashMap<Object, Integer> values) {
        double totalCount = 0;
        Set<?> keySet = values.keySet();
        int length = keySet.size();
        String[] name = new String[length];
        Double[] freq = new Double[length];
        int index = 0;
        for (Object key : values.keySet()) {
            name[index] = key.toString();
            freq[index] = Double.valueOf(values.get(key));
            totalCount += freq[index];
            index++;
        }
        for (int i = 0; i < length; i++)
            freq[i] /= totalCount;

        node.setMetaData(label + ".set", name);
        node.setMetaData(label + ".set.prob", freq);
    } // annotateFrequencyAttribute

    private void annotateRangeAttribute(Node node, String label, double[] values) {
        double min = DiscreteStatistics.min(values);
        double max = DiscreteStatistics.max(values);
        node.setMetaData(label, new Object[]{min, max});
    } // annotateRangeAttribute

    private void annotateHPDAttribute(Node node, String label, double hpd, double[] values) {
        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(hpd * values.length);
        for (int i = 0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex + diff - 1]];
        node.setMetaData(label, new Object[]{lower, upper});
    } // annotateHPDAttribute

    private String formattedLocation(double x) {
        return String.format("%5.2f", x);
    } // formattedLocation

    private void annotate2DHPDAttribute(Node node, String preLabel, String postLabel, double hpd, double[][] values) {
        if (USE_R) {

            // Uses R-Java interface, and the HPD routines from 'emdbook' and 'coda'

//                int N = 50;
//                if (rEngine == null) {
//
//                    if (!Rengine.versionCheck()) {
//                        throw new RuntimeException("JRI library version mismatch");
//                    }
//
//                    rEngine = new Rengine(rArgs, false, null);
//
//                    if (!rEngine.waitForR()) {
//                        throw new RuntimeException("Cannot load R");
//                    }
//
//                    for (String command : rBootCommands) {
//                        rEngine.eval(command);
//                    }
//                }
//
//                // todo Need a good method to pick grid size
//
//
//                REXP x = rEngine.eval("makeContour(" +
//                        makeRString(values[0]) + "," +
//                        makeRString(values[1]) + "," +
//                        hpd + "," +
//                        N + ")");
//
//                RVector contourList = x.asVector();
//                int numberContours = contourList.size();
//
//                if (numberContours > 1) {
//                    Log.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
//                    Log.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
//                }
//
//
//                node.setMetaData(preLabel + postLabel + "_modality", numberContours);
//
//                StringBuffer output = new StringBuffer();
//                for (int i = 0; i < numberContours; i++) {
//                    output.append("\n<" + CORDINATE + ">\n");
//                    RVector oneContour = contourList.at(i).asVector();
//                    double[] xList = oneContour.at(1).asDoubleArray();
//                    double[] yList = oneContour.at(2).asDoubleArray();
//                    StringBuffer xString = new StringBuffer("{");
//                    StringBuffer yString = new StringBuffer("{");
//                    for (int k = 0; k < xList.length; k++) {
//                        xString.append(formattedLocation(xList[k])).append(",");
//                        yString.append(formattedLocation(yList[k])).append(",");
//                    }
//                    xString.append(formattedLocation(xList[0])).append("}");
//                    yString.append(formattedLocation(yList[0])).append("}");
//
//                    node.setMetaData(preLabel + "1" + postLabel + "_" + (i + 1), xString);
//                    node.setMetaData(preLabel + "2" + postLabel + "_" + (i + 1), yString);
//                }


        } else { // do not use R


//                KernelDensityEstimator2D kde = new KernelDensityEstimator2D(values[0], values[1], N);
            //ContourMaker kde = new ContourWithSynder(values[0], values[1], N);
            boolean bandwidthLimit = false;

            ContourMaker kde = new ContourWithSynder(values[0], values[1], bandwidthLimit);

            ContourPath[] paths = kde.getContourPaths(hpd);

            node.setMetaData(preLabel + postLabel + "_modality", paths.length);

            if (paths.length > 1) {
                Log.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
                Log.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
            }

            StringBuffer output = new StringBuffer();
            int i = 0;
            for (ContourPath p : paths) {
                output.append("\n<" + CORDINATE + ">\n");
                double[] xList = p.getAllX();
                double[] yList = p.getAllY();
                StringBuffer xString = new StringBuffer("{");
                StringBuffer yString = new StringBuffer("{");
                for (int k = 0; k < xList.length; k++) {
                    xString.append(formattedLocation(xList[k])).append(",");
                    yString.append(formattedLocation(yList[k])).append(",");
                }
                xString.append(formattedLocation(xList[0])).append("}");
                yString.append(formattedLocation(yList[0])).append("}");

                node.setMetaData(preLabel + "1" + postLabel + "_" + (i + 1), xString);
                node.setMetaData(preLabel + "2" + postLabel + "_" + (i + 1), yString);
                i++;

            }
        }
    } // annotate2DHPDAttribute

    boolean setTreeHeightsByCA(Tree targetTree, Target targetOption) throws IOException {
        progressStream.println("Setting node heights...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int reportStepSize = totalTreesUsed / 60;
        if (reportStepSize < 1) reportStepSize = 1;
        int reported = 0;

        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        CladeSystem cladeSystem = new ScsCladeSystem(hasTrunk, targetTree);
        final int clades = cladeSystem.getCladeMap().size();

        // allocate posterior tree nodes order once
        int[] postOrderList;
        BitSet[] ctarget = new BitSet[clades];
        BitSet[] ctree = new BitSet[clades];

        for (int k = 0; k < clades; ++k) {
            ctarget[k] = new BitSet();
            ctree[k] = new BitSet();
        }

        cladeSystem.getTreeCladeCodes(targetTree, ctarget);

        // temp collecting heights inside loop allocated once
        double[][] hs;

        // heights total sum from posterior trees
        double[] ths;

        if (hasTrunk) {
            postOrderList = new int[clades + 1];
            hs = new double[clades + 1][treeSet.totalTrees];
            ths = new double[clades + 1];
        } else {
            postOrderList = new int[clades];
            hs = new double[clades][treeSet.totalTrees];
            ths = new double[clades];
        }

        int totalTreesUsed = 0;

        int counter = 0;
        treeSet.reset();
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            TreeUtils.preOrderTraversalList(tree, postOrderList);
            cladeSystem.getTreeCladeCodes(tree, ctree);
            for (int k = 0; k < clades; ++k) {
                int j = postOrderList[hasTrunk ? k + 1 : k];
                for (int i = 0; i < clades; ++i) {
                    if (CollectionUtils.isSubSet(ctarget[i], ctree[j])) {
                        hs[i][counter] = tree.getNode(j).getHeight();
                    }
                }
            }
            if (hasTrunk) {
                hs[clades][counter] = tree.getNode(postOrderList[0]).getHeight();
            }

            for (int k = 0; k < (hasTrunk ? clades + 1 : clades); ++k) {
                ths[k] += hs[k][counter];
            }
            totalTreesUsed += 1;
            if (counter > 0 && counter % reportStepSize == 0 && reported < 61) {
                while (1000 * reported < 61000 * (counter + 1) / this.totalTreesUsed) {
                    progressStream.print("*");
                    reported++;
                }
                progressStream.flush();
            }
            counter++;

        }

        if (targetOption != Target.USER_TARGET_TREE)
            targetTree.initAndValidate();

        cladeSystem.removeClades(targetTree.getRoot(), true);
        for (int k = 0; k < (hasTrunk ? clades + 1 : clades); ++k) {
            ths[k] /= totalTreesUsed;
            final Node node = targetTree.getNode(k);
            node.setHeight(ths[k]);
            String attributeName = "CAheight";
            double[] values = hs[k];
            double min = values[0];
            double max = values[0];
            for (double d : values) {
                min = Math.min(d, min);
                max = Math.max(d, max);
            }
            if (Math.abs(min - max) > 1e-10) {
                annotateMeanAttribute(node, attributeName + "_mean", values);
                annotateMedianAttribute(node, attributeName + "_median", values);
                annotateHPDAttribute(node, attributeName + "_95%_HPD", 0.95, values);
                annotateRangeAttribute(node, attributeName + "_range", values);
            }
        }

        assert (totalTreesUsed == this.totalTreesUsed);
        this.totalTreesUsed = totalTreesUsed;
        progressStream.println();
        progressStream.println();

        return true;
    } // setTreeHeightsByCA

    private void processMetaData(Node node) {
        for (Node child : node.getChildren()) {
            processMetaData(child);
        }
        Set<String> metaDataNames = node.getMetaDataNames();
        if (metaDataNames != null && !metaDataNames.isEmpty()) {
            String metadata = "";
            for (String name : metaDataNames) {
                Object value = node.getMetaData(name);
                metadata += name + "=";
                if (value instanceof Object[]) {
                    Object[] values = (Object[]) value;
                    metadata += "{";
                    for (int i = 0; i < values.length; i++) {
                        metadata += values[i].toString();
                        if (i < values.length - 1) {
                            metadata += ",";
                        }
                    }
                    metadata += "}";
                } else {
                    metadata += value.toString();
                }
                metadata += ",";
            }
            metadata = metadata.substring(0, metadata.length() - 1);
            node.metaDataString = metadata;
        }
    } // processMetaData


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("scstreeannotator", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: scstreeannotator test.trees out.txt");
        progressStream.println("  Example: scstreeannotator -burnin 10 -heights mean test.trees out.txt");
        progressStream.println("  Example: scstreeannotator -burnin 20 -target map.tree test.trees out.txt");
        progressStream.println();
    } // printUsage


    //***********************************************
    //*                 Main method                 *
    //***********************************************

    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        String targetTreeFileName = null;
        String inputFileName = null;
        String outputFileName = null;
        String simpleOutputFileName = null;

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

            // TODO: set up welcome message
            final String versionString = version.getVersionString();
            String nameString = "ScsTreeAnnotator, based on TreeAnnotator " + versionString;
            String aboutString = "<html><center><p>" + versionString + ", " + version.getDateString() + "</p>" +
                    "<p>by<br>" +
                    "Andrew Rambaut and Alexei J. Drummond</p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
                    "</center></html>";

            new ConsoleApplication(nameString, aboutString, icon, true);
            Log.info = System.out;
            Log.err = System.err;

            // The ConsoleApplication will have overridden System.out so set progressStream
            // to capture the output to the window:
            progressStream = System.out;

            printTitle();

            ScsTreeAnnotatorDialog dialog = new ScsTreeAnnotatorDialog(false, new JFrame(), null);

            if (!dialog.showDialog("ScsTreeAnnotator " + versionString)) {
                return;
            }

            int burninPercentage = dialog.getBurninPercentage();
            if (burninPercentage < 0) {
                Log.warning.println("burnin percentage is " + burninPercentage + " but should be non-negative. Setting it to zero");
                burninPercentage = 0;
            }
            if (burninPercentage >= 100) {
                Log.err.println("burnin percentage is " + burninPercentage + " but should be less than 100.");
                return;
            }
            double posteriorLimit = dialog.getPosteriorLimit();
            double hpd2D = 0.80;
            Target targetOption = dialog.getTargetOption();
            HeightsSummary heightsOption = dialog.getHeightsOption();

            targetTreeFileName = dialog.getTargetFileName();
            if (targetOption == Target.USER_TARGET_TREE && targetTreeFileName == null) {
                Log.err.println("No target file specified");
                return;
            }

            inputFileName = dialog.getInputFileName();
            if (inputFileName == null) {
                Log.err.println("No input file specified");
                return;
            }

            outputFileName = dialog.getOutputFileName();
            if (outputFileName == null) {
                Log.err.println("No output file specified");
                return;
            }

            boolean lowMem = dialog.useLowMem();
            saveSimpleTree = dialog.saveSimpleTree();
            hasTrunk = dialog.hasTrunk();

            if (saveSimpleTree) {
                String delimiter = String.valueOf(File.separatorChar);
                String[] tmp1 = outputFileName.split(delimiter);
                String[] tmp2 = tmp1[tmp1.length - 1].split("\\.");
                tmp1[tmp1.length - 1] = tmp2[0] + "_simple." + tmp2[1];
                simpleOutputFileName = String.join(delimiter, tmp1);
            }

            try {
                new ScsTreeAnnotator(burninPercentage,
                        lowMem,
                        heightsOption,
                        posteriorLimit,
                        hpd2D,
                        targetOption,
                        targetTreeFileName,
                        inputFileName,
                        outputFileName,
                        simpleOutputFileName);

            } catch (Exception ex) {
                Log.err.println("Exception: " + ex.getMessage());
            }

            progressStream.println("Finished - Quit program to exit.");
            while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("prefix", "output_file_prefix", "specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL"),

                        //new Arguments.StringOption("target", new String[] { "maxclade", "maxtree" }, false, "an option of 'maxclade' or 'maxtree'"),
                        new Arguments.StringOption("heights", new String[]{"keep", "median", "mean", "ca"}, false,
                                "an option of 'keep' (default), 'median', 'mean' or 'ca' -> OPTIONAL"),
                        new Arguments.IntegerOption("burnin", 0, 99, "the percentage of states to be considered as 'burn-in' -> MANDATORY, the same as -b\""),
                        // allow -b as burnin option, just like other apps
                        new Arguments.IntegerOption("b", 0, 99, "the percentage of states to be considered as 'burn-in' -> MANDATORY, the same as -burnin"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annotated -> OPTIONAL"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated -> OPTIONAL"),
                        new Arguments.Option("help", "option to print this message -> OPTIONAL"),
                        new Arguments.Option("forceDiscrete", "forces integer traits to be treated as discrete traits -> OPTIONAL"),
                        new Arguments.Option("lowMem", "use less memory, which is a bit slower -> OPTIONAL"),
                        new Arguments.RealOption("hpd2D", "the HPD interval to be used for the bivariate traits -> OPTIONAL"),
                        new Arguments.Option("nohpd2D", "suppress calculation of HPD intervals for the bivariate traits -> OPTIONAL"),
                        new Arguments.Option("noSA", "interpret the tree set as begin from a not being from a sampled ancestor analysis, even if there are zero branch lengths in the tree set -> OPTIONAL"),

                        // some added arguments
                        new Arguments.Option("simpleTree", "simple output tree only containing tree heights (the output file will be labeled with 'simple'), default disabled -> OPTIONAL"),
                        new Arguments.Option("noTrunk", "trees processed do not have a trunk connecting the root of samples and the root of entire tree, default disabled -> OPTIONAL")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("prefix")) {
            System.setProperty("tree.annotation.file.prefix", arguments.getStringOption("prefix").trim());
        }

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

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        boolean lowMem = false;
        if (arguments.hasOption("lowMem")) {
            lowMem = true;
        }

        HeightsSummary heights = HeightsSummary.CA_HEIGHTS;
        if (arguments.hasOption("heights")) {
            String value = arguments.getStringOption("heights");
            if (value.equalsIgnoreCase("mean")) {
                heights = HeightsSummary.MEAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("median")) {
                heights = HeightsSummary.MEDIAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("ca")) {
                heights = HeightsSummary.CA_HEIGHTS;
                Log.info.println("Please cite: Heled and Bouckaert: Looking for trees in the forest:\n" +
                        "summary tree from posterior samples. BMC Evolutionary Biology 2013 13:221.");
            }
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        } else if (arguments.hasOption("b")) {
            burnin = arguments.getIntegerOption("b");
        }
        if (burnin >= 100) {
            Log.err.println("burnin percentage is " + burnin + " but should be less than 100.");
            System.exit(1);
        }

        double posteriorLimit = 0.0;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }

        double hpd2D = 0.80;
        if (arguments.hasOption("hpd2D")) {
            hpd2D = arguments.getRealOption("hpd2D");
            if (hpd2D <= 0 || hpd2D >= 1) {
                Log.err.println("hpd2D is a fraction and should be in between 0.0 and 1.0.");
                System.exit(1);
            }
            processBivariateAttributes = true;
        }

        Target target = Target.MAX_CLADE_CREDIBILITY;
        if (arguments.hasOption("target")) {
            target = Target.USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
        }

        if (arguments.hasOption("simpleTree")) {
            saveSimpleTree = true;
        }

        if (arguments.hasOption("noTrunk")) {
            Log.warning.println("are you sure your sampled trees do not contain an extra trunk?");
            hasTrunk = false;
        }

        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 2:
                outputFileName = args2[1];
                // fall to
            case 1:
                inputFileName = args2[0];
                break;
            default: {
                Log.err.println("Unknown option: " + args2[2]);
                Log.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        if (saveSimpleTree) {
            assert outputFileName != null;
            simpleOutputFileName = FileNameProcessor.getRootPath(outputFileName) +
                    FileNameProcessor.getBaseName(outputFileName) +
                    "_simple" +
                    FileNameProcessor.getSuffix(outputFileName);
        }

        try {
            new ScsTreeAnnotator(
                    burnin,
                    lowMem,
                    heights,
                    posteriorLimit,
                    hpd2D,
                    target,
                    targetTreeFileName,
                    inputFileName,
                    outputFileName,
                    simpleOutputFileName
            );
        } catch (IOException e) {
            throw e;
        } catch (Exception e) {
            e.printStackTrace();
        }
    } // main

}
