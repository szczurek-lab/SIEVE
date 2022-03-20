package beast.app.variantcallingmodeladaptor;

import beast.app.tools.LogCombiner;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.app.variantcaller.AllelicSeqLogProcessor;
import beast.app.variantcaller.EstimatesTypeCollection;
import beast.app.variantcallingmodeladaptor.seqcovmodelmodifier.EstimatedSharedAllelicSeqCovModelModifier;
import beast.app.variantcallingmodeladaptor.seqcovmodelmodifier.ExploredSharedAllelicSeqCovModelModifier;
import beast.app.variantcallingmodeladaptor.seqcovmodelmodifier.PrivateAllelicSeqCovModelModifier;
import beast.app.variantcallingmodeladaptor.seqcovmodelmodifier.SeqCovModelModifier;
import beast.core.Description;
import beast.core.util.Log;
import beast.util.FileNameProcessor;
import jam.console.ConsoleApplication;
import org.jetbrains.annotations.NotNull;
import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.swing.*;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * The original model configuration document should be modified in the following aspects:
 * <p>
 * 1. Removal:
 * 1.1 Anything under 'run' tag including child nodes and attributes except:
 * 1.1.1 'state' tag
 * 1.1.2 'distribution' tag
 * 1.1.3 'id' attribute
 * 1.1.4 'spec' attribute
 * <p>
 * 2. Modification:
 * 2.1 Replace the value of 'tree' tag under 'state' tag with the path of the input best tree.
 * 2.2 Add a new attribute 'prioritizeMedianRate' to 'tree' tag if it is of 'ScsTree' type
 * 2.3 Replace the value of other 'parameter' tags under 'state' tag with the estimates of MCMC samples.
 * 2.4 'seqCovModel' tag:
 * 2.4.1 Set or add 'inVariantCallingMode' with value 'true'.
 * 2.4.2 'ReducedSeqCovModel' type: Nothing else needs to be done.
 * 2.4.3 'ExploredSharedAllelicSeqCovModel' type:
 * 2.4.3.1 Add an attribute 'allelicSeqCovArray' with the value of 'allelicSeqCov' tag.
 * 2.4.3.2 Add an attribute 'allelicSeqCovRawVarArray' with the value of 'allelicSeqCovRawVar' tag.
 * 2.4.3.3 Remove 'allelicSeqCov' and 'allelicSeqCovRawVar' tags along with their priors, if applicable.
 * 2.4.4 'PrivateAllelicSeqCovModel' type:
 * 2.4.4.1 Remove 'allelicSeqCov' and 'allelicSeqCovRawVar' attributes.
 * 2.4.4.2 Add 'allelicSeqCovArray' and 'allelicSeqCovRawVarArray' attributes with values collected from the estimates
 * (comma-separated within a matrix and semicolon-separated between matrices).
 * 2.5 Edit 'run' tag:
 * 2.5.1 Reset 'spec' attribute to 'beast.core.VariantCall';
 * 2.5.2 Reset 'id' attribute to 'VariantCall';
 * 2.6 Add a new child named 'variantsInfo' to 'treeLikelihood' tag.
 * 2.6.1 if 'treeLikelihood' is of type 'ThreadedScsThreeLikelihood', then 'variantsInfo' should be of type
 * 'beast.evolution.variantsinfo.ThreadedVariantsInfoVCF'.
 * 2.6.2 if 'treeLikelihood' is of type 'ScsThreeLikelihood', then 'variantsInfo' should be of type
 * 'beast.evolution.variantsinfo.VariantsInfoVCF'.
 * 2.7 Add attributes 'useOnlyBranchLength' and 'meanRate' to 'treeLikelihood' tag.
 * 2.8 Add a new child node named 'logger' to 'run' tag, and take 'variantsInfo' tag as input.
 * 2.9 Make sure the 'data' tag has no bias correction used, i.e., having value of 'ascertained' attribute to be 'none'.
 * <p>
 * 3. Key:
 * 3.1 Get values of 'lociInfo' tag in order to match the estimates of allelic sequencing information.
 */
@Description("Adapt the model configuration to variant calling")
public class VariantCallingModelAdaptor {

    final static private Pattern ID_REF_PATTERN = Pattern.compile("@(.+)");


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private DocumentBuilderFactory dbFac;
    private Document doc;

    private int numOfAdaptedConfigFiles = 0;
    private List<String> adaptedConfigFileNames = null;
    private List<String> vcFileNames = null;
    private List<EstimatesTypeCollection> adaptedConfigEstimateTypes = null;

    public static PrintStream progressStream = Log.err;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    /**
     * Called from GUI or CMD independently.
     *
     * @param configFileName      configuration file name
     * @param inputTreeFileName   input best tree file name
     * @param estimatesFileName   cached estimates file name
     * @param outputFileName      output file name (can be null)
     * @param saveDetails         whether to save ternary files or not
     * @param cellThreshold       number of mutated cells as a threshold to filter a invariant site
     * @param useMeanRate         whether to use mean estimates of rate or median
     * @param useOnlyBranchLength whether to only use branch length from the input tree or to use both branch length and rate
     */
    public VariantCallingModelAdaptor(
            final String configFileName,
            final String inputTreeFileName,
            final String estimatesFileName,
            final String outputFileName,
            final boolean saveDetails,
            final int cellThreshold,
            final boolean useMeanRate,
            final boolean useOnlyBranchLength
    ) {
        Map<String, Map<String, Double>> mcmcSamplesEstimates = new HashMap<>();
        AllelicSeqLogProcessor allelicSeqLog = null;
        EstimatesTypeCollection estimatesTypeCollection = new EstimatesTypeCollection();

        // Get cached estimates
        try {
            allelicSeqLog = new AllelicSeqLogProcessor(
                    estimatesFileName,
                    System.err,
                    mcmcSamplesEstimates,
                    estimatesTypeCollection
            );
        } catch (Exception e) {
            e.printStackTrace();
            progressStream.println("Error parsing cached estimates file: " + e.getMessage());
            return;
        }

        // List the output file
        listOutputFileNames(configFileName, outputFileName, 1, estimatesTypeCollection);

        // Load configuration file
        this.dbFac = DocumentBuilderFactory.newInstance();
        this.dbFac.setIgnoringElementContentWhitespace(true);
        try {
            this.doc = this.dbFac.newDocumentBuilder().parse(new File(configFileName));
        } catch (SAXException | IOException | ParserConfigurationException e) {
            e.printStackTrace();
            progressStream.println("Error parsing configuration file: " + e.getMessage());
            return;
        }
        this.doc.getDocumentElement().normalize();

        // Remove some tags
        removeTags(
                "run",
                Stream.of("state", "distribution").collect(Collectors.toList()),
                Stream.of("id", "spec", "name").collect(Collectors.toList())
        );

        // Get 'lociInfo' tag
        String[] lociInfo = getTagValue("lociInfo").trim().split(";");

        // Modify tags
        try {
            modifyTags(
                    inputTreeFileName,
                    lociInfo,
                    mcmcSamplesEstimates,
                    allelicSeqLog.getAllelicSeqEstimates(),
                    allelicSeqLog.getSitesMap(),
                    saveDetails,
                    cellThreshold,
                    useMeanRate,
                    useOnlyBranchLength
            );
        } catch (TransformerException e) {
            e.printStackTrace();
            progressStream.println("Something is wrong when modifying configuration file.");
        }
    }

    /**
     * Called by other classes as a component.
     *
     * @param configFileName          configuration file name
     * @param inputTreeFileName       input best tree file name
     * @param outputFileName          output file name (can be null)
     * @param estimatesTypeCollection median, mean, mode, or all?
     * @param mcmcSamplesEstimates    estimates of model parameters from mcmc samples
     * @param allelicSeqInfo     estimates of allelic sequencing information
     * @param saveDetails             whether to save to ternary or not
     * @param cellThreshold           number of mutated cells as a threshold to filter a invariant site
     * @param useMeanRate         whether to use mean estimates of rate or median
     * @param useOnlyBranchLength whether to only use branch length from the input tree or to use both branch length and rate
     */
    public VariantCallingModelAdaptor(
            final String configFileName,
            final String inputTreeFileName,
            final String outputFileName,
            final EstimatesTypeCollection estimatesTypeCollection,
            final Map<String, Map<String, Double>> mcmcSamplesEstimates,
            final AllelicSeqLogProcessor allelicSeqInfo,
            final boolean saveDetails,
            final int cellThreshold,
            final boolean useMeanRate,
            final boolean useOnlyBranchLength
    ) {
        // List the output file
        listOutputFileNames(configFileName, outputFileName, 0, estimatesTypeCollection);

        // Load configuration file
        this.dbFac = DocumentBuilderFactory.newInstance();
        this.dbFac.setIgnoringElementContentWhitespace(true);
        try {
            this.doc = this.dbFac.newDocumentBuilder().parse(new File(configFileName));
        } catch (SAXException | IOException | ParserConfigurationException e) {
            e.printStackTrace();
            progressStream.println("Error parsing configuration file: " + e.getMessage());
            return;
        }
        this.doc.getDocumentElement().normalize();

        // Remove some tags
        removeTags(
                "run",
                Stream.of("state", "distribution").collect(Collectors.toList()),
                Stream.of("id", "spec", "name").collect(Collectors.toList())
        );

        // Get 'lociInfo' tag
        String[] lociInfo = getTagValue("lociInfo").trim().split(";");

        // Modify tags
        try {
            modifyTags(
                    inputTreeFileName,
                    lociInfo,
                    mcmcSamplesEstimates,
                    allelicSeqInfo == null ? null : allelicSeqInfo.getAllelicSeqEstimates(),
                    allelicSeqInfo == null ? null : allelicSeqInfo.getSitesMap(),
                    saveDetails,
                    cellThreshold,
                    useMeanRate,
                    useOnlyBranchLength
            );
        } catch (Exception e) {
            e.printStackTrace();
            progressStream.println("Something is wrong when modifying configuration file.");
        }
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * Prepare output file names.
     *
     * @param configFileName          template configuration file name
     * @param outputFileName          output configuration file name, which defines the root path of output file names if not null
     * @param baseNameSource          where does base name come from? 0 from {configFileName}; 1 from {outputFileName}
     * @param estimatesTypeCollection apparently
     */
    public void listOutputFileNames(
            final String configFileName,
            final String outputFileName,
            final int baseNameSource,
            final EstimatesTypeCollection estimatesTypeCollection
    ) {
        // Sanity check
        if (baseNameSource != 0 && baseNameSource != 1)
            throw new IllegalArgumentException("'baseNameSource' should either be 0 or be 1.");

        if (baseNameSource == 1 && outputFileName == null)
            throw new IllegalArgumentException("'baseNameSource' cannot come from an undefined 'outputFileName'.");

        if (this.adaptedConfigFileNames == null)
            this.adaptedConfigFileNames = new ArrayList<>();

        if (this.adaptedConfigEstimateTypes == null)
            this.adaptedConfigEstimateTypes = new ArrayList<>();

        if (this.vcFileNames == null)
            this.vcFileNames = new ArrayList<>();

        final String prefix = (baseNameSource == 0 ?
                FileNameProcessor.getBaseName(configFileName) : FileNameProcessor.getBaseName(outputFileName)) +
                "_";
        final String suffix = FileNameProcessor.getSuffix(configFileName);
        final String rootPath = FileNameProcessor.getRootPath(outputFileName);
        String tmp;

        if (estimatesTypeCollection.getEstimatesType() == EstimatesTypeCollection.EstimatesType.ALL) {

            this.adaptedConfigEstimateTypes.add(new EstimatesTypeCollection(EstimatesTypeCollection.EstimatesType.MEDIAN, null));
            tmp = prefix + EstimatesTypeCollection.EstimatesType.MEDIAN.toString().toLowerCase();
            this.adaptedConfigFileNames.add(rootPath + tmp + suffix);
            this.vcFileNames.add(tmp + ".vcf");
            this.numOfAdaptedConfigFiles++;

            this.adaptedConfigEstimateTypes.add(new EstimatesTypeCollection(EstimatesTypeCollection.EstimatesType.MEAN, null));
            tmp = prefix + EstimatesTypeCollection.EstimatesType.MEAN.toString().toLowerCase();
            this.adaptedConfigFileNames.add(rootPath + tmp + suffix);
            this.vcFileNames.add(tmp + ".vcf");
            this.numOfAdaptedConfigFiles++;

            for (EstimatesTypeCollection.ModeKDEType i : EstimatesTypeCollection.ModeKDEType.values()) {
                this.adaptedConfigEstimateTypes.add(new EstimatesTypeCollection(EstimatesTypeCollection.EstimatesType.MODE, i));
                tmp = prefix + "mode_" + i.toString().toLowerCase();
                this.adaptedConfigFileNames.add(rootPath + tmp + suffix);
                this.vcFileNames.add(tmp + ".vcf");
                this.numOfAdaptedConfigFiles++;
            }

        } else if (estimatesTypeCollection.getEstimatesType() == EstimatesTypeCollection.EstimatesType.MODE) {

            this.adaptedConfigEstimateTypes.add(new EstimatesTypeCollection(estimatesTypeCollection.getEstimatesType(), estimatesTypeCollection.getModeKDEType()));
            tmp = prefix + "_mode_" + estimatesTypeCollection.getModeKDEType().toString().toLowerCase();
            this.adaptedConfigFileNames.add(rootPath + tmp + suffix);
            this.vcFileNames.add(tmp + ".vcf");
            this.numOfAdaptedConfigFiles = 1;

        } else {

            this.adaptedConfigEstimateTypes.add(new EstimatesTypeCollection(estimatesTypeCollection.getEstimatesType(), null));
            tmp = prefix + estimatesTypeCollection.getEstimatesType().toString().toLowerCase();
            this.adaptedConfigFileNames.add(rootPath + tmp + suffix);
            this.vcFileNames.add(tmp + ".vcf");
            this.numOfAdaptedConfigFiles = 1;

        }

    } // listOutputFileNames

    /**
     * Remove all children under 'name' tag except for those in exclusion.
     *
     * @param name               parent tag
     * @param excludedTags       exceptions of child nodes
     * @param excludedAttributes exceptions of attributes
     */
    public void removeTags(String name, List<String> excludedTags, List<String> excludedAttributes) {
        Node target;
        NodeList children;

        // Get the 'name' tag and its children
        NodeList candidates = this.doc.getElementsByTagName(name);
        if (candidates == null || candidates.getLength() == 0) {
            progressStream.println("No '" + name + "' tags found.");
            return;
        } else if (candidates.getLength() > 1) {
            progressStream.println("More than one '" + name + "' tags found, while only one expected.");
            return;
        } else {
            target = candidates.item(0);
            children = target.getChildNodes();
        }

        // Process children tags one by one
        Set<String> removalNodes = new HashSet<>();
        for (int i = 0; i < children.getLength(); i++) {
            Node child = children.item(i);
            if ((excludedTags == null || excludedTags.size() == 0 || !excludedTags.contains(child.getNodeName())) &&
                    child.getNodeType() == Node.ELEMENT_NODE)
                removalNodes.add(child.getNodeName());
        }
        for (String i : removalNodes)
            removeTags(i, true);

        // Process attributes of 'name' tag one by one
        Set<String> removalAttrs = new HashSet<>();
        NamedNodeMap attributes = target.getAttributes();
        for (int i = 0; i < attributes.getLength(); i++) {
            Node attr = attributes.item(i);
            if (excludedAttributes == null || excludedAttributes.size() == 0 ||
                    !excludedAttributes.contains(attr.getNodeName()))
                removalAttrs.add(attr.getNodeName());
        }
        for (String i : removalAttrs)
            ((Element) target).removeAttribute(i);
    } // removeTags

    public void removeTags(final String name, boolean includeItself) {
        Set<Node> removal = new HashSet<>();

        // Search for 'name' tags
        NodeList nodes = this.doc.getElementsByTagName(name);
        if (nodes == null || nodes.getLength() == 0) {
            progressStream.println("No '" + name + "' tags found.");
        } else {

            // Collect 'name' tags and all their children
            for (int i = 0; i < nodes.getLength(); i++) {
                removal.addAll(collectNodesUntilPreviousElement(nodes.item(i)));
                removal.addAll(collectChildren(nodes.item(i)));

                if (!includeItself) {
                    removal.remove(nodes.item(i));
                }
            }

            // Remove 'name' tags and all their children
            if (removal.size() > 0) {
                for (Node i : removal) {
                    i.getParentNode().removeChild(i);
                }
            }

        }
    } // removeTags

    public String getTagValue(String name) {
        NodeList nodes = this.doc.getElementsByTagName(name);
        if (nodes == null || nodes.getLength() == 0) {
            progressStream.println("No '" + name + "' found.");
            return null;
        } else {
            if (nodes.getLength() == 1) {

                // Search for 'value' attribute
                NamedNodeMap attributes = nodes.item(0).getAttributes();
                for (int i = 0; i < attributes.getLength(); i++) {
                    if (attributes.item(i).getNodeName().equals("value")) {
                        return attributes.item(i).getNodeValue();
                    }
                }

                // Use child node as value
                NodeList children = nodes.item(0).getChildNodes();
                if (children.getLength() == 0) {
                    progressStream.println("No value for '" + name + "' found.");
                    return null;
                } else if (children.getLength() > 1) {
                    progressStream.println("No value for '" + name + "' found.");
                    return null;
                } else {
                    return children.item(0).getNodeValue();
                }
            } else {
                progressStream.println("More than one '" + name + "' found.");
                return null;
            }
        }
    } // getTagValue

    /**
     * Include the argument itself.
     *
     * @param node target node
     * @return a list of children nodes including the argument itself
     */
    public Set<Node> collectChildren(final Node node) {
        Set<Node> results = new HashSet<>();

        if (node != null) {
            results.add(node);

            NodeList children = node.getChildNodes();
            if (children != null && children.getLength() > 0) {
                for (int i = 0; i < children.getLength(); i++) {
                    results.addAll(collectChildren(children.item(i)));
                }
            }
        }

        return results;
    } // collectChildren

    /**
     * Collect comments and empty tags between @node and its most recent element sibling.
     *
     * @param node apparently
     * @return a list of nodes
     */
    public Set<Node> collectNodesUntilPreviousElement(Node node) {
        Set<Node> results = new HashSet<>();

        while (node.getPreviousSibling() != null && node.getPreviousSibling().getNodeType() != Node.ELEMENT_NODE) {
            node = node.getPreviousSibling();
            results.add(node);
        }

        return results;
    } // collectPreviousComments

    private void modifyTags(
            final String inputTreeFileName,
            String[] lociInfo,
            final Map<String, Map<String, Double>> mcmcSamplesEstimates,
            final Map<String, Map<String, double[][]>> allelicSeqEstimates,
            final Map<String, String> sitesMap,
            final boolean saveDetails,
            final int cellThreshold,
            final boolean useMeanRate,
            final boolean useOnlyBranchLength
    ) throws TransformerException {
        SeqCovModelModifier seqCovModelModifier = null;
        List<Node> paraNodes = new ArrayList<>();

        // Modify 'run' tag by redefining 'id' and 'spec' attributes.
        NodeList runTags = this.doc.getElementsByTagName("run");
        Node runTag;
        if (runTags == null || runTags.getLength() == 0) {
            progressStream.println("No 'run' tags found.");
            return;
        } else if (runTags.getLength() > 1) {
            progressStream.println("More than one 'run' tags found, while only one expected.");
            return;
        } else {
            runTag = runTags.item(0);
            NamedNodeMap runAttributes = runTag.getAttributes();
            for (int i = 0; i < runAttributes.getLength(); i++) {
                Node runAttr = runAttributes.item(i);
                if (runAttr.getNodeName().equals("id"))
                    runAttr.setNodeValue("VariantCall");
                else if (runAttr.getNodeName().equals("spec"))
                    runAttr.setNodeValue("beast.core.VariantCall");
            }
        }

        // Modifying 'state' tag under 'run' tag by only keeping 'id', 'spec' and/or 'name' attributes.
        NodeList stateTags = this.doc.getElementsByTagName("state");
        if (stateTags == null || stateTags.getLength() == 0) {
            progressStream.println("No 'state' tags found.");
            return;
        } else if (runTags.getLength() > 1) {
            progressStream.println("More than one 'state' tags found, while only one expected.");
            return;
        } else {
            Node stateTag = stateTags.item(0);

            if (!stateTag.getParentNode().getNodeName().equals("run")) {
                progressStream.println("'state' tag should be a child of 'run' tag.");
                return;
            } else {
                NamedNodeMap stateAttributes = stateTag.getAttributes();
                for (int i = 0; i < stateAttributes.getLength(); i++) {
                    Node stateAttr = stateAttributes.item(i);
                    if (!stateAttr.getNodeName().equals("id") &&
                            !stateAttr.getNodeName().equals("name") &&
                            !stateAttr.getNodeName().equals("spec"))
                        ((Element) stateTag).removeAttribute(stateAttr.getNodeName());
                }
            }
        }

        // Modify 'tree' tag.
        // Only keep attributes of 'id', 'name', 'nodeType' and 'spec', if any exist.
        // Add a new attribute 'prioritizeMedianRate'.
        NodeList treeNodes = this.doc.getElementsByTagName("tree");
        if (treeNodes == null || treeNodes.getLength() == 0) {
            progressStream.println("No 'tree' tag found.");
            return;
        } else if (treeNodes.getLength() > 1) {
            progressStream.println("Error! More than one 'tree' tags found. Only one allowed.");
            return;
        } else {
            Element treeNode = (Element) treeNodes.item(0);

            // Remove useless attributes
            NamedNodeMap treeAttributes = treeNode.getAttributes();
            for (int i = 0; i < treeAttributes.getLength(); i++) {
                Node attr = treeAttributes.item(i);
                if (!attr.getNodeName().equals("id") &&
                        !attr.getNodeName().equals("name") &&
                        !attr.getNodeName().equals("nodetype") &&
                        !attr.getNodeName().equals("spec"))
                    treeNode.removeAttribute(attr.getNodeName());
            }

            // Remove all children
            removeTags("tree", false);

            // Add a new attribute 'treeFileName'
            treeNode.setAttribute("treeFileName", inputTreeFileName);

            // Add a new attribute 'prioritizeMedianRate'
            treeNode.setAttribute(
                    "prioritizeMedianRate",
                    useMeanRate ? "false" : "true"
            );
        }

        // Collect 'parameter' tags under 'state' tag
        NodeList stateNodes = this.doc.getElementsByTagName("state");
        if (stateNodes == null || stateNodes.getLength() == 0) {
            progressStream.println("No 'state' tag found.");
            return;
        } else if (stateNodes.getLength() > 1) {
            progressStream.println("Error! More than one 'state' tags found. Only one allowed.");
            return;
        } else {
            NodeList stateChildren = stateNodes.item(0).getChildNodes();
            for (int i = 0; i < stateChildren.getLength(); i++) {
                Node child = stateChildren.item(i);
                if (child.getNodeType() == Node.ELEMENT_NODE &&
                        child.getNodeName().equals("parameter"))
                    paraNodes.add(child);
            }
        }

        // Pre-modify 'seqCovModel' tag
        NodeList seqCovModelNodes = this.doc.getElementsByTagName("seqCovModel");
        Element seqCovModelNode;
        if (seqCovModelNodes == null || seqCovModelNodes.getLength() == 0) {
            progressStream.println("No 'seqCovModel' tag found.");
            return;
        } else if (seqCovModelNodes.getLength() > 1) {
            progressStream.println("Error! More than one 'seqCovModel' tags found. Only one allowed.");
            return;
        } else {
            seqCovModelNode = (Element) seqCovModelNodes.item(0);
            final String seqCovModelType = seqCovModelNode.getAttribute("spec");

            if (seqCovModelType.equals("PrivateAllelicSeqCovModel"))
                seqCovModelModifier = new PrivateAllelicSeqCovModelModifier(
                        this.doc,
                        seqCovModelNode
                );

            if (seqCovModelType.equals("EstimatedSharedAllelicSeqCovModel"))
                seqCovModelModifier = new EstimatedSharedAllelicSeqCovModelModifier();

            if (seqCovModelType.equals("ExploredSharedAllelicSeqCovModel"))
                seqCovModelModifier = new ExploredSharedAllelicSeqCovModelModifier(
                        this.doc,
                        seqCovModelNode,
                        paraNodes
                );

            if (seqCovModelModifier != null)
                seqCovModelModifier.preModify();
        }

        // Get the 'substModel' tag.
//        NodeList substModelNodes = this.doc.getElementsByTagName("substModel");
//        Element substModelNode;
//        if (substModelNodes == null || substModelNodes.getLength() == 0) {
//            progressStream.println("No 'substModel' tag found.");
//            return;
//        } else if (substModelNodes.getLength() > 1) {
//            progressStream.println("Error! More than one 'substModel' tags found. Only one allowed.");
//            return;
//        } else {
//            substModelNode = (Element) substModelNodes.item(0);
//
//            if (!substModelNode.getParentNode().getNodeName().equals("siteModel")) {
//                progressStream.println("Error! 'substModel' tag should be a child of 'siteModel' tag.");
//                return;
//            }
//        }

        // Add a sibling to 'rawReadCountsModel' tag (or a child to its parent tag).
        Node treeLikelihoodNode = seqCovModelNode.getParentNode().getParentNode();
        String treeLikelihoodClassName = null;
        NamedNodeMap treeLikelihoodAttributes = treeLikelihoodNode.getAttributes();
        for (int i = 0; i < treeLikelihoodAttributes.getLength(); i++) {
            Node tmp = treeLikelihoodAttributes.item(i);
            if (tmp.getNodeName().trim().equals("spec")) {
                treeLikelihoodClassName = tmp.getNodeValue().trim();
                break;
            }
        }
        assert treeLikelihoodClassName != null;
        Element variantInfoNode = this.doc.createElement("variantsInfo");
        variantInfoNode.setAttribute("id", "variantsInfo");
        variantInfoNode.setAttribute("seqCovModel", "@" + seqCovModelNode.getAttribute("id"));
        if (((Element) treeLikelihoodNode).hasAttribute("scsData"))
            variantInfoNode.setAttribute("scsData", ((Element) treeLikelihoodNode).getAttribute("scsData"));
        else
            throw new IllegalArgumentException("Error: 'treeLikelihood' tag should contain 'scsData' attribute.");
        if (((Element) treeLikelihoodNode).hasAttribute("tree"))
            variantInfoNode.setAttribute("tree", ((Element) treeLikelihoodNode).getAttribute("tree"));
        else
            throw new IllegalArgumentException("Error: 'treeLikelihood' tag should contain 'tree' attribute.");
        if (treeLikelihoodClassName.contains("ThreadedScsTreeLikelihood"))
            variantInfoNode.setAttribute("spec", "beast.evolution.variantsinfo.ThreadedVariantsInfoVCF");
        else if (treeLikelihoodClassName.contains("ScsTreeLikelihood"))
            variantInfoNode.setAttribute("spec", "beast.evolution.variantsinfo.VariantsInfoVCF");
        else
            throw new IllegalArgumentException("Error: 'treeLikelihood' in the configuration file is not of type " +
                    "'ThreadedScsTreeLikelihood' or 'ScsTreeLikelihood'.");
//        if (substModelNode.hasAttribute("id"))
//            variantInfoNode.setAttribute("substModel", "@" + substModelNode.getAttribute("id"));
//        else
//            throw new IllegalArgumentException("Error: 'substModel' in the configuration file should have an " +
//                    "attribute named 'id'.");
        variantInfoNode.setAttribute("cellThreshold", String.valueOf(cellThreshold));
        treeLikelihoodNode.appendChild(variantInfoNode);

        // Add attribute 'useOnlyBranchLength' to 'treeLikelihood' tag.
        ((Element) treeLikelihoodNode).setAttribute(
                "useOnlyBranchLength",
                useOnlyBranchLength ? "true" : "false"
        );

        // Add attribute 'meanRate' to 'treeLikelihood' tag.
        NodeList branchRateModelNodes = this.doc.getElementsByTagName("branchRateModel");
        Node branchRateModelNode;
        if (branchRateModelNodes == null || branchRateModelNodes.getLength() == 0) {
            progressStream.println("No 'branchRateModel' tag found.");
            return;
        } else if (branchRateModelNodes.getLength() > 1) {
            progressStream.println("Error! More than one 'branchRateModel' tags found. Only one allowed.");
            return;
        } else {
            branchRateModelNode = branchRateModelNodes.item(0);
            final String meanRateID = getAttribute(branchRateModelNode, "clock.rate");
            if (meanRateID != null)
                ((Element) treeLikelihoodNode).setAttribute("meanRate", "@" + meanRateID);
        }

        // Add a new 'logger' child tag to 'run' tag
        Element variantsLoggerNode = this.doc.createElement("variantsLogger");
        variantsLoggerNode.setAttribute("id", "variantsLogger");
        variantsLoggerNode.setAttribute("spec", "beast.core.variantslogger.VariantsLogger");
        variantsLoggerNode.setAttribute("log", "@variantsInfo");
        variantsLoggerNode.setAttribute("saveDetails", saveDetails ? "true" : "false");
        runTag.appendChild(variantsLoggerNode);

        // Modify tags
        for (int i = 0; i < this.numOfAdaptedConfigFiles; i++) {

            // Modify or add 'fileName' attribute of or to 'logger' tag
            variantsLoggerNode.setAttribute("fileName", this.vcFileNames.get(i));

            // Modify 'parameter' tags
            for (Node para : paraNodes) {
                String id = ((Element) para).getAttribute("id");

                if (!mcmcSamplesEstimates.containsKey(id)) continue;

                Map<String, Double> samples = mcmcSamplesEstimates.get(id);
                String value;
                if (adaptedConfigEstimateTypes.get(i).getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEDIAN ||
                        adaptedConfigEstimateTypes.get(i).getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEAN)
                    value = String.valueOf(samples.get(adaptedConfigEstimateTypes.get(i).getEstimatesType().toString().toLowerCase()));
                else
                    value = String.valueOf(samples.get(adaptedConfigEstimateTypes.get(i).getModeKDEType().toString().toLowerCase()));

                // Set value through 'value' attribute or child node
                if (((Element) para).hasAttribute("value"))
                    ((Element) para).setAttribute("value", value);
                else {
                    NodeList children = para.getChildNodes();
                    if (children.getLength() == 1 && children.item(0).getNodeType() == Node.TEXT_NODE)
                        children.item(0).setNodeValue(value);
                    else {
                        // Remove all children
                        for (int j = 0; j < children.getLength(); j++)
                            para.removeChild(children.item(j));

                        // Add new text node
                        Text text = this.doc.createTextNode(value);
                        para.appendChild(text);
                    }
                }
            }

            // Modify 'seqCovModel' tag
            if (seqCovModelModifier != null)
                seqCovModelModifier.postModify(
                        lociInfo,
                        allelicSeqEstimates,
                        sitesMap,
                        adaptedConfigEstimateTypes.get(i)
                );

            // Save modified configuration file to disk
            Transformer tf = TransformerFactory.newInstance().newTransformer();
            tf.setOutputProperty(OutputKeys.INDENT, "yes");
            DOMSource src = new DOMSource(this.doc);
            if (System.getProperty("variant.calling.file.prefix") != null) {
                adaptedConfigFileNames.set(i, System.getProperty("variant.calling.file.prefix") + adaptedConfigFileNames.get(i));
            }
            File outputConfigFile = new File(adaptedConfigFileNames.get(i));
            if (adaptedConfigFileNames.get(i).contains("/")) {
                outputConfigFile.getParentFile().mkdirs();
            }
            StreamResult sr = new StreamResult(outputConfigFile);
            tf.transform(src, sr);

            // Print some information
            progressStream.print(">>> ");
            if (adaptedConfigEstimateTypes.get(i).getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEDIAN ||
                    adaptedConfigEstimateTypes.get(i).getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEAN) {
                progressStream.print(adaptedConfigEstimateTypes.get(i).getEstimatesType().toString() + " estimates ");
            } else {
                progressStream.print("Mode estimates with " + adaptedConfigEstimateTypes.get(i).getModeKDEType().toString().toLowerCase() + " kernel ");
            }
            progressStream.println("updated configuration file has been saved to '" + adaptedConfigFileNames.get(i) + "'.");
        }
    } // modifyTags

    /**
     * Get the value of attribute {@param attrName} from all attributes and children of {@param node}.
     * <p>
     * if the attribute does not exist, null is returned.
     *
     * @param node     tag node object
     * @param attrName attribute name
     * @return value of attribute; if the attribute does not exist, null is returned.
     */
    private String getAttribute(
            @NotNull Node node,
            @NotNull String attrName
    ) {
        // check attributes first
        NamedNodeMap attrs = node.getAttributes();
        for (int i = 0; i < attrs.getLength(); i++) {
            Node attr = attrs.item(i);

            if (attrName.equals(attr.getNodeName()))
                return getString(attr.getNodeValue(), ID_REF_PATTERN);
        }

        // check all children and their attribute 'name'
        NodeList children = node.getChildNodes();
        for (int i = 0; i < children.getLength(); i++) {
            Node child = children.item(i);
            attrs = child.getAttributes();

            if (attrName.equals(child.getNodeName()))
                return getString(attrs.getNamedItem("id").getNodeValue(), ID_REF_PATTERN);

            if (attrs != null && attrName.equals(attrs.getNamedItem("name").getNodeValue()))
                return getString(attrs.getNamedItem("id").getNodeValue(), ID_REF_PATTERN);
        }

        return null;
    } // getBranchRateModelMeanRateID

    /**
     * If {@param value} matches {@param pattern}, the first group is returned. Otherwise {@param value} is returned.
     *
     * @param value   a string to be checked
     * @param pattern a given pattern
     * @return result
     */
    private String getString(String value, Pattern pattern) {
        final Matcher matcher = pattern.matcher(value);

        if (matcher.matches())
            return matcher.group(1);
        else
            return value;
    } // getString

    public List<String> getAdaptedConfigFileNames() {
        return this.adaptedConfigFileNames;
    } // getAdaptedConfigFileNames

    public List<EstimatesTypeCollection> getAdaptedConfigEstimateTypes() {
        return this.adaptedConfigEstimateTypes;
    } // getAdaptedConfigEstimateTypes

    public int getNumOfAdaptedConfigFiles() {
        return this.numOfAdaptedConfigFiles;
    } // getNumOfAdaptedConfigFiles

    public List<String> getVcFileNames() {
        return this.vcFileNames;
    } // getVcFileNames


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("variantcallingmodeladaptor", "");
    } // printUsage


    //***********************************************
    //*                 Main method                 *
    //***********************************************


    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // Define variables
        String configFileName = null;
        String inputTreeFileName = null;
        String estimatesFileName = null;
        String outputFileName = null;
        boolean saveDetails;
        int cellThreshold = 1;
        boolean useMeanRate = true;
        boolean useOnlyBranchLength = true;

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
            System.out.println("VariantCallingModelAdaptor");

            VariantCallingModelAdaptorDialog dialog = new VariantCallingModelAdaptorDialog(new JFrame());
            if (!dialog.showDialog("VariantCallingModelAdaptor")) {
                return;
            }

            // Get parameters
            configFileName = dialog.getConfigFileName();
            if (configFileName == null) {
                Log.err.println("No model configuration file specified!");
                return;
            }

            inputTreeFileName = dialog.getInputTreeFileName();
            if (inputTreeFileName == null) {
                Log.err.println("No input tree file specified!");
                return;
            }

            estimatesFileName = dialog.getEstimatesFileName();
            if (estimatesFileName == null) {
                Log.err.println("No cached estimates file specified!");
                return;
            }

            outputFileName = dialog.getOutputFileName();

            saveDetails = dialog.getSaveDetails();

            cellThreshold = dialog.getCellThreshold();

            useMeanRate = dialog.isUseMedianRate();

            useOnlyBranchLength = dialog.isUseBranchLengthAndRate();

            try {
                new VariantCallingModelAdaptor(
                        configFileName,
                        inputTreeFileName,
                        estimatesFileName,
                        outputFileName,
                        saveDetails,
                        cellThreshold,
                        useMeanRate,
                        useOnlyBranchLength
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
                            new Arguments.StringOption("prefix", "output_file_prefix", "specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL"),
                            new Arguments.StringOption("config", "config_file", "specifies the configuration file used to performing phylogenetic analisys (either xml or json) -> MANDATORY"),
                            new Arguments.StringOption("tree", "tree_file", "specifies the tree summarized from ScsTreeAnnotator -> MANDATORY"),
                            new Arguments.StringOption("cached", "cached_estimates_file", "specifies the cached estimates file -> MANDATORY"),
                            new Arguments.StringOption("out", "output_file", "specifies the output modified configuration file (if more than one file to be generated, '_EstimatesType' will be appended) -> OPTIONAL"),
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

            // Set prefix
            if (arguments.hasOption("prefix")) {
                System.setProperty("variant.calling.file.prefix", arguments.getStringOption("prefix").trim());
            }

            // Set configFileName
            if (arguments.hasOption("config")) {
                configFileName = arguments.getStringOption("config");
            } else {
                printUsage(arguments);
                progressStream.println("Option -config is missing.");
                System.exit(1);
            }

            // Set inputTreeFileName
            if (arguments.hasOption("tree")) {
                inputTreeFileName = arguments.getStringOption("tree");
            } else {
                printUsage(arguments);
                progressStream.println("Option -tree is missing.");
                System.exit(1);
            }

            // Set estimatesFileName
            if (arguments.hasOption("cached")) {
                estimatesFileName = arguments.getStringOption("cached");
            } else {
                printUsage(arguments);
                progressStream.println("Option -cached is missing.");
                System.exit(1);
            }

            // Set outputFileName
            if (arguments.hasOption("out")) {
                outputFileName = arguments.getStringOption("out");
            }

            // Set saveTernary
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
                new VariantCallingModelAdaptor(
                        configFileName,
                        inputTreeFileName,
                        estimatesFileName,
                        outputFileName,
                        saveDetails,
                        cellThreshold,
                        useMeanRate,
                        useOnlyBranchLength
                );
            } catch (Exception e) {
                e.printStackTrace();
                progressStream.println("Exception: " + e.getMessage());
            }

            progressStream.println();
            progressStream.println("Successful!");
        }

    } // main

}
