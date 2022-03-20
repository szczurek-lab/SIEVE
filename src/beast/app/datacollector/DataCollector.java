package beast.app.datacollector;

import beast.app.tools.LogCombiner;
import beast.app.util.Arguments;
import beast.app.util.Utils;
import beast.app.utils.ChromosomeLabel;
import beast.core.Description;
import beast.core.util.Log;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.alignment.ScsBackgroundInfo;
import beast.evolution.alignment.ScsLociInfo;
import beast.evolution.alignment.VariantSiteInfo;
import beast.evolution.alignment.sequence.CovSupSeq;
import beast.evolution.alignment.sequence.FullSupsCovSeq;
import beast.evolution.branchratemodel.ScsRandomLocalClockModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.math.util.MathFunctions;
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
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Description("Collect data from a tsv document and write to a configuration xml document")
public class DataCollector {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    // constant strings
    private final static String LOCIINFO_CLASS = ScsLociInfo.class.getName();
    private final static String BACKGROUND_CLASS = ScsBackgroundInfo.class.getName();
    private final static String ALIGNMENT_CLASS = ScsAlignment.class.getName();
    private final static String STRICT_MOLECULAR_CLOCK_MODEL_CLASS = StrictClockModel.class.getName();
    private final static String RELAXED_MOLECULAR_CLOCK_MODEL_CLASS = UCRelaxedClockModel.class.getName();
    private final static String RANDOM_LOCAL_CLOCK_MODEL_CLASS = ScsRandomLocalClockModel.class.getName();

    private final static Pattern ID_REF_PATTERN = Pattern.compile("@(.+)");

    private final static List<String> RELAXED_MOLECULAR_CLOCK_MODEL_ATTRS = Stream.of("rates", "rateCategories").collect(Collectors.toList());
    private final static List<String> RANDOM_LOCAL_CLOCK_MODEL_ATTRS = Stream.of("rates", "indicators").collect(Collectors.toList());

    // Attributes of data tag which must be updated based on setups
    private final static List<String> UPDATED_ATTRS = Stream.of("dataType").collect(Collectors.toList());

    private static String MUTATIONS_CLASS;
    private static String DATATYPE;
    private static String[] BACKGROUND_NAMES;

    private final static String nrOfSamplesFlag = "=numSamples=";
    private final static String nrOfMutatedSitesFlag = "=numCandidateMutatedSites=";
    private final static String nrOfBackgroundSitesFlag = "=numBackgroundSites=";
    private final static String mutationsFlag = "=mutations=";
    private final static String backgroundFlag = "=background=";

    // structures to store read data
    protected List<String>[] dataMutations;
    protected List<String> dataMutationsInfo;
    protected List<String> dataBackground;

    // variables used to load data
    protected int nrOfCells;
    protected List<String> cellNames;
    protected List<String> excludedCellNames;

    public static PrintStream progressStream = Log.err;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    public DataCollector(
            final String cellNamesFileName,
            final String excludedCellNamesFileName,
            boolean compatibleWithSciPhi,
            boolean useSex,
            double missingDataThreshold,
            final String filteredCellNamesFileName,
            final String dataFileName,
            final DataCollectorDialog.DataType datatype,
            final String templateFileName,
            final String outputFileName,
            final String outputBaseName,
            final int[] sample,
            final int[] backgroundNameOrder
    ) throws RuntimeException {

        // 1. Parse cell names
        try {
            parseCellNames(cellNamesFileName, excludedCellNamesFileName, compatibleWithSciPhi);
        } catch (Exception e) {
            e.printStackTrace();
            progressStream.println("Error reading cell names from '" + cellNamesFileName + "'");
            throw new RuntimeException(e);
        }

        // 2. Read data from dataFileName
        try {
            readData(dataFileName, datatype, useSex, missingDataThreshold, sample);
        } catch (Exception e) {
            e.printStackTrace();
            progressStream.println("Error reading data from '" + dataFileName + "'");
            throw new RuntimeException(e);
        }

        // 3. Process tags
        try {
            processTemplateDoc(
                    filteredCellNamesFileName,
                    templateFileName,
                    outputFileName,
                    backgroundNameOrder,
                    outputBaseName
            );
        } catch (Exception e) {
            e.printStackTrace();
            progressStream.println("Error processing the template configuration file '" + templateFileName + "'");
            throw new RuntimeException(e);
        }

    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public void parseCellNames(
            @NotNull final String cellNamesFileName,
            final String excludedCellNamesFileName,
            boolean compatibleWithSciPhi
    ) throws IOException, IllegalArgumentException {
        BufferedReader fin = new BufferedReader(new FileReader(new File(cellNamesFileName)));

        cellNames = new ArrayList<>();

        while (fin.ready()) {
            String[] tmp;

            if (compatibleWithSciPhi) {
                tmp = fin.readLine().trim().split("\t");

                if (tmp.length == 2) {
                    if (tmp[1].trim().equals("CT")) {
                        String[] tmp2 = tmp[0].trim().split("/");
                        cellNames.add(tmp2[tmp2.length - 1].trim().split("\\.")[0].trim());
                    }
                } else {
                    throw new IOException("Error! The input cell names file is incompatible with SCIPhI although it is claimed so.");
                }
            } else {
                tmp = fin.readLine().trim().split("\\s+");

                if (tmp.length > 0) {
                    cellNames.addAll(Arrays.asList(tmp));
                }
            }

        }

        fin.close();

        if (excludedCellNamesFileName != null) {

            fin = new BufferedReader(new FileReader(excludedCellNamesFileName));

            excludedCellNames = new ArrayList<>();

            while (fin.ready()) {
                String[] tmp = fin.readLine().trim().split("\\s+");

                if (tmp.length > 0) {
                    for (String i : tmp) {
                        if (cellNames.contains(i)) {
                            excludedCellNames.add(i);
                        } else {
                            throw new IllegalArgumentException(cellNamesFileName + "does not contain a name: " + i +
                                    " appearing in " + excludedCellNamesFileName);
                        }
                    }
                }
            }

            fin.close();

        }

    } // parseCellNames

    /**
     * load data from .tsv document into memory
     */
    public void readData(
            final String dataFileName,
            final DataCollectorDialog.DataType datatype,
            boolean useSex,
            double missingDataThreshold,
            final int[] sample
    ) throws IOException, RuntimeException {
        List<String[]> mutationLines = new ArrayList<>();

        int numMutatedSites = 0;
        long numBackgroundSites = 0;
        int mutatedSitesCounter = 0;

        BufferedReader fin = new BufferedReader(new FileReader(dataFileName));
        String str;

        int flag = -1;
        while (fin.ready()) {
            str = fin.readLine().trim();

            // meet numSamples flag
            if (str.equals(nrOfSamplesFlag)) {
                flag = 0;
                continue;
            }

            // meet numCandidateMutatedSites flag
            if (str.equals(nrOfMutatedSitesFlag)) {
                flag = 1;
                continue;
            }

            // meet numBackgroundSites flag
            if (str.equals(nrOfBackgroundSitesFlag)) {
                flag = 2;
                continue;
            }

            // meet mutations flag
            if (str.equals(mutationsFlag)) {
                flag = 3;
                continue;
            }

            // meet background flag
            if (str.equals(backgroundFlag)) {
                if (sample[1] == -1 && numMutatedSites != 0 && mutatedSitesCounter != numMutatedSites)
                    throw new IOException("Error! Unmatched number of candidate mutated sites detected. " +
                            numMutatedSites + " specified, " + mutatedSitesCounter + " counted under " +
                            mutationsFlag + " section.");

                // remove mutation sites to meet the threshold requirement of missing data percentage
                if (missingDataThreshold < 1)
                    removeMutationLines(
                            mutationLines,
                            missingDataThreshold
                    );

                // process mutations
                sortAndParseMutationLines(
                        mutationLines,
                        datatype,
                        sample[0]
                );

                flag = 4;
                continue;
            }

            if (flag == 0) {
                // starting with numOfCells

                nrOfCells = Integer.parseInt(str);

                if (cellNames.size() != nrOfCells)
                    throw new RuntimeException("The number of cells and the number of cell names provided do not " +
                            "match: " + nrOfCells + " cells specified, but " + cellNames.size() + " cell " +
                            "names provided. (" + this.getClass().getName() + ")");

                flag = -1;

            } else if (flag == 1) {

                numMutatedSites = Integer.parseInt(str);
                flag = -1;

            } else if (flag == 2) {

                numBackgroundSites = Long.parseLong(str);
                flag = -1;

            } else if (flag == 3) {
                mutatedSitesCounter++;

                if (sample[1] > -1 && mutatedSitesCounter >= sample[1])
                    continue;

                final String[] comp = str.split("\t");

                final ChromosomeLabel chr = new ChromosomeLabel(comp[0].trim());
                if (chr.isOnMitochondrial() || (!useSex && chr.isOnSex()))
                    continue;

                mutationLines.add(comp);

            } else if (flag == 4) {

                long numCellTimesNumBGSites = 0;

                if (dataBackground == null)
                    dataBackground = new ArrayList<>();

                String[] fullStr = str.split("\t");

                if (datatype == DataCollectorDialog.DataType.CovSup) {
                    if (fullStr.length % 2 == 0) {
                        String[] partStr = new String[fullStr.length / 2];

                        for (int i = 0; i < fullStr.length; i += 2) {
                            numCellTimesNumBGSites += Long.parseLong(fullStr[i + 1].trim());

                            partStr[i / 2] = fullStr[i].trim() + "," + fullStr[i + 1].trim();
                        }

                        dataBackground.add(String.join(";", partStr));
                    } else
                        throw new RuntimeException("Background data length is not in an even number!");
                }

                if (datatype == DataCollectorDialog.DataType.FullSupsCov) {
                    for (String tmp : fullStr) {
                        tmp = tmp.trim();

                        numCellTimesNumBGSites += Long.parseLong(tmp.split(",")[1].trim());

                    }

                    dataBackground.add(String.join(";", fullStr));
                }

                if (numBackgroundSites == 0)
                    numBackgroundSites = numCellTimesNumBGSites / nrOfCells;
                else if (numBackgroundSites != numCellTimesNumBGSites / nrOfCells)
                    throw new RuntimeException("Please make sure the number of background sites are consistent in " +
                            "every line of background information.");

            }

        }

        fin.close();
    } // readData

    private void initializeVariables(
            int numCellsPerLine,
            int sampleCell
    ) {
        // check whether the size of data match the number of samples or not
        if (nrOfCells != 0) {
            if (numCellsPerLine != nrOfCells)
                throw new RuntimeException("The size of data and the number of samples do not match!");
        } else {
            nrOfCells = numCellsPerLine;
        }

        // allocate memory for dataMutations
        if (dataMutations == null) {
            dataMutations = new List[sampleCell > -1 ? Math.min(sampleCell, nrOfCells) : nrOfCells];

            for (int i = 0; i < dataMutations.length; i++)
                dataMutations[i] = new ArrayList<>();
        }
    } // initializeVariables

    private void processCovSupDatatype(
            @NotNull final String[] mutationLine,
            final int sampleCell
    ) {
        // process mutations, (coverage, alternative reads)
        if ((mutationLine.length - 4) % 2 == 0) {

            initializeVariables((mutationLine.length - 4) / 2, sampleCell);

            // load data of necessary length
            for (int i = 4;
                 i < (sampleCell > -1
                         ? Math.min(4 + 2 * sampleCell, mutationLine.length)
                         : mutationLine.length);
                 i += 2) {
                dataMutations[(i - 4) / 2].add(mutationLine[i].trim() + "," + mutationLine[i + 1].trim());
            }
        } else {
            throw new RuntimeException("Mutation data length is not an even number!");
        }
    } // processCovSupDatatype

    private void processFullSupsCovDataType(
            @NotNull final String[] mutationLine,
            final int sampleCell
    ) {
        initializeVariables(mutationLine.length - 4, sampleCell);

        // load data of necessary length
        for (int i = 4;
             i < (sampleCell > -1
                     ? Math.min(4 + sampleCell, mutationLine.length)
                     : mutationLine.length);
             i++) {
            dataMutations[i - 4].add(String.join(",", mutationLine[i].trim().split(";")));
        }
    } // processFullSupsCovDataType

    private void removeMutationLines(
            @NotNull List<String[]> mutationLines,
            final double missingDataThreshold
    ) {
        if (missingDataThreshold <= 0 || missingDataThreshold >= 1) return;

        int[] indicesToKeptCells = new int[cellNames.size() - (excludedCellNames == null ? 0 : excludedCellNames.size())];
        int j = 0;
        for (int i = 0; i < cellNames.size(); i++) {
            if (excludedCellNames == null || !excludedCellNames.contains(cellNames.get(i))) {
                indicesToKeptCells[j] = i;
                j++;
            }
        }

        MissingDataComparator comp = new MissingDataComparator(indicesToKeptCells);

        mutationLines.sort(comp);

        int entryNum = comp.getEntryNum();
        int missingEntryNum = comp.getMissingEntryNum();

        double percentage = ((double) missingEntryNum) / entryNum;
        progressStream.printf("Percentage of missing data (with some cells excluded, if specified) in the original data = %.3f\n", percentage);
        if (percentage <= missingDataThreshold) {
            progressStream.println("The specified threshold of missing data percentage (" + missingDataThreshold + ") is already met. Skip filtering.");
            return;
        }

        int numSites = mutationLines.size();
        progressStream.println("The number of sites in the original data = " + numSites);

        Iterator<String[]> iter = mutationLines.iterator();
        while (iter.hasNext()) {
            final String[] line = iter.next();
            iter.remove();
            numSites--;

            final int[] tmp = comp.getEntryNum(line);

            missingEntryNum -= tmp[0];
            entryNum -= tmp[1];

            percentage = ((double) missingEntryNum) / entryNum;
            if (percentage <= missingDataThreshold) {
                progressStream.printf("Current percentage of missing data = %.3f\n", percentage);
                progressStream.println("Current number of sites = " + numSites);
                break;
            }
        }
    } // removeMutationLines

    private void sortAndParseMutationLines(
            @NotNull List<String[]> mutationLines,
            final DataCollectorDialog.DataType datatype,
            final int sampleCell
    ) {
        mutationLines.sort(new LociInfoComparator());

        for (String[] fullStr : mutationLines) {
            // process mutation information, in default (chromosome, location, ref, alts)
            if (dataMutationsInfo == null)
                dataMutationsInfo = new ArrayList<>();

            StringBuilder infoStr = new StringBuilder(fullStr[0].trim());
            for (int i = 1; i < 4; i++)
                infoStr.append(",").append(fullStr[i].trim());
            dataMutationsInfo.add(infoStr.toString());

            // process mutations, (coverage,alternative reads)
            if (datatype == DataCollectorDialog.DataType.CovSup)
                processCovSupDatatype(fullStr, sampleCell);

            // process mutations, (variant1,variant2,variant3;reads1,reads2,reads3,coverage)
            if (datatype == DataCollectorDialog.DataType.FullSupsCov)
                processFullSupsCovDataType(fullStr, sampleCell);
        }
    } // sortAndParseMutationLines

    /**
     * get the default attributes of data label
     *
     * @return a map containing (key, value) pairs
     */
    public HashMap<String, String> getDataDefaultAttrs() {
        HashMap<String, String> attrs = new HashMap<>();
        attrs.put("id", "data");
        attrs.put("spec", ALIGNMENT_CLASS);
        attrs.put("dataType", DATATYPE);
        attrs.put("ascertained", "none");
        attrs.put("bgInfoSanityCheck", "true");
        attrs.put("bgSitesNum", "-1");
        attrs.put("meanAscBiasCorrection", "mean");
        return attrs;
    } // getDataDefaultAttrs

    /**
     * make lociInfo label
     *
     * @param doc parsed object of source .xml
     * @return a lociInfo element
     */
    public Element makeLociInfoNode(@NotNull Document doc) {
        Element lociInfoElem = doc.createElement("lociInfo");
        lociInfoElem.setAttribute("spec", LOCIINFO_CLASS);
        lociInfoElem.setAttribute("loci", "true");
        lociInfoElem.setTextContent(String.join(";", dataMutationsInfo));
        return lociInfoElem;
    } // makeLociInfoNode

    /**
     * make sequence label
     *
     * @param doc parsed object of source .xml
     * @return an array of sequence element
     */
    public Element[] makeMutationsNode(
            final String filteredCellNamesFileName,
            Document doc
    ) {
        Element[] mutationsElem;
        if (excludedCellNames != null) {
            mutationsElem = new Element[cellNames.size() - excludedCellNames.size()];
        } else {
            mutationsElem = new Element[cellNames.size()];
        }

        List<String> filteredCellNames = new ArrayList<>();

        List<String> cellNamesCopy = new ArrayList<>(cellNames);

        int[] indices = new int[this.nrOfCells];
        Arrays.fill(indices, -1);

        cellNamesCopy.sort((s, t1) -> {
            if (s.compareTo(t1) > 0) {
                return 1;
            } else if (s.compareTo(t1) < 0) {
                return -1;
            }

            throw new RuntimeException("Duplicate cell names found: " + s + ", " + t1);
        });

        for (int i = 0; i < nrOfCells; i++) {
            for (int j = 0; j < nrOfCells; j++) {
                if (cellNamesCopy.get(i).equals(cellNames.get(j))) {
                    indices[i] = j;
                }
            }
        }

        int j = 0;
        for (int i = 0; i < nrOfCells; i++) {
            if (excludedCellNames == null || !excludedCellNames.contains(cellNamesCopy.get(i))) {
                filteredCellNames.add(cellNamesCopy.get(i));

                mutationsElem[j] = doc.createElement("sequence");
                mutationsElem[j].setAttribute("spec", MUTATIONS_CLASS);
                mutationsElem[j].setAttribute("taxon", cellNamesCopy.get(i));
                mutationsElem[j].setTextContent(String.join(";", dataMutations[indices[i]]));
                j++;
            }
        }

        if (filteredCellNamesFileName != null) {
            PrintStream out;
            try {
                out = new PrintStream(filteredCellNamesFileName);

                for (int i = 0; i < filteredCellNames.size(); i++) {
                    out.print(filteredCellNames.get(i));

                    if (i < filteredCellNames.size() - 1) {
                        out.println();
                    }
                }

                out.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        return mutationsElem;
    } // makeMutationsNode

    /**
     * make backgroundInfo label
     *
     * @param doc                 parsed object of source .xml
     * @param backgroundNameOrder
     * @return an array of backgroundInfo element
     */
    public Element[] makeBackgroundNode(Document doc, final int[] backgroundNameOrder) {
        Element[] backgroundElem = new Element[backgroundNameOrder.length];
        for (int i = 0; i < backgroundNameOrder.length; i++) {
            backgroundElem[i] = doc.createElement("backgroundInfo");
            backgroundElem[i].setAttribute("spec", BACKGROUND_CLASS);
            backgroundElem[i].setAttribute("backgroundName", BACKGROUND_NAMES[backgroundNameOrder[i]]);
            backgroundElem[i].setTextContent(String.join(";", dataBackground.get(i)));
        }
        return backgroundElem;
    } // makeBackgroundNode

    /**
     * modify the data node
     *
     * @param filteredCellNamesFileName a file storing filtered cell names
     * @param doc                       parsed object of source .xml
     * @param dataNode                  a data node to be modified
     * @param backgroundNameOrder       apparently
     */
    private void buildDataNode(
            final String filteredCellNamesFileName,
            Document doc,
            @NotNull Element dataNode,
            final int[] backgroundNameOrder
    ) {
        Map<String, String> dataDefaultAttrs = getDataDefaultAttrs();

        // process attributes
        if (dataNode.getAttributes().getLength() == 0) {
            for (String key : dataDefaultAttrs.keySet()) {
                dataNode.setAttribute(key, dataDefaultAttrs.get(key));
            }
        } else {
            List<String> uselessAttrs = new ArrayList<>();
            List<String> processedAttrs = new ArrayList<>();

            // replace some existing attributes (not "id")
            for (int i = 0; i < dataNode.getAttributes().getLength(); i++) {
                final String attr = dataNode.getAttributes().item(i).getNodeName();
                if (dataDefaultAttrs.containsKey(attr)) {
                    processedAttrs.add(attr);

                    if (UPDATED_ATTRS.contains(attr))
                        dataNode.setAttribute(attr, dataDefaultAttrs.get(attr));
                } else
                    uselessAttrs.add(attr);
            }

            // add missing attributes
            if (processedAttrs.size() != dataDefaultAttrs.keySet().size()) {
                for (String key : dataDefaultAttrs.keySet()) {
                    if (!processedAttrs.contains(key))
                        dataNode.setAttribute(key, dataDefaultAttrs.get(key));
                }
            }

            // remove useless attributes
            if (uselessAttrs.size() > 0) {
                for (String i : uselessAttrs) {
                    dataNode.removeAttribute(i);
                }
            }
        }

        // process its children
        // first, cleanup all the existing children
        if (dataNode.getChildNodes().getLength() > 0) {
            int size = dataNode.getChildNodes().getLength();
            for (int i = 0; i < size; i++) {
                dataNode.removeChild(dataNode.getChildNodes().item(0));
            }
        }

        // second, add default children
        // lociInfo
        if (dataMutationsInfo != null)
            dataNode.appendChild(makeLociInfoNode(doc));

        // sequences
        for (Element i : makeMutationsNode(filteredCellNamesFileName, doc)) {
            dataNode.appendChild(i);
        }

        // backgroundInfo
        if (dataBackground != null) {
            for (Element i : makeBackgroundNode(doc, backgroundNameOrder)) {
                dataNode.appendChild(i);
            }
        }
    } // buildDataNode

    /**
     * change the value of "fileName" attribute of logger tags to the base name of input template configuration xml document
     *
     * @param loggerNodes    logger tags
     * @param outputBaseName base name of output file
     */
    private void processLogs(
            List<Node> loggerNodes,
            String outputBaseName
    ) {
        for (Node loggerNode : loggerNodes) {
            NamedNodeMap attributes = loggerNode.getAttributes();
            for (int i = 0; i < attributes.getLength(); i++) {
                Node attr = attributes.item(i);

                if (attr.getNodeName().equals("fileName")) {
                    String[] subStrs = attr.getNodeValue().split("\\.");
                    if (subStrs.length > 1) {
                        attr.setNodeValue(String.join(".", outputBaseName, subStrs[subStrs.length - 1]));
                    } else {
                        attr.setNodeValue(outputBaseName);
                    }
                }
            }
        }
    } // processLogs

    /**
     * update certain nodes if not using strict molecular clock
     *
     * @param doc xml object
     */
    private void updateNonStrictClockModel(@NotNull Document doc) {
        NodeList branchRateModelNodes = doc.getElementsByTagName("branchRateModel");

        if (branchRateModelNodes == null || branchRateModelNodes.getLength() == 0) {
            progressStream.println("No 'branchRateModel' tags found.");
            return;
        } else if (branchRateModelNodes.getLength() > 1) {
            progressStream.println("More than one 'branchRateModel' tags found, while only one is expected.");
            return;
        }

        Node branchRateModelNode = branchRateModelNodes.item(0);
        NamedNodeMap branchRateModelAttributes = branchRateModelNode.getAttributes();
        Node specAttr = branchRateModelAttributes.getNamedItem("spec");

        final String[] branchRateModelName = specAttr.getNodeValue().split("\\.");

        // strict molecular clock model
        if (STRICT_MOLECULAR_CLOCK_MODEL_CLASS.equals(specAttr.getNodeValue()) ||
                STRICT_MOLECULAR_CLOCK_MODEL_CLASS.contains(branchRateModelName[branchRateModelName.length - 1]))
            return;

        // relaxed molecular clock model
        if (RELAXED_MOLECULAR_CLOCK_MODEL_CLASS.equals(specAttr.getNodeValue()) ||
                RELAXED_MOLECULAR_CLOCK_MODEL_CLASS.contains(branchRateModelName[branchRateModelName.length - 1]))
            updateNonStrictClockModelParamDim(doc, branchRateModelAttributes, RELAXED_MOLECULAR_CLOCK_MODEL_ATTRS);

        // random local clock model
        if (RANDOM_LOCAL_CLOCK_MODEL_CLASS.equals(specAttr.getNodeValue()) ||
                RANDOM_LOCAL_CLOCK_MODEL_CLASS.contains(branchRateModelName[branchRateModelName.length - 1]))
            updateNonStrictClockModelParamDim(doc, branchRateModelAttributes, RANDOM_LOCAL_CLOCK_MODEL_ATTRS);
    } // updateNonStrictClockModel

    /**
     * update parameter dimension based on the number of tips.
     *
     * @param doc                  xml object
     * @param branchRateModelAttrs attributes of non-strict branch rate model
     * @param attrs                names of attributes to be updated
     */
    private void updateNonStrictClockModelParamDim(
            @NotNull Document doc,
            @NotNull NamedNodeMap branchRateModelAttrs,
            @NotNull final List<String> attrs
    ) {
        // collect state nodes
        NodeList stateTags = doc.getElementsByTagName("state");
        if (stateTags == null || stateTags.getLength() == 0) {
            progressStream.println("No 'state' tags found.");
            return;
        } else if (stateTags.getLength() > 1) {
            progressStream.println("More than one 'state' tags found, while only one is expected.");
            return;
        }
        NodeList stateNodes = stateTags.item(0).getChildNodes();

        for (int i = 0; i < branchRateModelAttrs.getLength(); i++) {
            final Node attr = branchRateModelAttrs.item(i);

            if (attrs.contains(attr.getNodeName())) {
                final Matcher matcher = ID_REF_PATTERN.matcher(attr.getNodeValue());
                if (matcher.matches())
                    updateParamDim(stateNodes, matcher.group(1));
                else
                    throw new IllegalArgumentException("Error! Invalid format for the value of " + attr.getNodeName());
            }
        }
    } // updateRelaxedMolecularClockModel

    private void updateParamDim(
            @NotNull NodeList stateNodes,
            @NotNull String id
    ) {
        boolean found = false;
        for (int i = 0; i < stateNodes.getLength(); i++) {
            NamedNodeMap attrs = stateNodes.item(i).getAttributes();

            if (attrs != null && attrs.getNamedItem("id") != null && attrs.getNamedItem("id").getNodeValue().equals(id)) {
                if (found)
                    throw new IllegalArgumentException("Error! Duplicate state node found with id: " + id);

                found = true;

                Node dimAttr = attrs.getNamedItem("dimension");

                final int excludedCellNr = excludedCellNames == null ? 0 : excludedCellNames.size();
                final String dim = String.valueOf(2 * (cellNames.size() - excludedCellNr) - 1);

                if (dimAttr == null)
                    ((Element) stateNodes.item(i)).setAttribute("dimension", dim);
                else
                    dimAttr.setNodeValue(dim);
            }
        }
    } // updateParamDim

    /**
     * load source .xml document, modify data node, and output to a new .xml document
     */
    public void processTemplateDoc(
            final String filteredCellNamesFileName,
            final String templateFileName,
            final String outputFileName,
            final int[] backgroundNameOrder,
            String outputBaseName
    ) throws TransformerException, ParserConfigurationException, IOException, SAXException {
        DocumentBuilderFactory dbFac = DocumentBuilderFactory.newInstance();
        dbFac.setIgnoringElementContentWhitespace(true);
        Document doc;

        doc = dbFac.newDocumentBuilder().parse(new File(templateFileName));
        doc.getDocumentElement().normalize();

        // filter "data" element whose parent node is "beast"
        NodeList candidateDataNodes = doc.getElementsByTagName("data");
        if (candidateDataNodes == null || candidateDataNodes.getLength() == 0) {
            throw new RuntimeException("No tag named 'data' is found! Please check your input template. (" +
                    this.getClass().getName() + ")");
        }
        Node dataNode = null;
        for (int i = 0; i < candidateDataNodes.getLength(); i++) {
            if (candidateDataNodes.item(i).getParentNode().getNodeName().equals("beast")) {
                dataNode = candidateDataNodes.item(i);
                break;
            }
        }

        // update data element
        if (dataNode == null) {
            // no data element exist; create a new one

            Element newDataNode = doc.createElement("data");
            Node root;

            if (doc.getElementsByTagName("beast").getLength() == 1) {
                root = doc.getElementsByTagName("beast").item(0);
            } else {
                throw new RuntimeException("Only one 'beast' label is supposed to exist, but " + doc.getElementsByTagName("beast").getLength() + " found.");
            }

            buildDataNode(filteredCellNamesFileName, doc, newDataNode, backgroundNameOrder);
            root.appendChild(newDataNode);

        } else
            buildDataNode(filteredCellNamesFileName, doc, (Element) dataNode, backgroundNameOrder);

        // filter "logger" elements whose parent node is "run" and have an attribute of "spec" equivalent to "Logger"
        NodeList candidateLoggerNodes = doc.getElementsByTagName("logger");
        if (candidateLoggerNodes == null || candidateLoggerNodes.getLength() == 0) {
            throw new RuntimeException("No tags named 'logger' are found! Please check your input template. (" +
                    this.getClass().getName() + ")");
        }
        List<Node> loggerNodes = new ArrayList<>();
        for (int i = 0; i < candidateLoggerNodes.getLength(); i++) {
            if (candidateLoggerNodes.item(i).getParentNode().getNodeName().equals("run")) {
                for (int j = 0; j < candidateLoggerNodes.item(i).getAttributes().getLength(); j++) {
                    if (candidateLoggerNodes.item(i).getAttributes().item(j).getNodeName().equals("spec") &&
                            candidateLoggerNodes.item(i).getAttributes().item(j).getNodeValue().equals("Logger")) {
                        loggerNodes.add(candidateLoggerNodes.item(i));
                    }
                }
            }
        }
        processLogs(loggerNodes, outputBaseName);

        // update certain nodes if not using strict molecular clock
        updateNonStrictClockModel(doc);

        // save the xml document to outputTarget
        Transformer tf = TransformerFactory.newInstance().newTransformer();
        tf.setOutputProperty(OutputKeys.INDENT, "yes");
        DOMSource domS = new DOMSource(doc);
        File outputConfig = new File(outputFileName);

        if (outputFileName.contains("/"))
            outputConfig.getParentFile().mkdirs();

        StreamResult sr = new StreamResult(outputConfig);
        tf.transform(domS, sr);
    } // processTargetDoc


    //**********************************************
    //*               Static methods               *
    //**********************************************

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("datacollector", "");
        progressStream.println("IMPORTANT: please make sure that in the input xml template the section containing " +
                "alignment has a tag named 'data', and the log sections have a tag named 'logger'.");
    } // printUsage

    public static void setConstants(DataCollectorDialog.DataType dataType) {
        if (dataType == DataCollectorDialog.DataType.CovSup) {
            MUTATIONS_CLASS = CovSupSeq.class.getName();
            DATATYPE = "coverage-support";
            BACKGROUND_NAMES = new String[]{"coverage", "variant", "normal"};
        } else if (dataType == DataCollectorDialog.DataType.FullSupsCov) {
            MUTATIONS_CLASS = FullSupsCovSeq.class.getName();
            DATATYPE = "full supports-coverage";
            BACKGROUND_NAMES = new String[]{"variant1", "variant2", "variant3", "normal", "coverage"};
        }
    } // setConstants


    //************************************************
    //*                Nested classes                *
    //************************************************

    static class LociInfoComparator implements Comparator<String[]> {
        @Override
        public int compare(String[] s1, String[] s2) throws RuntimeException {
            if (s1.length != s2.length) {
                throw new RuntimeException("Error: Compared string arrays do not share the same length." +
                        Arrays.toString(s1) + ", " + Arrays.toString(s2));
            }

            // chromosome label
            ChromosomeLabel cl1 = new ChromosomeLabel(s1[0]);
            ChromosomeLabel cl2 = new ChromosomeLabel(s2[0]);
            final int chrResult = cl1.compareTo(cl2);
            if (chrResult != 0) {
                return chrResult;
            }

            // position
            final long pos1 = Long.parseLong(s1[1]);
            final long pos2 = Long.parseLong(s2[1]);
            if (pos1 > pos2) {
                return 1;
            } else if (pos1 < pos2) {
                return -1;
            }

            // alt nucleotide
            final int alt = s1[3].trim().toUpperCase().compareTo(s2[3].trim().toUpperCase());
            if (alt > 0) {
                return 1;
            } else if (alt < 0) {
                return -1;
            }

            // ref nucleotide
            assert s1[2].trim().equalsIgnoreCase(s2[2].trim());

            throw new RuntimeException("Duplicate sites information detected: " + s1[0] + ", " + s1[1] + ", " +
                    s1[2] + ", " + s1[3]);
        } // compare
    } // class LociInfoComparator

    static class MissingDataComparator implements Comparator<String[]> {

        private final int[] indicesToKeptCells;
        private final List<VariantSiteInfo> varInfo;

        // The number of total entries of cells to be kept across all sites.
        private int entryNum;

        // The number of entries of missing data of cells to be kept across all sites.
        private int missingEntryNum;

        MissingDataComparator(final int[] indicesToKeptCells) {
            this.indicesToKeptCells = indicesToKeptCells;
            this.varInfo = new ArrayList<>();
        } // MissingDataComparator

        public int[] getEntryNum(final String[] o) {
            int[] ret = {0, o.length - 4};

            for (int i : this.indicesToKeptCells) {
                final String[] comp = o[i + 4].split(",");
                if (Integer.parseInt(comp[comp.length - 1]) == 0) ret[0]++;
            }

            return ret;
        } // getEntryNum

        @Override
        public int compare(final String[] o1, final String[] o2) {
            if (o1.length != o2.length) {
                throw new RuntimeException("Error: Compared string arrays do not share the same length." +
                        Arrays.toString(o1) + ", " + Arrays.toString(o2));
            }

            if (this.indicesToKeptCells.length > o1.length - 4 || MathFunctions.max(Arrays.stream(this.indicesToKeptCells).boxed().toArray(Integer[]::new)) > o1.length - 4) {
                throw new RuntimeException("Error: indices exceed valid range.");
            }

            final VariantSiteInfo var1 = new VariantSiteInfo(
                    o1[0],
                    Long.parseLong(o1[1]),
                    o1[2],
                    o1[3].split(",")
            );

            final VariantSiteInfo var2 = new VariantSiteInfo(
                    o2[0],
                    Long.parseLong(o2[1]),
                    o2[2],
                    o2[3].split(",")
            );

            final int[] entryNum1 = getEntryNum(o1);
            final int[] entryNum2 = getEntryNum(o2);

            if (!this.varInfo.contains(var1)) {
                this.varInfo.add(var1);
                this.missingEntryNum += entryNum1[0];
                this.entryNum += entryNum1[1];
            }

            if (!this.varInfo.contains(var2)) {
                this.varInfo.add(var2);
                this.missingEntryNum += entryNum2[0];
                this.entryNum += entryNum2[1];
            }

            return Integer.compare(entryNum2[0], entryNum1[0]);
        } // compare

        public int getEntryNum() {
            return entryNum;
        } // getEntryNum

        public int getMissingEntryNum() {
            return missingEntryNum;
        } // getMissingEntryNum
    } // class MissingDataComparator


    //**********************************************
    //*                    Main                    *
    //**********************************************

    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        // file names
        String cellNamesFileName = null;
        String dataFileName = null;
        String templateFileName = null;
        String excludedCellNamesFileName = null;
        String filteredCellNamesFileName = null;
        String outputFileName = null;
        String cellNameBaseName; // base name of cell names file
        String rootPath;

        // whether cell names file is compatible with SciPhi or not
        boolean compatibleWithSciPhi = false;

        // whether using candidate sites from sex chromosomes or not
        boolean useSex = true;

        DataCollectorDialog.DataType datatype = DataCollectorDialog.DataType.FullSupsCov;

        // variables optionally defined
        int[] sample = {-1, -1}; // cells, loci
        int[] backgroundNameOrder = {0, 1, 2, 3, 4};
        double missingDataThreshold = 1.0;

        // No arguments provided, launch GUI
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
            System.out.println("DataCollector");

            DataCollectorDialog dialog = new DataCollectorDialog(new JFrame());
            if (!dialog.showDialog("DataCollector"))
                return;

            // Get parameters
            cellNamesFileName = dialog.getCellNamesFileName();
            if (cellNamesFileName == null) {
                Log.err.println("No cell names file specified!");
                return;
            }

            compatibleWithSciPhi = dialog.isCellNamesCompatibleWithSCIPhI();

            useSex = dialog.isUseSex();

            dataFileName = dialog.getDataFileName();
            if (dataFileName == null) {
                Log.err.println("No data file specified!");
                return;
            }

            datatype = dialog.getDataType();

            templateFileName = dialog.getTemplateFileName();
            if (templateFileName == null) {
                Log.err.println("No template configuration file specified!");
                return;
            }

            excludedCellNamesFileName = dialog.getExcludedCellNamesFile();

            outputFileName = dialog.getOutputFileName();

            sample = dialog.getSamplesControl();
            backgroundNameOrder = dialog.getBackgroundInfoOrder();

        } else {

            // TODO: print some information here
            Arguments arguments = new Arguments(
                    new Arguments.Option[]{
                            new Arguments.Option("help", "option to print this message -> OPTIONAL"),
                            new Arguments.StringOption("prefix", "output_file_prefix", "specifies the prefix of output files (a folder must be ended with '/') -> OPTIONAL"),
                            new Arguments.StringOption("cell", "cell_names_file", "specifies a blank spaces separated document containing cell names -> MANDATORY"),
                            new Arguments.Option("sciphi", "specifies whether the file of cell names is compatible with SCIPhI (default incompatible) -> OPTIONAL"),
                            new Arguments.StringOption("data", "data_file", "provides a data document ending with .tsv -> MANDATORY"),
                            new Arguments.IntegerOption("datatype", 0, 1, "specifies input datatype (0: Coverage-Support, 1: Full supports-Coverage; default: 1) -> OPTIONAl"),
                            new Arguments.Option("ignoreSex", "ignores candidate mutated sites from sex chromosomes; usually turned on for males -> OPTIONAL"),
                            new Arguments.StringOption("template", "template_configuration_file", "provides a configuration document ending with .xml -> MANDATORY"),
                            new Arguments.StringOption("exclude", "excluded_cell_names", "specifies names of the cells to be excluded -> OPTIONAL"),
                            new Arguments.StringOption("out", "output_file", "specifies the configuration file integrating with the input data -> OPTIONAL"),
                            new Arguments.IntegerArrayOption("sample", 2, -1, 100000, "samples a part of the input data for test purposes; the first number defines the number of sampled cells, and the second number defines the number of sampled loci; by default, all the data will be loaded -> OPTIONAL"),
                            new Arguments.IntegerArrayOption("bgcs", 3, 0, 2, "specifies the order of background information for \"Coverage-Support\" datatype w.r.t. (coverage 0, variant 1, normal 2); default: 0 1 2; working with \"-datatype 0\" if specified -> OPTIONAL"),
                            new Arguments.IntegerArrayOption("bgfsc", 5, 0, 4, "specifies the order of background information for \"Full support-Coverage\" datatype w.r.t. (0 - variant1, 1 - variant2, 2 - variant3, 3 - normal, 4 - coverage); default: 0 1 2 3 4; working with \"-datatype 1\" if specified -> OPTIONAL"),
                            new Arguments.RealOption("miss", 0.001, 1, "specifies the threshold of missing data percentage; set to 1 to turn off the selection; the programme will choose as many sites as possible without exceeding the specified threshold -> OPTIONAL")
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
            if (arguments.hasOption("prefix"))
                System.setProperty("data.collection.file.prefix", arguments.getStringOption("prefix").trim());

            // Set cellNamesFileName
            if (arguments.hasOption("cell"))
                cellNamesFileName = arguments.getStringOption("cell");
            else {
                printUsage(arguments);
                Log.err.println("Option -cell is missing.");
                System.exit(1);
            }

            // Set compatibleWithSciPhi
            if (arguments.hasOption("sciphi"))
                compatibleWithSciPhi = true;

            // Set useSex
            if (arguments.hasOption("ignoreSex"))
                useSex = false;

            // Set dataFileName
            if (arguments.hasOption("data"))
                dataFileName = arguments.getStringOption("data");
            else {
                printUsage(arguments);
                Log.err.println("Option -data is missing.");
                System.exit(1);
            }

            // Set datatype
            if (arguments.hasOption("datatype") && arguments.getIntegerOption("datatype") == 0)
                datatype = DataCollectorDialog.DataType.CovSup;

            // Set templateFileName
            if (arguments.hasOption("template"))
                templateFileName = arguments.getStringOption("template");
            else {
                printUsage(arguments);
                Log.err.println("Option -template is missing.");
                System.exit(1);
            }

            // Set excludedCellNamesFileName
            if (arguments.hasOption("exclude"))
                excludedCellNamesFileName = arguments.getStringOption("exclude");

            // Set outputFileName
            if (arguments.hasOption("out"))
                outputFileName = arguments.getStringOption("out");

            // Set sample
            if (arguments.hasOption("sample"))
                sample = arguments.getIntegerArrayOption("sample");

            // Set background
            if (arguments.hasOption("bgcs")) {
                if (datatype != DataCollectorDialog.DataType.CovSup) {
                    printUsage(arguments);
                    Log.err.println("Option -bgcs should be specified together with \"-datatype 0\".");
                    System.exit(1);
                }

                backgroundNameOrder = arguments.getIntegerArrayOption("bgcs");
            }

            if (arguments.hasOption("bgfsc")) {
                if (datatype != DataCollectorDialog.DataType.FullSupsCov) {
                    printUsage(arguments);
                    Log.err.println("Option -bgfsc should be specified together with \"-datatype 1\".");
                    System.exit(1);
                }

                backgroundNameOrder = arguments.getIntegerArrayOption("bgfsc");
            }

            if (arguments.hasOption("bgcs") && arguments.hasOption("bgfsc")) {
                printUsage(arguments);
                Log.err.println("Either \"-bgcs\" or \"-bgfsc\" can be specified, not both.");
                System.exit(1);
            }

            // Set missing data threshold
            if (arguments.hasOption("miss"))
                missingDataThreshold = arguments.getRealOption("miss");
        }

        cellNameBaseName = FileNameProcessor.getBaseName(cellNamesFileName);
        rootPath = FileNameProcessor.getRootPath(cellNamesFileName);

        if (outputFileName == null) {
            outputFileName = FileNameProcessor.getBaseName(templateFileName) + "_updated" + FileNameProcessor.getSuffix(templateFileName);

            if (System.getProperty("data.collection.file.prefix") != null)
                outputFileName = System.getProperty("data.collection.file.prefix") + outputFileName;
            else
                outputFileName = rootPath + outputFileName;
        } else {
            if (System.getProperty("data.collection.file.prefix") != null)
                outputFileName = System.getProperty("data.collection.file.prefix") + outputFileName;
        }

        if (excludedCellNamesFileName != null) {
            if (System.getProperty("data.collection.file.prefix") != null)
                filteredCellNamesFileName = System.getProperty("data.collection.file.prefix");
            else
                filteredCellNamesFileName = rootPath;

            filteredCellNamesFileName = filteredCellNamesFileName + "filtered_" + cellNameBaseName;

            String tmp = FileNameProcessor.getSuffix(cellNamesFileName);
            if (tmp.length() > 0)
                filteredCellNamesFileName = filteredCellNamesFileName + tmp;
        }

        setConstants(datatype);

        try {
            new DataCollector(
                    cellNamesFileName,
                    excludedCellNamesFileName,
                    compatibleWithSciPhi,
                    useSex,
                    missingDataThreshold,
                    filteredCellNamesFileName,
                    dataFileName,
                    datatype,
                    templateFileName,
                    outputFileName,
                    FileNameProcessor.getBaseName(outputFileName),
                    sample,
                    backgroundNameOrder
            );
        } catch (Exception e) {
            e.printStackTrace();
            progressStream.println("Exception: " + e.getMessage());
            return;
        }

        progressStream.println();
        progressStream.println("Successful!");

    } // main

}
