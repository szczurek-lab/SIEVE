package beast.evolution.variantsinfo;

import beast.core.Input;
import beast.core.variantslogger.VariantsLoggable;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.rawreadcountsmodel.seqcovmodel.SeqCovModelInterface;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.variantsinfo.vcfentry.VCFEntry;
import org.jetbrains.annotations.NotNull;

import java.io.PrintStream;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.List;

import static beast.util.TreeUtils.processMetaData;

public abstract class GenericVariantsInfoVCF extends GenericVariantsInfo.Base implements VariantsLoggable {

    public final static String META_DATA_GENOTYPES = "genotypes";
    protected final static String MISSING_VALUE_IN_VCF = ".";


    //**********************************************
    //*                   Inputs                   *
    //**********************************************

    final public Input<ScsAlignment> scsDataInput = new Input<>(
            "scsData",
            "single-cell sequence data for the beast.tree",
            Input.Validate.REQUIRED
    );

    final public Input<SeqCovModelInterface.Base> seqCovModelInput = new Input<>(
            "seqCovModel",
            "model of sequencing coverage",
            Input.Validate.REQUIRED
    );

    final public Input<TreeInterface> treeInput = new Input<>(
            "tree",
            "phylogenetic beast.tree with sequence data in the leaves",
            Input.Validate.REQUIRED
    );

    final public Input<Integer> cellThresholdInput = new Input<>(
            "cellThreshold",
            "the number of mutated cells required for a site being variant (default 1)",
            Input.Validate.OPTIONAL
    );

    final public Input<AmbiguityGenotypesStrategy> ambiguityGenotypesStrategyInput = new Input<>(
            "ambiguity",
            "what to do when a site has ambiguous variant calling genotypes",
            AmbiguityGenotypesStrategy.Default,
            AmbiguityGenotypesStrategy.toArray()
    );


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    protected int cellThreshold;

    /**
     * genotypes of all nodes of the tree for each SNV after variant calling
     * [*nodes]
     */
    protected List<Integer>[] genotypes;

    protected double[] sizeFactors;

    // [#matrices] * [#loci (all)]
    protected double[][] allelicSeqCovArray = null;
    protected double[][] allelicSeqCovRawVarArray = null;

    protected boolean logSiteWiseAllelicInfo = false;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {

        scsData = scsDataInput.get();

        if (cellThresholdInput.get() != null)
            cellThreshold = cellThresholdInput.get();
        else
            cellThreshold = 1;

        ambiguityGenotypesStrategy = ambiguityGenotypesStrategyInput.get();

        sizeFactors = seqCovModelInput.get().getSizeFactors();

        final String[] seqCovStr = seqCovModelInput.get().allelicSeqCovArrayInput.get().split(";");
        allelicSeqCovArray = new double[seqCovStr.length][];
        final String[] seqCovRawVarStr = seqCovModelInput.get().allelicSeqCovRawVarArrayInput.get().split(";");
        allelicSeqCovRawVarArray = new double[seqCovRawVarStr.length][];

        for (int i = 0; i < seqCovStr.length; i++)
            allelicSeqCovArray[i] = Arrays.stream(seqCovStr[i].split(",")).mapToDouble(Double::parseDouble).toArray();
        for (int i = 0; i < seqCovRawVarStr.length; i++)
            allelicSeqCovRawVarArray[i] = Arrays.stream(seqCovRawVarStr[i].split(",")).mapToDouble(Double::parseDouble).toArray();

        if (allelicSeqCovArray.length > 1 || allelicSeqCovArray[0].length > 1)
            logSiteWiseAllelicInfo = true;

    }

    protected void initHeader(
            @NotNull PrintStream vcfOut,
            PrintStream cellNamesOut,
            PrintStream probsOut,
            PrintStream treeOut,
            PrintStream allelicInfoOut
    ) {
        // print hearer information of vcf file
        vcfOut.println("##fileformat=VCFv4.3");
        vcfOut.println("##fileDate=" + DateTimeFormatter.ofPattern("yyyyMMdd HH:mm:ss").format(LocalDateTime.now()));
        vcfOut.println("##source=SIEVE/VariantCaller");
        vcfOut.println("##FILTER=<ID=LowQual,Description=\"Low quality\">");
        vcfOut.println("##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
        vcfOut.println("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
        vcfOut.println("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
        vcfOut.println("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">");
        vcfOut.println("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
        vcfOut.println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">");
        vcfOut.println("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">");
        vcfOut.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        vcfOut.println("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, phred-scaled likelihoods for genotypes as defined in the VCF specification; in the following order: " + substModel.getAllGenotypes(",") + "\">");
        vcfOut.println("##FORMAT=<ID=ADO,Number=1,Type=Integer,Description=\"Allelic dropout (ADO) state corresponding to the number of ADOs.\">");
        vcfOut.println("##FORMAT=<ID=ADOQ,Number=1,Type=Integer,Description=\"Allelic dropout quality\">");
        vcfOut.println("##FORMAT=<ID=ADOP,Number=D,Type=Integer,Description=\"Normalized, phred-scaled likelihoods for each possible ADO state\">");
        vcfOut.println(VCFEntry.getCombinedHeader());

        // print cell names
        if (cellNamesOut != null) {
            for (int i = 0; i < this.numOfTips; i++) {
                cellNamesOut.print(this.sortedTipNames[i]);

                if (i < this.numOfTips - 1)
                    cellNamesOut.println();
            }
        }

        if (probsOut != null)
            probsOut.println("#Genotypes order: " + substModel.getAllGenotypes(","));

        if (treeOut != null) {
            ((Tree) this.tree).init(treeOut);
            treeOut.println();
        }

        if (allelicInfoOut != null) {
            if (logSiteWiseAllelicInfo)
                allelicInfoOut.println("#for a specific site in each line where information for each cell is tab-separated: allelic sequencing coverage of a cell,allelic sequencing coverage raw variance of a cell");
            else
                allelicInfoOut.println("#for a specific cell in each line where information is tab-separated: allelic sequencing coverage\tallelic sequencing coverage raw variance");
        }
    } // initHeader

    public void logTree(PrintStream out) {
        // update meta data string and length meta data string
        processMetaData(this.tree.getRoot());

        out.print("tree TREE1 = ");
        String newick = this.tree.getRoot().toSortedNewick(new int[1], true);
        out.print(newick);
        out.println(";");
        ((Tree) this.tree).close(out);
        out.println();
    }

    /**
     * Not `Loggable`.
     *
     * @param out apparently
     */
    public void init(PrintStream out) {
        throw new IllegalArgumentException("Unsupported function.");
    } // init

    /**
     * Not `Loggable`.
     *
     * @param out apparently
     */
    public void log(PrintStream out, String siteSeparator, String cellSeparator) {
        throw new IllegalArgumentException("Unsupported function.");
    } // log

    /**
     * Not `Loggable`.
     *
     * @param out apparently
     */
    public void cloze(PrintStream out) {
        throw new IllegalArgumentException("Unsupported function.");
    } // close

}
