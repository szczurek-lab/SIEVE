package beast.core.variantslogger;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.util.Log;
import beast.util.FileNameProcessor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

public class VariantsLogger extends BEASTObject {


    //**********************************************
    //*                   Inputs                   *
    //**********************************************

    final public Input<String> fileNameInput = new Input<>("fileName", "Name of the vcf file",
            Input.Validate.REQUIRED);

    final public Input<VariantsLoggable> loggerInput = new Input<>("log",
            "Element in a log. This can be any plug in that is VariantsLoggable.",
            Input.Validate.REQUIRED, VariantsLoggable.class);

    final public Input<Boolean> saveDetailsInput = new Input<>(
            "saveDetails",
            "whether to save detailed information, such as inferred maximum likelihood ADO states (*.ado), " +
                    "inferred maximum likelihood genotypes at tips (*.genotypes) and their corresponding ternary " +
                    "matrix (*.ternary) and probabilities (*.probs), with cell names (*.cell_names) and loci " +
                    "information (*.loci_info) in separate files (default false). An annotated intermediate tree " +
                    "(*.intermediate_tree) with maximum likelihood genotypes of all nodes labeled will also be " +
                    "saved for gene annotation purpose in the following step. Besides, for each cell, size factor " +
                    "(*.size_factors) as well as allelic sequencing coverage and raw variance (*.allelic_info) " +
                    "will also be saved.",
            Input.Validate.OPTIONAL
    );


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private String vcFileName;
    private String adoFileName = null;
    private String genotypesFileName = null;
    private String ternaryFileName = null;
    private String cellNamesFileName = null;
    private String lociInfoFileName = null;
    private String probsFileName = null;
    private String treeFileName = null;
    private String allelicInfoFileName = null;
    private String sizeFactorFileName = null;


    private VariantsLoggable logger;

    private boolean saveDetails;

    protected PrintStream vcfOut;
    protected PrintStream adoOut;
    protected PrintStream genotypesOut;
    protected PrintStream ternaryOut;
    protected PrintStream cellNamesOut;
    protected PrintStream lociInfoOut;
    protected PrintStream probsOut;
    protected PrintStream treeOut;
    protected PrintStream allelicInfoOut;
    protected PrintStream sizeFactorOut;



    //**********************************************
    //*                Constructors                *
    //**********************************************

    public VariantsLogger() {
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {

        this.vcFileName = fileNameInput.get();

        this.logger = loggerInput.get();

        if (saveDetailsInput.get() != null)
            this.saveDetails = saveDetailsInput.get();
        else
            this.saveDetails = false;

        if (this.saveDetails) {
            String prefix = FileNameProcessor.getRootPath(this.vcFileName) +
                    FileNameProcessor.getBaseName(this.vcFileName);

            this.adoFileName = prefix + ".ado";
            this.genotypesFileName = prefix + ".genotypes";
            this.ternaryFileName = prefix + ".ternary";
            this.cellNamesFileName = prefix + ".cell_names";
            this.lociInfoFileName = prefix + ".loci_info";
            this.probsFileName = prefix + ".probs";
            this.treeFileName = prefix + ".intermediate_tree";
            this.allelicInfoFileName = prefix + ".allelic_info";
            this.sizeFactorFileName = prefix + ".size_factors";

        }

    } // initAndValidate

    protected void openFiles() throws FileNotFoundException {
        if (System.getProperty("variant.calling.file.prefix") != null)
            this.vcFileName = System.getProperty("variant.calling.file.prefix") + this.vcFileName;
        final File vcFile = new File(this.vcFileName);
        if (vcFile.exists())
            Log.err.println("Overwriting " + this.vcFileName + "...");
        if (this.vcFileName.contains("/"))
            vcFile.getParentFile().mkdirs();
        this.vcfOut = new PrintStream(vcFile);

        if (this.saveDetails) {
            if (System.getProperty("variant.calling.file.prefix") != null)
                this.adoFileName = System.getProperty("variant.calling.file.prefix") + this.adoFileName;
            final File adoFile = new File(this.adoFileName);
            if (adoFile.exists())
                Log.err.println("Overwriting " + this.adoFileName + "...");
            this.adoOut = new PrintStream(this.adoFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.genotypesFileName = System.getProperty("variant.calling.file.prefix") + this.genotypesFileName;
            final File genotypesFile = new File(this.genotypesFileName);
            if (genotypesFile.exists())
                Log.err.println("Overwriting " + this.genotypesFileName + "...");
            this.genotypesOut = new PrintStream(this.genotypesFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.ternaryFileName = System.getProperty("variant.calling.file.prefix") + this.ternaryFileName;
            final File ternaryFile = new File(this.ternaryFileName);
            if (ternaryFile.exists())
                Log.err.println("Overwriting " + this.ternaryFileName + "...");
            this.ternaryOut = new PrintStream(this.ternaryFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.cellNamesFileName = System.getProperty("variant.calling.file.prefix") + this.cellNamesFileName;
            final File cellNamesFile = new File(this.cellNamesFileName);
            if (cellNamesFile.exists())
                Log.err.println("Overwriting " + this.cellNamesFileName + "...");
            this.cellNamesOut = new PrintStream(this.cellNamesFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.lociInfoFileName = System.getProperty("variant.calling.file.prefix") + this.lociInfoFileName;
            final File lociInfoFile = new File(this.lociInfoFileName);
            if (lociInfoFile.exists())
                Log.err.println("Overwriting " + this.lociInfoFileName + "...");
            this.lociInfoOut = new PrintStream(this.lociInfoFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.probsFileName = System.getProperty("variant.calling.file.prefix") + this.probsFileName;
            final File probsFile = new File(this.probsFileName);
            if (probsFile.exists())
                Log.err.println("Overwriting " + this.probsFileName + "...");
            this.probsOut = new PrintStream(this.probsFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.treeFileName = System.getProperty("variant.calling.file.prefix") + this.treeFileName;
            final File treeFile = new File(this.treeFileName);
            if (treeFile.exists())
                Log.err.println("Overwriting " + this.treeFileName + "...");
            this.treeOut = new PrintStream(this.treeFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.allelicInfoFileName = System.getProperty("variant.calling.file.prefix") + this.allelicInfoFileName;
            final File allelicInfoFile = new File(this.allelicInfoFileName);
            if (allelicInfoFile.exists())
                Log.err.println("Overwriting " + this.allelicInfoFileName + "...");
            this.allelicInfoOut = new PrintStream(this.allelicInfoFileName);

            if (System.getProperty("variant.calling.file.prefix") != null)
                this.sizeFactorFileName = System.getProperty("variant.calling.file.prefix") + this.sizeFactorFileName;
            final File sizeFactorFile = new File(this.sizeFactorFileName);
            if (sizeFactorFile.exists())
                Log.err.println("Overwriting " + this.sizeFactorFileName + "...");
            this.sizeFactorOut = new PrintStream(this.sizeFactorFileName);
        }
    } // openFiles

    public void init() {

        try {
            openFiles();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }

        logger.init(
                this.vcfOut,
                this.cellNamesOut,
                this.probsOut,
                this.treeOut,
                this.allelicInfoOut,
                this.sizeFactorOut
        );

    } // init

    public void log() {
        logger.log(this.vcfOut, this.adoOut, this.genotypesOut, this.ternaryOut, this.lociInfoOut, this.probsOut);

        if (this.treeOut != null) {
            logger.logTree(this.treeOut);
        }
    } // log

    public void close() {
        logger.close(this.vcfOut);

        this.vcfOut.close();

        if (this.cellNamesOut != null)
            this.cellNamesOut.close();

        if (this.ternaryOut != null)
            this.ternaryOut.close();

        if (this.lociInfoOut != null)
            this.lociInfoOut.close();

        if (this.probsOut != null)
            this.probsOut.close();

        if (this.treeOut != null)
            this.treeOut.close();

        if (this.allelicInfoOut != null)
            this.allelicInfoOut.close();

        if (this.sizeFactorOut != null)
            this.sizeFactorOut.close();
    } // close

}
