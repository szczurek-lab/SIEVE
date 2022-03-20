package beast.core.variantslogger;

import java.io.PrintStream;

public interface VariantsLoggable {

    /**
     * write header information
     *
     * @param vcfOut         vcf log stream
     * @param cellNamesOut   cell names log stream
     * @param probsOut       posterior probability log stream
     * @param treeOut        annotated tree log stream
     * @param allelicInfoOut allelic information log stream
     */
    void init(
            PrintStream vcfOut,
            PrintStream cellNamesOut,
            PrintStream probsOut,
            PrintStream treeOut,
            PrintStream allelicInfoOut
    );

    /**
     * log each entry to PrintStream
     *
     * @param vcfOut       vcf log stream
     * @param adoOut       ado state log stream
     * @param genotypesOut inferred genotypes at tips log stream
     * @param ternaryOut   ternary log stream
     * @param lociInfoOut  loci information log stream
     * @param probsOut     posterior probability log stream
     */
    void log(
            PrintStream vcfOut,
            PrintStream adoOut,
            PrintStream genotypesOut,
            PrintStream ternaryOut,
            PrintStream lociInfoOut,
            PrintStream probsOut
    );

    /**
     * log tree with meta data
     *
     * @param treeOut annotated tree log stream
     */
    void logTree(PrintStream treeOut);

    /**
     * close log.
     *
     * @param vcfOut vcf log stream
     */
    void close(PrintStream vcfOut);

}
