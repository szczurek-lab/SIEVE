package beast.evolution.likelihood;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

public class GenotypeAdoStateLogger extends BEASTObject implements Loggable {

    private final static String SITE_SEPARATOR = "\t";
    private final static String CELL_SEPARATOR = ";";


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    final public Input<ScsGenericTreeLikelihood> treeLikelihoodInput = new Input<>(
            "treeLikelihood",
            "intermediate variables to be logged",
            Input.Validate.REQUIRED
    );


    //*********************************************
    //*            Implemented methods            *
    //*********************************************

    @Override
    public void initAndValidate() {
        if (!treeLikelihoodInput.get().isTraceMLGenotypes())
            throw new IllegalArgumentException("Error! Maximum likelihood genotypes and ADO states are not " +
                    "being traced. Set attribute `traceMLGenotypes` of `treeLikelihood` to true.");
    } // initAndValidate

    @Override
    public void init(PrintStream out) {
        treeLikelihoodInput.get().initSiteHeader(1, out);
        out.println();
    } // init

    @Override
    public void log(long sample, PrintStream out) {
        treeLikelihoodInput.get().logGenotypeAndAdo(out, SITE_SEPARATOR, CELL_SEPARATOR);
        out.println();
    } // log

    @Override
    public void close(PrintStream out) {
        out.println("#Sampled genotypes and ADO states for all sites and cells.");
        out.println("#Format: for each line (sample), sites are tab-separated, and cells are semi colon-separated with such information: genotype,ADO state");
        out.println("#Cells: " + treeLikelihoodInput.get().getSortedCellNamesFromVarInfo(","));
        out.println("#Sites map (chromosome number,locus number,reference nucleotide,alternative nucleotides)");
        treeLikelihoodInput.get().closeGenotypeAndAdo(1, out);
    } // close

}
