package beast.evolution.rawreadcountsmodel;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.likelihood.ScsGenericTreeLikelihood;

import java.io.PrintStream;

public class CoverageRawVarianceLogger extends BEASTObject implements Loggable {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    final public Input<ScsGenericTreeLikelihood> treeLikelihoodInput = new Input<>("treeLikelihood",
            "intermediate variables to be logged", Input.Validate.REQUIRED);


    private boolean updateSeqCovModel;

    //*********************************************
    //*            Implemented methods            *
    //*********************************************

    @Override
    public void initAndValidate() {
        updateSeqCovModel = treeLikelihoodInput.get().updateSeqCovModel();
    } // initAndValidate

    @Override
    public void init(PrintStream out) {
        if (!updateSeqCovModel) return;

        treeLikelihoodInput.get().initSiteHeader(1, out);
        out.println();
    } // init

    @Override
    public void log(long sample, PrintStream out) {
        if (!updateSeqCovModel) return;

        treeLikelihoodInput.get().logCovar(out);
    } // log

    @Override
    public void close(PrintStream out) {
        if (!updateSeqCovModel) return;

        out.println("#Sampled allelic sequencing coverage and raw variance.");
        out.println("#Format (for each semicolon separated site-wise rate variation category, if applicable, " +
                "and tab separated site): allelic sequencing coverage,raw variance");
        out.println("#Map (chromosome number,locus number,reference nucleotide,alternative nucleotide)");
        treeLikelihoodInput.get().closeCovar(1, out);
    } // close

}
