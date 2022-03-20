package beast.evolution.variantsinfo;

import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import org.jetbrains.annotations.NotNull;

import java.io.PrintStream;

public abstract class GenericVariantsInfoLog extends GenericVariantsInfo.Base {


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initialise(int numOfThreads, ScsSubstitutionModelBase substModel) {
        throw new IllegalArgumentException("Unsupported function.");
    } // initialise

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
    @Override
    public void cloze(@NotNull PrintStream out) {
        throw new IllegalArgumentException("Unsupported function.");
    } // close

}
