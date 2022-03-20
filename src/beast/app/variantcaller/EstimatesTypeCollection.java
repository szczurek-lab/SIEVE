package beast.app.variantcaller;

public class EstimatesTypeCollection {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public enum EstimatesType {
        MEDIAN("Median"),
        MEAN("Mean"),
        MODE("Mode (Maximum a posteriori)"),
        ALL("All (Consider all situations)");

        String desc;

        EstimatesType(String s) {
            desc = s;
        }

        @Override
        public String toString() {
            return desc;
        }
    }

    public enum ModeKDEType {
        GAUSSIAN("Gaussian")/*,
        EPANECHNIKOV("Epanechnikov"),
        RECTANGULAR("Rectangular"),
        TRIANGULAR("Triangular"),
        BIWEIGHT("Biweight"),
        COSINE("Cosine"),
        OPTCOSINE("Optcosine")*/;

        String desc;

        ModeKDEType(String s) {
            desc = s;
        }

        @Override
        public String toString() {
            return desc;
        }
    }

    private EstimatesType estimatesType;
    private ModeKDEType modeKDEType;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public EstimatesTypeCollection() {
        this.estimatesType = null;
        this.estimatesType = null;
    }

    public EstimatesTypeCollection(EstimatesType estimatesType, ModeKDEType modeKDEType) {
        this.estimatesType = estimatesType;
        this.modeKDEType = modeKDEType;
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public void setEstimateType(EstimatesType estimateType) {
        this.estimatesType = estimateType;
    }

    public void setModeKDEType(ModeKDEType modeKDEType) {
        this.modeKDEType = modeKDEType;
    }

    public EstimatesType getEstimatesType() {
        return this.estimatesType;
    }

    public ModeKDEType getModeKDEType() {
        return this.modeKDEType;
    }

}
