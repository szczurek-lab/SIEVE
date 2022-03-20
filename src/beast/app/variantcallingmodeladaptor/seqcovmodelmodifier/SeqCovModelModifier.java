package beast.app.variantcallingmodeladaptor.seqcovmodelmodifier;

import beast.app.variantcaller.EstimatesTypeCollection;

import java.util.Map;

public interface SeqCovModelModifier {

    public void preModify();

    public void postModify(
            final String[] lociInfo,
            final Map<String, Map<String, double[][]>> allelicSeqEstimates,
            final Map<String, String> sitesMap,
            final EstimatesTypeCollection adaptedConfigEstimateType
    );

}
