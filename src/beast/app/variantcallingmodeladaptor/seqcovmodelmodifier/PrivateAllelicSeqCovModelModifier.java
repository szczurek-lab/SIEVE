package beast.app.variantcallingmodeladaptor.seqcovmodelmodifier;

import beast.app.variantcaller.EstimatesTypeCollection;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.Text;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PrivateAllelicSeqCovModelModifier implements SeqCovModelModifier {

    Document doc;
    Element seqCovModelNode;

    Node covNode;
    Node varNode;


    public PrivateAllelicSeqCovModelModifier(
            Document doc,
            Element seqCovModelNode
    ) {
        this.doc = doc;
        this.seqCovModelNode = seqCovModelNode;
    }

    @Override
    public void preModify() {
        seqCovModelNode.setAttribute("inVariantCallingMode", "true");

        seqCovModelNode.removeAttribute("allelicSeqCov");
        seqCovModelNode.removeAttribute("allelicSeqCovRawVar");

        covNode = this.doc.createElement("allelicSeqCovArray");
        varNode = this.doc.createElement("allelicSeqCovRawVarArray");

        Text covText = this.doc.createTextNode("");
        Text varText = this.doc.createTextNode("");

        covNode.appendChild(covText);
        varNode.appendChild(varText);

        seqCovModelNode.appendChild(covNode);
        seqCovModelNode.appendChild(varNode);
    }

    @Override
    public void postModify(
            final String[] lociInfo,
            final Map<String, Map<String, double[][]>> allelicSeqEstimates,
            final Map<String, String> sitesMap,
            final EstimatesTypeCollection adaptedConfigEstimateType
    ) {
        List<String> cov = new ArrayList<>();
        List<String> var = new ArrayList<>();
        int j = 0, numOfMatrices = 0;
        do {
            StringBuilder covStr = new StringBuilder();
            StringBuilder varStr = new StringBuilder();

            int k = 0;
            for (String site : lociInfo) {
                Map<String, double[][]> estimates = allelicSeqEstimates.get(sitesMap.get(site));
                double[][] data;
                if (adaptedConfigEstimateType.getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEDIAN ||
                        adaptedConfigEstimateType.getEstimatesType() == EstimatesTypeCollection.EstimatesType.MEAN)
                    data = estimates.get(adaptedConfigEstimateType.getEstimatesType().toString().toLowerCase());
                else
                    data = estimates.get(adaptedConfigEstimateType.getModeKDEType().toString().toLowerCase());

                if (numOfMatrices == 0)
                    numOfMatrices = data.length;
                else {
                    if (numOfMatrices != data.length)
                        throw new IllegalArgumentException("The number of matrices of all sites in allelic sequencing information are not consistent. Check the provided data.");
                }

                covStr.append(data[j][0]);
                varStr.append(data[j][1]);

                if (k < lociInfo.length - 1) {
                    covStr.append(",");
                    varStr.append(",");
                }

                k++;
            }

            cov.add(covStr.toString());
            var.add(varStr.toString());

            j++;
        } while (j != numOfMatrices);

        String covStr = String.join(";", cov);
        String varStr = String.join(";", var);

        // Add a 'text' type child node to 'allelicSeqCovArray' and 'allelicSeqCovRawVarArray', respectively
        assert covNode != null;
        assert varNode != null;
        covNode.getFirstChild().setTextContent(covStr);
        varNode.getFirstChild().setTextContent(varStr);
    }

}
