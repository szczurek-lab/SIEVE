package beast.app.variantcallingmodeladaptor.seqcovmodelmodifier;

import beast.app.variantcaller.EstimatesTypeCollection;
import org.w3c.dom.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ExploredSharedAllelicSeqCovModelModifier implements SeqCovModelModifier {

    final static private Pattern ID_REF_PATTERN = Pattern.compile("@(.+)");

    Document doc;
    Element seqCovModelNode;
    final List<Node> paraNodes;

    Map<String, Node> renewedAttributes;

    double initialCovVal, initialVarVal;


    public ExploredSharedAllelicSeqCovModelModifier(
            Document doc,
            Element seqCovModelNode,
            List<Node> paraNodes
    ) {
        this.doc = doc;
        this.seqCovModelNode = seqCovModelNode;

        this.paraNodes = new ArrayList<>();

        String initCov = null, initVar = null;

        // Get children of `seqCovModelNode`.
        NodeList childrenNodes = seqCovModelNode.getChildNodes();

        // Find parameters related to `allelicSeqCov` and `allelicSeqCovRawVar`.
        if (seqCovModelNode.hasAttribute("allelicSeqCov"))
            initCov = seqCovModelNode.getAttribute("allelicSeqCov");
        else {
            for (int i = 0; i < childrenNodes.getLength(); i++) {
                Node child = childrenNodes.item(i);

                if (child.getNodeName().equals("allelicSeqCov")) {
                    initCov = child.getNodeValue();
                    break;
                } else if (child.getNodeName().equals("parameter")) {
                    boolean matched = false;
                    NamedNodeMap attrs = child.getAttributes();

                    for (int j = 0; j < attrs.getLength(); j++) {
                        Node attr = attrs.item(j);

                        if (attr.getNodeName().equals("name") && attr.getNodeValue().equals("allelicSeqCov")) {
                            initCov = child.getTextContent();

                            seqCovModelNode.removeChild(child);

                            matched = true;
                            break;
                        }
                    }

                    if (matched)
                        break;
                }
            }
        }

        if (seqCovModelNode.hasAttribute("allelicSeqCovRawVar"))
            initVar = seqCovModelNode.getAttribute("allelicSeqCovRawVar");
        else {
            for (int i = 0; i < childrenNodes.getLength(); i++) {
                Node child = childrenNodes.item(i);

                if (child.getNodeName().equals("allelicSeqCovRawVar")) {
                    initVar = child.getNodeValue();
                    break;
                } else if (child.getNodeName().equals("parameter")) {
                    boolean matched = false;
                    NamedNodeMap attrs = child.getAttributes();

                    for (int j = 0; j < attrs.getLength(); j++) {
                        Node attr = attrs.item(j);

                        if (attr.getNodeName().equals("name") && attr.getNodeValue().equals("allelicSeqCovRawVar")) {
                            initVar = child.getTextContent();

                            seqCovModelNode.removeChild(child);

                            matched = true;
                            break;
                        }
                    }

                    if (matched)
                        break;
                }
            }
        }

        try {
            initialCovVal = Double.parseDouble(initCov);

            if (initialCovVal < 0)
                throw new IllegalArgumentException("Error! Initial value of 'allelicSeqCov' should be non-negative.");
        } catch (NumberFormatException e) {
            initialCovVal = -1;
        }

        try {
            initialVarVal = Double.parseDouble(initVar);

            if (initialCovVal < 0)
                throw new IllegalArgumentException("Error! Initial value of 'allelicSeqCovRawVar' should be non-negative.");
        } catch (NumberFormatException e) {
            initialVarVal = -1;
        }

        Map<String, String> refs = new HashMap<>();
        if (initialCovVal == -1)
            refs.put("allelicSeqCovArray", getIDFromRef(initCov));
        if (initialVarVal == -1)
            refs.put("allelicSeqCovRawVarArray", getIDFromRef(initVar));

        if (refs.size() > 0) {
            renewedAttributes = new HashMap<>();

            for (String i : refs.keySet()) {

                // Search for matched parameter node.
                for (Node j : paraNodes) {
                    boolean matched = false;

                    NamedNodeMap attrs = j.getAttributes();
                    for (int k = 0; k < attrs.getLength(); k++) {
                        Node attr = attrs.item(k);
                        if (attr.getNodeName().equals("id") && attr.getNodeValue().equals(refs.get(i))) {
                            renewedAttributes.put(i, j);

                            // Remove matched parameter node from `doc`.
                            j.getParentNode().removeChild(j);

                            matched = true;
                            break;
                        }
                    }

                    if (matched)
                        break;
                }

                // Remove matched prior node.
                NodeList priorNodes = this.doc.getElementsByTagName("prior");

                for (int j = 0; j < priorNodes.getLength(); j++) {
                    boolean matched = false;
                    Node priorNode = priorNodes.item(j);

                    NamedNodeMap attrs = priorNode.getAttributes();
                    for (int k = 0; k < attrs.getLength(); k++) {
                        Node attr = attrs.item(k);
                        if (attr.getNodeName().equals("x") && getIDFromRef(attr.getNodeValue()).equals(refs.get(i))) {
                            // Remove matched prior node from `doc`.
                            priorNode.getParentNode().removeChild(priorNode);

                            matched = true;
                        }
                    }

                    if (matched)
                        break;
                }
            }
        }
    }

    private String getIDFromRef(String val) {
        final Matcher matcher = ID_REF_PATTERN.matcher(val);

        if (matcher.matches())
            return matcher.group(1);
        else
            return val;
    } // getIDFromRef

    @Override
    public void preModify() {
        seqCovModelNode.setAttribute("inVariantCallingMode", "true");

        seqCovModelNode.removeAttribute("allelicSeqCov");
        seqCovModelNode.removeAttribute("allelicSeqCovRawVar");

        seqCovModelNode.setAttribute(
                "allelicSeqCovArray",
                initialCovVal == -1 ? "" : String.valueOf(initialCovVal)
        );

        seqCovModelNode.setAttribute(
                "allelicSeqCovRawVarArray",
                initialVarVal == -1 ? "" : String.valueOf(initialVarVal)
        );
    }

    @Override
    public void postModify(
            final String[] lociInfo,
            final Map<String, Map<String, double[][]>> allelicSeqEstimates,
            final Map<String, String> sitesMap,
            final EstimatesTypeCollection adaptedConfigEstimateType
    ) {
        if (renewedAttributes == null || renewedAttributes.size() == 0)
            return;

        for (String i : renewedAttributes.keySet()) {
            Node para = renewedAttributes.get(i);

            if (((Element) para).hasAttribute("value"))
                seqCovModelNode.setAttribute(i, ((Element) para).getAttribute("value"));
            else {
                // Get the only child.
                seqCovModelNode.setAttribute(i, para.getChildNodes().item(0).getNodeValue());
            }
        }
    }

}
