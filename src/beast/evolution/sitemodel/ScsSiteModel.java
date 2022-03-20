package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.ReadCounts;
import beast.evolution.likelihood.ScsTreeLikelihood;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;

@Description(value = "An extended site model including an error model applied to leaves of a tree")
public class ScsSiteModel extends SiteModel {

    ReadCounts m_scsDataType;

    @Override
    public void initAndValidate() {
        useBeast1StyleGamma = true; //useBeast1StyleGammaInput.get();
        muParameter = muParameterInput.get();
        if (muParameter == null) {
            muParameter = new RealParameter("1.0");
        }
        shapeParameter = shapeParameterInput.get();
        invarParameter = invarParameterInput.get();
        if (invarParameter == null) {
            invarParameter = new RealParameter("0.0");
            invarParameter.setBounds(Math.max(0.0, invarParameter.getLower()), Math.min(1.0, invarParameter.getUpper()));
        }

        //if (muParameter != null) {
        muParameter.setBounds(Math.max(muParameter.getLower(), 0.0), Math.min(muParameter.getUpper(), Double.POSITIVE_INFINITY));
        //}
        if (shapeParameter != null) {
            // The quantile calculator fails when the shape parameter goes much below
            // 1E-3 so we have put a hard lower bound on it. If this is not there then
            // the category rates can go to 0 and cause a -Inf likelihood (whilst this
            // is not a problem as the state will be rejected, it could mask other issues
            // and this seems the better approach.
            shapeParameter.setBounds(Math.max(shapeParameter.getLower(), 1.0E-3), Math.min(shapeParameter.getUpper(), 1.0E3));
        }


        if (/*invarParameter != null && */(invarParameter.getValue() < 0 || invarParameter.getValue() > 1)) {
            throw new IllegalArgumentException("proportion invariant should be between 0 and 1");
        }
        refresh();

        addCondition(muParameterInput);
        addCondition(invarParameterInput);
        addCondition(shapeParameterInput);
    }

    public void setDataType(final ReadCounts dataType) {
        m_scsDataType = dataType;
    }

    public boolean canSetSubstModel(final ScsSubstitutionModelBase substModel) {
        if (m_scsDataType == null) {
            // try to find out the data type from the data in a treelikelihood in an output
            for (Object beastObject : getOutputs()) {
                if (beastObject instanceof ScsTreeLikelihood) {
                    ScsTreeLikelihood likelihood = (ScsTreeLikelihood) beastObject;
                    m_scsDataType = likelihood.scsDataInput.get().getDataType();
                    break;
                }
            }
        } else {
            return substModel.canHandleDataType(m_scsDataType);
        }
        return true;
    }
}
