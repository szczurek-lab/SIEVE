package beast.evolution.alignment;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Description("Loci information of single sequence read counts data in an alignment")
public class ScsLociInfo extends BEASTObject {


    //**********************************************
    //*                   Inputs                   *
    //**********************************************

    final public Input<Boolean> lociInput = new Input<>("loci", "a flag representing whether " +
            "loci information existing or not.");

    final public Input<String> dataInput = new Input<>("value", "the information of each locus in a " +
            "default form of 'chromosome number,position number,reference nucleotide,the alternative nucleotide;" +
            "with the greatest number of read counts;' (note that comma is used to separate details at a locus and " +
            "semicolon is used to separate different loci).");


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    /**
     * a flag to show whether loci information is provided or not
     */
    protected boolean loci = false;

    /**
     * if loci information is offered, store it in lociList
     */
    protected List<VariantSiteInfo> lociList = null;

    /**
     * number of loci
     */
    protected int lociNr = 0;


    //**********************************************
    //*                Constructors                *
    //**********************************************

    /**
     * constructor for testing purpose
     */
    public ScsLociInfo() {
    }

    /**
     * constructor for testing purpose
     *
     * @param info
     * @param loci
     */
    public ScsLociInfo(boolean loci, String info) {
        lociInput.setValue(loci, this);
        dataInput.setValue(info, this);
        initAndValidate();
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        if (this.lociInput.get() != null) {
            this.loci = lociInput.get();
            if (this.loci) initLoci();
        }
    } // initAndValidate

    /**
     * save the loci information into an ArrayList consisting of ImmutableList
     */
    private void initLoci() {
        lociList = new ArrayList<>();

        String data = dataInput.get().trim().replaceAll("\\s+", "");

        // separate different loci
        String[] liInfo = data.split(";");

        for (String str : liInfo) {
            String[] lsInfo = str.split(",");

            if (lsInfo.length < 4)
                throw new IllegalArgumentException("Error! Insufficient amount of locus information or wrong format " +
                        "to properly define a variant site: " + str);

            // save locus information
            lociList.add(
                    new VariantSiteInfo(
                            lsInfo[0],
                            Long.parseLong(lsInfo[1]),
                            lsInfo[2],
                            Arrays.copyOfRange(lsInfo, 3, lsInfo.length)
                    )
            );
        }

        lociNr = lociList.size();
    } // initLoci

    /**
     * getter method for outside calls
     *
     * @return loci
     */
    public boolean isLoci() {
        return loci;
    }

    /**
     * getter method for outside calls
     *
     * @return lociList
     */
    public List<VariantSiteInfo> getLociList() {
        return lociList;
    }

    /**
     * getter method for outside calls
     *
     * @return lociNr
     */
    public int getLociNr() {
        return lociNr;
    }

    @Override
    public String toString() {
        if (loci) {
            return dataInput.get();
        } else {
            return null;
        }
    } // toString

} // class ScsLociInfo
