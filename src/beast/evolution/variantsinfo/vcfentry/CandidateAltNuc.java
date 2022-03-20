package beast.evolution.variantsinfo.vcfentry;

public class CandidateAltNuc {

    private final char nuc;
    private int count;
    private int numOfCells;

    public CandidateAltNuc(char nuc) {
        this.nuc = nuc;
        this.count = 0;
        this.numOfCells = 0;
    }

    public void addCount(int value) {
        assert value >= 0;
        this.count += value;
    } // addCount

    public void addNumOfCells(int value) {
        assert value >= 0;
        this.numOfCells += value;
    } // addNumOfCells

    public char getNuc() {
        return this.nuc;
    } // getNuc

    public int getCount() {
        return this.count;
    } // getCount

    public int getNumOfCells() {
        return this.numOfCells;
    } // getNumOfCells

}
