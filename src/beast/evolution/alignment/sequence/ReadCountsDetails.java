package beast.evolution.alignment.sequence;

public class ReadCountsDetails {

    final char[] nucs;
    final int[] reads;
    final int refReads;
    final int cov;

    public ReadCountsDetails(char[] nucs, int[] reads, int refReads, int cov) {
        if (nucs.length != reads.length)
            throw new RuntimeException("Length should be the same.");

        this.nucs = new char[nucs.length];
        this.reads = new int[reads.length];
        System.arraycopy(nucs, 0, this.nucs, 0, nucs.length);
        System.arraycopy(reads, 0, this.reads, 0, reads.length);

        this.refReads = refReads;
        this.cov = cov;
    }

    public char getAltNuc(int pos) {
        return this.nucs[pos];
    }

    public char[] getAltNucs() {
        return this.nucs;
    }

    public int getAltReads(int pos) {
        return this.reads[pos];
    }

    public int getRefReads() {
        return this.refReads;
    }

    public int getCoverage() {
        return this.cov;
    }

}
