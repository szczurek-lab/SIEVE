package beast.evolution.alignment;

import beast.app.utils.ChromosomeLabel;
import org.apache.commons.lang.ArrayUtils;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VariantSiteInfo implements Comparable<VariantSiteInfo> {

    private final ChromosomeLabel chrom;
    private final long pos;
    private final char refNuc;
    private final char[] altNucs;
    private final char[] orderedAltNucs;

    public VariantSiteInfo(
            String chr,
            long pos,
            String refNuc,
            String[] altNucs
    ) {
        if (refNuc.trim().length() != 1 || !Character.isAlphabetic(refNuc.trim().charAt(0)))
            throw new IllegalArgumentException("Error! Illegal reference nucleotide detected: " + refNuc);

        List<Character> tmp = new ArrayList<>();
        for (String i : altNucs) {
            if (i.trim().length() != 1 || !Character.isAlphabetic(i.trim().charAt(0)))
                throw new IllegalArgumentException("Error! Illegal alternative nucleotide detected: " + i);

            if (isACGT(i.trim().charAt(0)))
                tmp.add(i.trim().charAt(0));
        }

        this.chrom = new ChromosomeLabel(chr);
        this.pos = pos;
        this.refNuc = refNuc.trim().charAt(0);
        this.altNucs = ArrayUtils.toPrimitive(tmp.toArray(new Character[0]));

        this.orderedAltNucs = new char[this.altNucs.length];
        System.arraycopy(this.altNucs, 0, this.orderedAltNucs, 0, this.altNucs.length);

        for (int i = 0; i < this.altNucs.length; i++) {
            this.orderedAltNucs[i] = Character.toUpperCase(this.orderedAltNucs[i]);
        }
        Arrays.sort(this.orderedAltNucs);
    }

    public VariantSiteInfo(
            String chr,
            long pos,
            char refNuc,
            char[] altNucs
    ) {
        this.chrom = new ChromosomeLabel(chr);
        this.pos = pos;
        this.refNuc = refNuc;
        this.altNucs = new char[altNucs.length];
        System.arraycopy(altNucs, 0, this.altNucs, 0, altNucs.length);

        this.orderedAltNucs = new char[altNucs.length];
        System.arraycopy(altNucs, 0, this.orderedAltNucs, 0, altNucs.length);

        for (int i = 0; i < altNucs.length; i++) {
            this.orderedAltNucs[i] = Character.toUpperCase(this.orderedAltNucs[i]);
        }
        Arrays.sort(this.orderedAltNucs);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VariantSiteInfo that = (VariantSiteInfo) o;

        if (!this.chrom.equals(that.chrom)) return false;

        if (this.pos != that.pos) return false;

        if (this.refNuc != that.refNuc) return false;

        if (this.altNucs.length != that.altNucs.length) return false;

        for (int i = 0; i < this.altNucs.length; i++) {
            if (this.orderedAltNucs[i] != that.orderedAltNucs[i]) return false;
        }

        return true;
    } // equal

    @Override
    public int hashCode() {
        int result = chrom.hashCode();
        result = 31 * result + (int) (pos ^ (pos >>> 32));
        result = 31 * result + (int) refNuc;
        result = 31 * result + Arrays.hashCode(orderedAltNucs);
        return result;
    }

    public String toString(String separator, boolean usePlaceHolder) {
        StringBuilder sb = new StringBuilder();

        sb.append(this.chrom.getFullLabel()).append(separator);
        sb.append(this.pos).append(separator);
        sb.append(this.refNuc);

        if (this.altNucs == null || this.altNucs.length == 0) {
            if (usePlaceHolder)
                sb.append(separator).append('N');
        } else {
            for (char i : this.altNucs)
                sb.append(separator).append(i);
        }

        return sb.toString();
    } // toString

    public String getChromosome() {
        return this.chrom.getFullLabel();
    }

    public ChromosomeLabel getChromosomeLabel() {
        return this.chrom;
    }

    public long getPosition() {
        return this.pos;
    }

    public char getRefNuc() {
        return this.refNuc;
    }

    public char[] getAltNucs() {
        return this.altNucs;
    }

    public int getAltNucsLength() {
        return this.altNucs.length;
    }

    public String getAltNucsAsString(String separator) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < this.altNucs.length; i++) {
            sb.append(this.altNucs[i]);

            if (i < this.altNucs.length - 1)
                sb.append(separator);
        }

        return sb.toString();
    }

    public static boolean isACGT(char val) {
        final char upper = Character.toUpperCase(val);
        return upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T';
    }

    @Override
    public int compareTo(@NotNull VariantSiteInfo v) {
        final int chrC = this.chrom.compareTo(v.chrom);
        if (chrC != 0) return chrC;

        if (this.pos > v.pos)
            return 1;
        else if (this.pos < v.pos)
            return -1;

        return 0;
    }

}
