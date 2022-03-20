package beast.app.utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChromosomeLabel implements Comparable<ChromosomeLabel> {

    private final static String FORMAT_ERROR = "The format of chromosome fullLabel should be pure chromosome names, " +
            "such as 1, 20, and X (case insensitive), with or without a prefix, which should not contain any " +
            "numbers. The detected fullLabel violates this rule: ";

    private final static String CHROMOSOME_FORMAT = "^([a-zA-Z]*)([1-9]|1\\d|2[012]|[XYxy]|MT|mt)$";
    private final static String CHROMOSOME_FORMAT_WITH_PREFIX = "^([a-zA-Z]+)([1-9]|1\\d|2[012]|[XYxy]|MT|mt)$";
    private final static String CHROMOSOME_FORMAT_WITHOUT_PREFIX = "^([1-9]|1\\d|2[012]|[XYxy]|MT|mt)$";
    private final static String CHROMOSOME_FORMAT_MT = "^([a-zA-Z]*)(MT|mt)$";
    private final static String CHROMOSOME_FORMAT_SEX = "^([a-zA-Z]*)([XYxy])$";

    private final static Pattern PATTERN = Pattern.compile(CHROMOSOME_FORMAT);
    private final static Pattern PATTERN_WITH_PREFIX = Pattern.compile(CHROMOSOME_FORMAT_WITH_PREFIX);
    private final static Pattern PATTERN_WITHOUT_PREFIX = Pattern.compile(CHROMOSOME_FORMAT_WITHOUT_PREFIX);
    private final static Pattern PATTERN_MT = Pattern.compile(CHROMOSOME_FORMAT_MT);
    private final static Pattern PATTERN_SEX = Pattern.compile(CHROMOSOME_FORMAT_SEX);

    private final String fullLabel;
    private final String label;
    private final boolean withPrefix;
    private final Matcher matcher;

    public ChromosomeLabel(String fullLabel) throws IllegalArgumentException {
        if (!matches(fullLabel))
            throw new IllegalArgumentException(FORMAT_ERROR + fullLabel);

        this.fullLabel = fullLabel;
        this.matcher = PATTERN.matcher(this.fullLabel);
        this.label = getChromosomeLabel();
        this.withPrefix = matchesWithPrefix(fullLabel);
    }

    public String getFullLabel() {
        return this.fullLabel;
    } // getFullLabel

    public String getLabel() {
        return this.label;
    } // getLabel

    public boolean isOnMitochondrial() {
        return PATTERN_MT.matcher(fullLabel).matches();
    } // isOnMitochondrial

    public boolean isOnSex() {
        return PATTERN_SEX.matcher(fullLabel).matches();
    } // isOnSex

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ChromosomeLabel that = (ChromosomeLabel) o;

        return this.label.equals(that.label);
    } // equals

    @Override
    public int hashCode() {
        return label.hashCode();
    } // hashCode

    public boolean isWithPrefix() {
        return this.withPrefix;
    } // isWithPrefix

    @Override
    public int compareTo(ChromosomeLabel cl1) {
        final int i;
        final int i1;
        try {
            i = Integer.parseInt(this.label);
            i1 = Integer.parseInt(cl1.label);
        } catch (NumberFormatException e) {
            return Integer.compare(this.label.compareTo(cl1.label), 0);
        }

        return Integer.compare(i, i1);
    } // compareTo

    private String getChromosomeLabel() {
        if (this.matcher.find(0))
            return this.matcher.group(2);
        else
            return null;
    } // getChromosomeLabel

    public static boolean matches(String fullLabel) {
        return PATTERN.matcher(fullLabel).matches();
    } // matches

    public static boolean matchesWithPrefix(String fullLabel) {
        return PATTERN_WITH_PREFIX.matcher(fullLabel).matches();
    } // matchesWithPrefix

    public static boolean matchesWithoutPrefix(String fullLabel) {
        return PATTERN_WITHOUT_PREFIX.matcher(fullLabel).matches();
    } // matchesWithoutPrefix

}
