package beast.evolution.variantsinfo;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.evolution.alignment.ScsAlignment;
import beast.evolution.substitutionmodel.ScsSubstitutionModelBase;
import beast.evolution.tree.TreeInterface;

import java.io.PrintStream;
import java.util.List;

@Description("Collection of variants information in variant calling")
public interface GenericVariantsInfo {

    abstract class Base extends CalculationNode implements GenericVariantsInfo {

        public enum AmbiguityGenotypesStrategy {
            Default,
            KeepAll;

            public static AmbiguityGenotypesStrategy[] toArray() {
                return new AmbiguityGenotypesStrategy[]{Default, KeepAll};
            }
        }


        //***********************************************
        //*                  Variables                  *
        //***********************************************

        protected ScsAlignment scsData;

        protected TreeInterface tree;

        protected int numOfTips;
        protected int numOfNodes;

        protected String[] sortedTipNames;
        protected int[] originalTipNamesIndices;

        protected ScsSubstitutionModelBase substModel;

        protected AmbiguityGenotypesStrategy ambiguityGenotypesStrategy;


        //**********************************************
        //*              Abstract methods              *
        //**********************************************

        public abstract void initialise(int value, ScsSubstitutionModelBase substModel);

        public abstract void deeplyInitialise();

        public abstract void addMLCategory(int patternIndex, List<Integer> values, boolean overwrite);

        public abstract void addMLCategory(int patternIndex, int value, boolean overwrite);

        public abstract void addVariantsInfo(int index, GenericVariantsInfo.Base value);

        public abstract List<Integer> getMLCategory(int patternIndex);

        public abstract void initialiseAPattern(int patternIndex);

        /**
         * Add maximum likelihood genotypes of tips.
         *
         * @param patternIndex which pattern?
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addMLGenotypesTips(int patternIndex, final int[] values, boolean overwrite);

        public abstract void addMLGenotypesTips(int patternIndex, final List<int[]> values, boolean overwrite);

        /**
         * Add maximum likelihood genotypes of all nodes.
         *
         * @param patternIndex which pattern?
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addMLGenotypesNodes(int patternIndex, final int[] values, boolean overwrite);

        public abstract void addMLGenotypesNodes(int patternIndex, final List<int[]> values, boolean overwrite);

        /**
         * Add maximum likelihood ADO state of tips.
         *
         * @param patternIndex which pattern?
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addMLAdo(int patternIndex, final int[] values, boolean overwrite);

        public abstract void addMLAdo(int patternIndex, final List<int[]> values, boolean overwrite);

        /**
         * Add log likelihood of a pattern being constant (matrix wise).
         *
         * @param patternIndex which pattern?
         * @param value        apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addGenotypeLogLikelihoodConstantPattern(int patternIndex, double value, boolean overwrite);

        public abstract void addGenotypeLogLikelihoodConstantPattern(int patternIndex, List<Double> values, boolean overwrite);

        /**
         * Add log likelihoods of tips.
         *
         * @param patternIndex which pattern?
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addGenotypeLogLikelihoodsTips(int patternIndex, final double[] values, boolean overwrite);

        public abstract void addGenotypeLogLikelihoodsTips(int patternIndex, final List<double[]> values, boolean overwrite);

        /**
         * Add log likelihoods of all nodes.
         *
         * @param patternIndex which pattern?
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addGenotypeLogLikelihoodsNodes(int patternIndex, final double[] values, boolean overwrite);

        public abstract void addGenotypeLogLikelihoodsNodes(int patternIndex, final List<double[]> values, boolean overwrite);

        /**
         * Add log likelihoods of all genotypes for a specific tip at a specific pattern.
         *
         * @param patternIndex which pattern?
         * @param tipIndex     which tip?
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addGenotypeLogLikelihoodsAll(int patternIndex, int tipIndex, final double[] values, boolean overwrite);

        public abstract void addGenotypeLogLikelihoodsAll(int patternIndex, int tipIndex, final List<double[]> values, boolean overwrite);

        /**
         * Add log transfer probabilities of a parent node from its maximum likelihood genotype to all genotypes of a child tip.
         *
         * @param tipIndex       which tip?
         * @param parentGenotype maximum likelihood genotype of parent node
         * @param numOfGenotypes the number of genotypes
         * @param values         apparently
         * @param overwrite      overwrite current values or not
         */
        public abstract void addGenotypeLogTransferProbabilitiesAll(int tipIndex, int parentGenotype, int numOfGenotypes, final double[] values, boolean overwrite);

        /**
         * Add log likelihoods of ADOs.
         *
         * @param patternIndex which pattern?
         * @param tipIndex     which tip?
         * @param numOfAdos    the number of allowed ADOs
         * @param values       apparently
         * @param overwrite    overwrite current values or not
         */
        public abstract void addAdoLogLikelihoods(int patternIndex, int tipIndex, int numOfAdos, double[] values, boolean overwrite);

        /**
         * Not `Loggable`.
         *
         * @param out apparently
         */
        public abstract void init(PrintStream out);

        /**
         * Not `Loggable`.
         *
         * @param out apparently
         */
        public abstract void log(PrintStream out, String siteSeparator, String cellSeparator);

        /**
         * Not `Loggable`.
         *
         * @param out apparently
         */
        public abstract void cloze(PrintStream out);


        //***********************************************
        //*                   Methods                   *
        //***********************************************

        @Override
        public void initAndValidate() {
        }

        public String getSortedTipNames(String separator) {
            StringBuilder ret = new StringBuilder();

            for (int i = 0; i < sortedTipNames.length; i++) {
                ret.append(sortedTipNames[i]);

                if (i < sortedTipNames.length - 1)
                    ret.append(separator);
            }

            return ret.toString();
        } // getSortedTipNames

        /**
         * adjust the number of ambiguous genotypes to log in the VCF file:
         * {@link AmbiguityGenotypesStrategy} Default: keep only the first item
         * {@link AmbiguityGenotypesStrategy} KeepAll: keep all items
         *
         * @param value number of ambiguous genotypes to log
         * @return adjusted number of ambiguous genotypes to log
         */
        protected int adjustAmbiguousGenotypes(int value) {
            assert value > 0;

            if (this.ambiguityGenotypesStrategy == AmbiguityGenotypesStrategy.Default)
                return 1;
            else if (this.ambiguityGenotypesStrategy == AmbiguityGenotypesStrategy.KeepAll)
                return value;
            else
                throw new IllegalArgumentException("Error! Unknown strategy to processing ambiguous genotypes.");
        } // adjustAmbiguousGenotypes

        /**
         * @param tipName the taxon name as a string
         * @param data    the alignment
         * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
         * or -1 if the taxon is not in the alignment.
         */
        protected int getTaxonIndex(String tipName, ScsAlignment data) {
            int taxonIndex = data.getTaxonIndex(tipName);
            if (taxonIndex == -1) {
                if (tipName.startsWith("'") || tipName.startsWith("\""))
                    taxonIndex = data.getTaxonIndex(tipName.substring(1, tipName.length() - 1));

                if (taxonIndex == -1)
                    throw new RuntimeException("Could not find sequence " + tipName + " in the ScsAlignment (" +
                            this.getClass().getName() + ")");
            }
            return taxonIndex;
        } // getTaxonIndex

    } // abstract class Base

}
