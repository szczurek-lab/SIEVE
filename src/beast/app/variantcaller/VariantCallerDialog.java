package beast.app.variantcaller;

import beast.app.util.WholeNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.io.File;
import java.util.Objects;

public class VariantCallerDialog {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private File cachedEstimatesFile = null;
    private File modelConfigFile = null;
    private File mcmcSamplesFile = null;
    private File allelicInfoFile = null;
    private File gtAdoSamplesFile = null;
    private File inputTreeFile = null;
    private File outputVCFile = null;
    // private File variantCallingLogFile = null;

    private final JFrame frame;

    private final OptionsPanel optionPanel;

    private final WholeNumberField threadsText = new WholeNumberField(1, 1000);

    private final WholeNumberField burninText = new WholeNumberField(0, Integer.MAX_VALUE);

    private final JCheckBox cachedEstimatesCheckBox = new JCheckBox();

    private final JComboBox<EstimatesTypeCollection.EstimatesType> estimateTypeCombo = new JComboBox<>(EstimatesTypeCollection.EstimatesType.values());

    private final JComboBox<EstimatesTypeCollection.ModeKDEType> modeKDETypeJCombo = new JComboBox<>(EstimatesTypeCollection.ModeKDEType.values());

    private final JCheckBox saveDetailsCheckBox = new JCheckBox();

    private final WholeNumberField cellThresholdText = new WholeNumberField(1, Integer.MAX_VALUE);

    private final JRadioButton useMeanRateRadioButton = new JRadioButton("Use mean rate", true);
    private final JRadioButton useMedianRateRadioButton = new JRadioButton("Use median rate", false);

    private final JRadioButton useBranchLengthRadioButton = new JRadioButton("Only use branch length", true);
    private final JRadioButton useBranchLengthAndRateRadioButton = new JRadioButton("Use both branch length and rate", false);

    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    /**
     * constructor
     *
     * @param combined    whether this is called from PostProcessor (combined with ScsTreeAnnotator; True) or alone (False)
     * @param frame       apparently
     * @param optionPanel apparently
     */
    public VariantCallerDialog(boolean combined, final JFrame frame, final OptionsPanel optionPanel) {
        this.frame = frame;
        this.optionPanel = new OptionsPanel(12, 12);

        if (combined) {
            JLabel title = new JLabel("Variant Calling Section");
            title.setFont(new Font("Default", Font.ITALIC, 14));
            this.optionPanel.addSpanningComponent(title);
        }

        // Input threads value
        threadsText.setColumns(12);
        threadsText.setValue(1);
        threadsText.setVisible(true);
        this.optionPanel.addComponentWithLabel("Threads: ", threadsText);

        // Input burnin value
        burninText.setColumns(12);
        burninText.setValue(0);
        if (combined) {
            burninText.setVisible(false);
        } else {
            burninText.setVisible(true);
            this.optionPanel.addComponentWithLabel("Burnin percentage: ", burninText);
        }

        this.optionPanel.addSeparator();

        // Load cached estimates or not
        this.optionPanel.addComponentWithLabel("Load cached estimates: ", cachedEstimatesCheckBox);

        // Choose cached estimates file
        final JButton cachedEstimatesFileButton = new JButton("Choose file...");
        final JTextField cachedEstimatesFileNameText = new JTextField("not selected", 30);

        cachedEstimatesFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select cached estimates file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            cachedEstimatesFile = new File(dialog.getDirectory(), dialog.getFile());
            cachedEstimatesFileNameText.setText(cachedEstimatesFile.getName());
        });
        cachedEstimatesFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(cachedEstimatesFileNameText, BorderLayout.CENTER);
        panel1.add(cachedEstimatesFileButton, BorderLayout.EAST);
        final JLabel label1 = this.optionPanel.addComponentWithLabel("Cached estimates file: ", panel1);
        label1.setEnabled(false);
        cachedEstimatesFileButton.setEnabled(false);
        cachedEstimatesFileNameText.setEnabled(false);

        // Choose estimates type
        final JLabel label2 = this.optionPanel.addComponentWithLabel("Estimates from MCMC samples: ", estimateTypeCombo);
        estimateTypeCombo.setEnabled(true);

        // Choose KDE distribution for mode estimates
        final JLabel label3 = this.optionPanel.addComponentWithLabel("KDE distribution for mode estimates: ", modeKDETypeJCombo);
        label3.setEnabled(false);
        modeKDETypeJCombo.setEnabled(false);

        estimateTypeCombo.addItemListener(itemEvent -> {
            boolean selected = Objects.requireNonNull(estimateTypeCombo.getSelectedItem()).toString().equals(EstimatesTypeCollection.EstimatesType.MODE.toString());
            label3.setEnabled(selected);
            modeKDETypeJCombo.setEnabled(selected);
        });

        // Choose MCMC samples log file
        final JButton mcmcSamplesFileButton = new JButton("Choose file...");
        final JTextField mcmcSamplesFileNameText = new JTextField("not selected", 30);

        mcmcSamplesFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select MCMC samples log file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            mcmcSamplesFile = new File(dialog.getDirectory(), dialog.getFile());
            mcmcSamplesFileNameText.setText(mcmcSamplesFile.getName());
        });
        mcmcSamplesFileNameText.setEditable(false);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.add(mcmcSamplesFileNameText, BorderLayout.CENTER);
        panel2.add(mcmcSamplesFileButton, BorderLayout.EAST);
        final JLabel label4 = this.optionPanel.addComponentWithLabel("MCMC samples log file: ", panel2);

        // Choose allelic sequencing log file
//        final JButton allelicInfoFileButton = new JButton("Choose file...");
//        final JTextField allelicInfoFileNameText = new JTextField("not selected", 30);
//
//        allelicInfoFileButton.addActionListener(ae -> {
//            FileDialog dialog = new FileDialog(frame,
//                    "Select allelic sequencing log file...",
//                    FileDialog.LOAD);
//
//            dialog.setVisible(true);
//            if (dialog.getFile() == null) {
//                return;
//            }
//
//            allelicInfoFile = new File(dialog.getDirectory(), dialog.getFile());
//            allelicInfoFileNameText.setText(allelicInfoFile.getName());
//        });
//        allelicInfoFileNameText.setEditable(false);
//
//        JPanel panel3 = new JPanel(new BorderLayout(0, 0));
//        panel3.add(allelicInfoFileNameText, BorderLayout.CENTER);
//        panel3.add(allelicInfoFileButton, BorderLayout.EAST);
//        final JLabel label5 = this.optionPanel.addComponentWithLabel("Allelic sequencing log file: ", panel3);

        // Choose genotypes and ado states log file
//        final JButton gtAdoSamplesFileButton = new JButton("Choose file...");
//        final JTextField gtAdoSamplesFileNameText = new JTextField("not selected", 30);
//
//        gtAdoSamplesFileButton.addActionListener(ae -> {
//            FileDialog dialog = new FileDialog(frame,
//                    "Select genotypes and ado states log file...",
//                    FileDialog.LOAD);
//
//            dialog.setVisible(true);
//            if (dialog.getFile() == null) {
//                return;
//            }
//
//            gtAdoSamplesFile = new File(dialog.getDirectory(), dialog.getFile());
//            gtAdoSamplesFileNameText.setText(gtAdoSamplesFile.getName());
//        });
//        gtAdoSamplesFileNameText.setEditable(false);
//
//        JPanel panel4 = new JPanel(new BorderLayout(0, 0));
//        panel4.add(gtAdoSamplesFileNameText, BorderLayout.CENTER);
//        panel4.add(gtAdoSamplesFileButton, BorderLayout.EAST);
//        final JLabel label6 = this.optionPanel.addComponentWithLabel("Genotypes and ado states log file: ", panel4);

        cachedEstimatesCheckBox.addItemListener(ItemEvent -> {
            boolean selected = cachedEstimatesCheckBox.isSelected();

            label1.setEnabled(selected);
            cachedEstimatesFileButton.setEnabled(selected);
            cachedEstimatesFileNameText.setEnabled(selected);

            cachedEstimatesFile = null;
            cachedEstimatesFileNameText.setText("not selected");

            label2.setEnabled(!selected);
            estimateTypeCombo.setEnabled(!selected);

            label3.setEnabled(false);
            modeKDETypeJCombo.setEnabled(false);

            label4.setEnabled(!selected);
            mcmcSamplesFileButton.setEnabled(!selected);
            mcmcSamplesFileNameText.setEnabled(!selected);

//            label5.setEnabled(!selected);
//            allelicInfoFileButton.setEnabled(!selected);
//            allelicInfoFileNameText.setEnabled(!selected);

//            label6.setEnabled(!selected);
//            gtAdoSamplesFileButton.setEnabled(!selected);
//            gtAdoSamplesFileNameText.setEnabled(!selected);

            estimateTypeCombo.setSelectedIndex(0);

            modeKDETypeJCombo.setSelectedIndex(0);

            mcmcSamplesFile = null;
            mcmcSamplesFileNameText.setText("not selected");

//            allelicInfoFile = null;
//            allelicInfoFileNameText.setText("not selected");

//            gtAdoSamplesFile = null;
//            gtAdoSamplesFileNameText.setText("not selected");
        });

        this.optionPanel.addSeparator();

        // Choose model configuration file
        final JButton modelConfigFileButton = new JButton("Choose file...");
        final JTextField modelConfigFileNameText = new JTextField("not selected", 30);

        modelConfigFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select model configuration file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            modelConfigFile = new File(dialog.getDirectory(), dialog.getFile());
            modelConfigFileNameText.setText(modelConfigFile.getName());
        });
        modelConfigFileNameText.setEditable(false);

        JPanel panel5 = new JPanel(new BorderLayout(0, 0));
        panel5.add(modelConfigFileNameText, BorderLayout.CENTER);
        panel5.add(modelConfigFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Model configuration file: ", panel5);

        // Choose best inferred tree
        String inputTreeFileName;
        if (combined) {
            inputTreeFileName = "Use the best tree summarized from ScsTreeAnnotator...";
        } else {
            inputTreeFileName = "not selected";
        }
        final JButton inputTreeFileButton = new JButton("Choose file...");
        final JTextField inputTreeFileNameText = new JTextField(inputTreeFileName, 30);

        inputTreeFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select best inferred tree file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            inputTreeFile = new File(dialog.getDirectory(), dialog.getFile());
            inputTreeFileNameText.setText(inputTreeFile.getName());
        });
        inputTreeFileNameText.setEditable(false);

        JPanel panel6 = new JPanel(new BorderLayout(0, 0));
        panel6.add(inputTreeFileNameText, BorderLayout.CENTER);
        panel6.add(inputTreeFileButton, BorderLayout.EAST);
        final JLabel label7 = this.optionPanel.addComponentWithLabel("Best inferred tree file: ", panel6);

        if (combined) {
            label7.setEnabled(false);
            inputTreeFileNameText.setEditable(false);
            inputTreeFileButton.setEnabled(false);
        }

        this.optionPanel.addSeparator();

        JLabel title = new JLabel("Optional");
        title.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title);

        // Choose output variant calling file
        final JButton outputVCFileButton = new JButton("Choose file...");
        final JTextField outputVCFileNameText = new JTextField("not selected", 30);

        outputVCFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select output variant calling file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            outputVCFile = new File(dialog.getDirectory(), dialog.getFile());
            outputVCFileNameText.setText(outputVCFile.getName());
        });
        outputVCFileNameText.setEditable(false);

        JPanel panel7 = new JPanel(new BorderLayout(0, 0));
        panel7.add(outputVCFileNameText, BorderLayout.CENTER);
        panel7.add(outputVCFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Output variant calling file: ", panel7);

        cellThresholdText.setColumns(12);
        cellThresholdText.setValue(1);
        cellThresholdText.setVisible(true);
        this.optionPanel.addComponentWithLabel("Mutated cell number threshold: ", cellThresholdText);

        this.optionPanel.addComponentWithLabel("Save details: ", this.saveDetailsCheckBox);

        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(useMeanRateRadioButton);
        buttonGroup1.add(useMedianRateRadioButton);
        JPanel panel8 = new JPanel(new BorderLayout(0, 0));
        panel8.add(useMeanRateRadioButton, BorderLayout.CENTER);
        panel8.add(useMedianRateRadioButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Branch rate estimates: ", panel8);

        ButtonGroup buttonGroup2 = new ButtonGroup();
        buttonGroup2.add(useBranchLengthRadioButton);
        buttonGroup2.add(useBranchLengthAndRateRadioButton);
        JPanel panel9 = new JPanel(new BorderLayout(0, 0));
        panel9.add(useBranchLengthRadioButton, BorderLayout.CENTER);
        panel9.add(useBranchLengthAndRateRadioButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Evolutionary time from input tree: ", panel9);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public boolean showDialog(String title) {

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                new String[]{"Run", "Quit"},
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, title);
        //dialog.setResizable(true);
        dialog.pack();

        dialog.setVisible(true);

        return optionPane.getValue().equals("Run");
    }

    public int getNumOfThreads() {
        return threadsText.getValue();
    }

    public int getBurninText() {
        return burninText.getValue();
    }

    public boolean isLoadCachedEstimates() {
        return cachedEstimatesCheckBox.isSelected();
    }

    public String getCachedEstimatesFileName() {
        return cachedEstimatesFile.getPath();
    }

    public String getModelConfigFileName() {
        if (modelConfigFile == null) return null;
        return modelConfigFile.getPath();
    }

    public String getMCMCSamplesFileName() {
        if (mcmcSamplesFile == null) return null;
        return mcmcSamplesFile.getPath();
    }

    public String getAllelicInfoFileName() {
        if (allelicInfoFile == null) return null;
        return allelicInfoFile.getPath();
    }

    public String getGtAdoSamplesFileName() {
        if (gtAdoSamplesFile == null) return null;
        return gtAdoSamplesFile.getPath();
    }

    public String getInputTreeFileName() {
        if (inputTreeFile == null) return null;
        return inputTreeFile.getPath();
    }

    public String getOutputVCFileName() {
        if (outputVCFile == null) return null;
        return outputVCFile.getPath();
    }

    public EstimatesTypeCollection.EstimatesType getEstimatesType() {
        return (EstimatesTypeCollection.EstimatesType) estimateTypeCombo.getSelectedItem();
    }

    public EstimatesTypeCollection.ModeKDEType getModeKDEType() {
        return (EstimatesTypeCollection.ModeKDEType) modeKDETypeJCombo.getSelectedItem();
    }

    public boolean getSaveDetails() {
        return this.saveDetailsCheckBox.isSelected();
    }

    public int getCellThreshold() {
        return this.cellThresholdText.getValue();
    }

    public boolean isUseMedianRate() {
        return useMedianRateRadioButton.isSelected();
    }

    public boolean isUseBranchLengthAndRate() {
        return useBranchLengthAndRateRadioButton.isSelected();
    }

}
