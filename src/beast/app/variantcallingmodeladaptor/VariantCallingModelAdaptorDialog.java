package beast.app.variantcallingmodeladaptor;

import beast.app.util.WholeNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.io.File;

public class VariantCallingModelAdaptorDialog {

    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private File configFile = null;
    private File inputTreeFile = null;
    private File estimatesFile = null;
    private File outputFile = null;

    private final JRadioButton useMeanRateRadioButton = new JRadioButton("Use mean rate", true);
    private final JRadioButton useMedianRateRadioButton = new JRadioButton("Use median rate", false);

    private final JRadioButton useBranchLengthRadioButton = new JRadioButton("Only use branch length", true);
    private final JRadioButton useBranchLengthAndRateRadioButton = new JRadioButton("Use both branch length and rate", false);

    private JCheckBox saveDetailsCheckBox = new JCheckBox();

    private final WholeNumberField cellThresholdText = new WholeNumberField(1, Integer.MAX_VALUE);


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public VariantCallingModelAdaptorDialog(final JFrame frame) {
        this.frame = frame;
        this.optionPanel = new OptionsPanel(10, 10);

        // Choose configuration template file
        final JButton configFileButton = new JButton("Choose file...");
        final JTextField configFileNameText = new JTextField("not selected", 30);

        configFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select configuration template file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            configFile = new File(dialog.getDirectory(), dialog.getFile());
            configFileNameText.setText(configFile.getName());
        });
        configFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(configFileNameText, BorderLayout.CENTER);
        panel1.add(configFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Configuration template file: ", panel1);

        // Choose input best tree file
        final JButton inputTreeFileButton = new JButton("Choose file...");
        final JTextField inputTreeFileNameText = new JTextField("not selected", 30);

        inputTreeFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select input tree file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            inputTreeFile = new File(dialog.getDirectory(), dialog.getFile());
            inputTreeFileNameText.setText(inputTreeFile.getName());
        });
        inputTreeFileNameText.setEditable(false);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.add(inputTreeFileNameText, BorderLayout.CENTER);
        panel2.add(inputTreeFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Input tree file: ", panel2);

        // Choose cached estimates file
        final JButton estimatesFileButton = new JButton("Choose file...");
        final JTextField estimatesFileNameText = new JTextField("not selected", 30);

        estimatesFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select cached estimates file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            estimatesFile = new File(dialog.getDirectory(), dialog.getFile());
            estimatesFileNameText.setText(estimatesFile.getName());
        });
        estimatesFileNameText.setEditable(false);

        JPanel panel3 = new JPanel(new BorderLayout(0, 0));
        panel3.add(estimatesFileNameText, BorderLayout.CENTER);
        panel3.add(estimatesFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Cached estimates file: ", panel3);

        this.optionPanel.addSeparator();

        JLabel title = new JLabel("Optional");
        title.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title);

        // Choose output file (can be empty)
        final JButton outputFileButton = new JButton("Choose file...");
        final JTextField outputFileNameText = new JTextField("not selected", 30);

        outputFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select output file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            outputFile = new File(dialog.getDirectory(), dialog.getFile());
            outputFileNameText.setText(outputFile.getName());
        });
        outputFileNameText.setEditable(false);

        JPanel panel4 = new JPanel(new BorderLayout(0, 0));
        panel4.add(outputFileNameText, BorderLayout.CENTER);
        panel4.add(outputFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Output file: ", panel4);

        cellThresholdText.setColumns(12);
        cellThresholdText.setValue(1);
        cellThresholdText.setVisible(true);
        this.optionPanel.addComponentWithLabel("Mutated cell number threshold: ", cellThresholdText);

        this.optionPanel.addComponentWithLabel("Save details: ", this.saveDetailsCheckBox);

        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(useMeanRateRadioButton);
        buttonGroup1.add(useMedianRateRadioButton);
        JPanel panel5 = new JPanel(new BorderLayout(0, 0));
        panel5.add(useMeanRateRadioButton, BorderLayout.CENTER);
        panel5.add(useMedianRateRadioButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Branch rate estimates: ", panel5);

        ButtonGroup buttonGroup2 = new ButtonGroup();
        buttonGroup2.add(useBranchLengthRadioButton);
        buttonGroup2.add(useBranchLengthAndRateRadioButton);
        JPanel panel6 = new JPanel(new BorderLayout(0, 0));
        panel6.add(useBranchLengthRadioButton, BorderLayout.CENTER);
        panel6.add(useBranchLengthAndRateRadioButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Evolutionary time from input tree: ", panel6);
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
    } // showDialog

    public String getConfigFileName() {
        if (configFile == null) return null;
        return configFile.getPath();
    }

    public String getInputTreeFileName() {
        if (inputTreeFile == null) return null;
        return inputTreeFile.getPath();
    }

    public String getEstimatesFileName() {
        if (estimatesFile == null) return null;
        return estimatesFile.getPath();
    }

    public String getOutputFileName() {
        if (outputFile == null) return null;
        return outputFile.getPath();
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
