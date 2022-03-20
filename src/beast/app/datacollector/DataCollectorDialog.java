package beast.app.datacollector;

import beast.app.util.WholeNumberField;
import beast.math.util.MathFunctions;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.io.File;
import java.util.Arrays;

public class DataCollectorDialog {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public enum DataType {
        CovSup,
        FullSupsCov
    }

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private final JCheckBox compatibleWithSciPhiCheckBox = new JCheckBox("Compatible with SCIPhI", false);
    private final JCheckBox useSexCheckBox = new JCheckBox("Use candidate sites from sex chromosomes", true);

    private final JRadioButton CSRadioButton = new JRadioButton("Coverage-Support", false);
    private final JRadioButton FSCRadioButton = new JRadioButton("Full supports-Coverage", true);

    private final JCheckBox sampleCheckBox = new JCheckBox("", true);
    private final WholeNumberField sampledCellsText = new WholeNumberField(-1, 100000);
    private final WholeNumberField sampledLociText = new WholeNumberField(-1, 100000);

    private final JCheckBox backgroundCheckBox = new JCheckBox("", true);
    private final JTextField backgroundText = new JTextField();

    private File cellNamesFile = null;
    private File excludedCellNamesFile = null;
    private File dataFile = null;
    private File templateFile = null;
    private File outputFile = null;


    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public DataCollectorDialog(JFrame frame) {
        this.frame = frame;
        this.optionPanel = new OptionsPanel(10, 10);

        // Choose cell names file
        final JButton cellNamesFileButton = new JButton("Choose file...");
        final JTextField cellNamesFileNameText = new JTextField("not selected", 35);

        cellNamesFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select cell names file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            cellNamesFile = new File(dialog.getDirectory(), dialog.getFile());
            cellNamesFileNameText.setText(cellNamesFile.getName());
        });
        cellNamesFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(cellNamesFileNameText, BorderLayout.WEST);
        panel1.add(cellNamesFileButton, BorderLayout.CENTER);
        panel1.add(compatibleWithSciPhiCheckBox, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Cell names file: ", panel1);

        // Choose data file
        final JButton dataFileButton = new JButton("Choose file...");
        final JTextField dataFileNameText = new JTextField("not selected", 35);

        dataFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select data file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            dataFile = new File(dialog.getDirectory(), dialog.getFile());
            dataFileNameText.setText(dataFile.getName());
        });
        dataFileNameText.setEditable(false);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.add(dataFileNameText, BorderLayout.WEST);
        panel2.add(dataFileButton, BorderLayout.CENTER);
        panel2.add(useSexCheckBox, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Data file: ", panel2);

        // Set data type
        CSRadioButton.setActionCommand("Coverage-Support");
        FSCRadioButton.setActionCommand("Full supports-Coverage");

        ButtonGroup buttonGroup = new ButtonGroup();
        buttonGroup.add(CSRadioButton);
        buttonGroup.add(FSCRadioButton);

        JPanel panel3 = new JPanel(new BorderLayout(0, 0));
        panel3.add(CSRadioButton, BorderLayout.WEST);
        panel3.add(FSCRadioButton, BorderLayout.CENTER);
        this.optionPanel.addComponentWithLabel("Input datatype", panel3);

        // Choose template file
        final JButton templateFileButton = new JButton("Choose file...");
        final JTextField templateFileNameText = new JTextField("not selected", 35);

        templateFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select template configuration file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            templateFile = new File(dialog.getDirectory(), dialog.getFile());
            templateFileNameText.setText(templateFile.getName());
        });
        templateFileNameText.setEditable(false);

        JPanel panel4 = new JPanel(new BorderLayout(0, 0));
        panel4.add(templateFileNameText, BorderLayout.CENTER);
        panel4.add(templateFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Template configuration file: ", panel4);

        /*JLabel comment = new JLabel("Note that in the template configuration file the section containing " +
                "alignment should have a tag named 'data', and the log sections should have at least a tag named 'logger'.");
        comment.setFont(new Font("Default", Font.ITALIC, 10));
        this.optionPanel.addSpanningComponent(comment);*/

        this.optionPanel.addSeparator();

        JLabel title1 = new JLabel("Optional: excluded cell names file");
        title1.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title1);

        // Choose healthy cell names file
        final JButton excludedCellNamesFileButton = new JButton("Choose file...");
        final JTextField excludedCellNamesFileNameText = new JTextField("not selected", 35);

        excludedCellNamesFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select healthy cell names file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            excludedCellNamesFile = new File(dialog.getDirectory(), dialog.getFile());
            excludedCellNamesFileNameText.setText(excludedCellNamesFile.getName());
        });
        excludedCellNamesFileNameText.setEditable(false);

        JPanel panel5 = new JPanel(new BorderLayout(0, 0));
        panel5.add(excludedCellNamesFileNameText, BorderLayout.CENTER);
        panel5.add(excludedCellNamesFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Excluded cell names file: ", panel5);

        this.optionPanel.addSeparator();

        JLabel title2 = new JLabel("Optional: output file");
        title2.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title2);

        // Choose output updated configuration file
        final JButton outputFileButton = new JButton("Choose file...");
        final JTextField outputFileNameText = new JTextField("not selected", 35);

        outputFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select output configuration file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            outputFile = new File(dialog.getDirectory(), dialog.getFile());
            outputFileNameText.setText(outputFile.getName());
        });
        outputFileNameText.setEditable(false);

        JPanel panel6 = new JPanel(new BorderLayout(0, 0));
        panel6.add(outputFileNameText, BorderLayout.CENTER);
        panel6.add(outputFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Output configuration file: ", panel6);

        this.optionPanel.addSeparator();

        JLabel title3 = new JLabel("Optional: sample from data; default: no sampling");
        title3.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title3);

        // Sample from data or not
        this.optionPanel.addComponentWithLabel("Use default: ", sampleCheckBox);

        // Add sample control
        sampledCellsText.setColumns(12);
        sampledCellsText.setValue(-1);
        sampledCellsText.setVisible(true);
        sampledCellsText.setEnabled(false);
        JLabel label1 = this.optionPanel.addComponentWithLabel("Number of cells to sample: ", sampledCellsText);
        label1.setEnabled(false);

        sampledLociText.setColumns(12);
        sampledLociText.setValue(-1);
        sampledLociText.setVisible(true);
        sampledLociText.setEnabled(false);
        JLabel label2 = this.optionPanel.addComponentWithLabel("Number of loci to sample: ", sampledLociText);
        label2.setEnabled(false);

        sampleCheckBox.addItemListener(ItemEvent -> {
            boolean selected = sampleCheckBox.isSelected();

            sampledCellsText.setValue(-1);
            sampledCellsText.setEnabled(!selected);
            label1.setEnabled(!selected);

            sampledLociText.setValue(-1);
            sampledLociText.setEnabled(!selected);
            label2.setEnabled(!selected);
        });

        this.optionPanel.addSeparator();

        JLabel title5 = new JLabel("<html>Optional: the order of background information appearing in data file " +
                "default: <br/>" +
                "For \"Full supports-Coverage\" datatype: 0 - variant1, 1 - variant2, 2 - variant3, " +
                "3 - normal, 4 - coverage" +
                "<br/>" +
                "For \"Coverage-Support\" datatype: 0 - coverage, 1 - variant, 2 - normal</html>");
        title5.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title5);

        // The order of background information
        this.optionPanel.addComponentWithLabel("Use default: ", backgroundCheckBox);

        // Add background coverage control
        backgroundText.setColumns(12);
        if (CSRadioButton.isSelected())
            backgroundText.setText("0 1 2");
        if (FSCRadioButton.isSelected())
            backgroundText.setText("0 1 2 3 4");
        backgroundText.setVisible(true);
        backgroundText.setEnabled(false);
        JLabel label4 = this.optionPanel.addComponentWithLabel("Order (whitespaces separated): ", backgroundText);
        label4.setEnabled(false);

        backgroundCheckBox.addItemListener(ItemEvent -> {
            boolean activated = backgroundCheckBox.isSelected();
            boolean csSelected = CSRadioButton.isSelected();
            boolean fscSelected = FSCRadioButton.isSelected();

            if (csSelected)
                backgroundText.setText("0 1 2");

            if (fscSelected)
                backgroundText.setText("0 1 2 3 4");

            backgroundText.setEnabled(!activated);
            label4.setEnabled(!activated);
        });

        CSRadioButton.addItemListener(ItemEvent -> {
            if (!label4.isEnabled())
                backgroundText.setText("0 1 2");
        });

        FSCRadioButton.addItemListener(ItemEvent -> {
            if (!label4.isEnabled())
                backgroundText.setText("0 1 2 3 4");
        });

    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    public boolean showDialog(String title) {

        JOptionPane optionPane = new JOptionPane(
                optionPanel,
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

    public String getCellNamesFileName() {
        if (cellNamesFile == null) return null;
        return cellNamesFile.getPath();
    }

    public boolean isCellNamesCompatibleWithSCIPhI() {
        return compatibleWithSciPhiCheckBox.isSelected();
    }

    public boolean isUseSex() {
        return useSexCheckBox.isSelected();
    }

    public String getDataFileName() {
        if (dataFile == null) return null;
        return dataFile.getPath();
    }

    public String getTemplateFileName() {
        if (templateFile == null) return null;
        return templateFile.getPath();
    }

    public String getExcludedCellNamesFile() {
        if (excludedCellNamesFile == null) return null;
        return excludedCellNamesFile.getPath();
    }

    public String getOutputFileName() {
        if (outputFile == null) return null;
        return outputFile.getPath();
    }

    public int[] getSamplesControl() {
        int[] results = new int[2];

        results[0] = sampledCellsText.getValue();
        results[1] = sampledLociText.getValue();

        return results;
    }

    public DataType getDataType() {
        if (CSRadioButton.isSelected())
            return DataType.CovSup;

        if (FSCRadioButton.isSelected())
            return DataType.FullSupsCov;

        return null;
    }

    public int[] getBackgroundInfoOrder() throws IllegalArgumentException {
        String[] tmp = backgroundText.getText().trim().split("\\s+");

        if (CSRadioButton.isSelected() && tmp.length != 3) {
            throw new IllegalArgumentException("The order of background information must be a combination of 0, 1, " +
                    "and 2, separated by whitespaces.");
        }

        if (FSCRadioButton.isSelected() && tmp.length != 5) {
            throw new IllegalArgumentException("The order of background information must be a combination of 0, 1, " +
                    "2, 3, and 4, separated by whitespaces.");
        }

        int[] results = new int[tmp.length];
        Arrays.fill(results, -1);

        for (int i = 0; i < tmp.length; i++) {
            int val = Integer.parseInt(tmp[i]);

            if (results[i] == -1)
                results[i] = val;
            else
                throw new IllegalArgumentException("Error: Duplicate background information index provided: " + val);
        }

        final int max = MathFunctions.max(
                Arrays.stream(results).boxed().toArray(Integer[]::new)
        );

        final int min = MathFunctions.min(
                Arrays.stream(results).boxed().toArray(Integer[]::new)
        );

        if (CSRadioButton.isSelected()) {
            if (max > 2 || min < 0)
                throw new IllegalArgumentException("Error: Only expecting combinations of 0, 1, and 2, but provided " +
                        Arrays.toString(tmp));

            return Arrays.copyOfRange(results, 0, 3);
        }

        if (FSCRadioButton.isSelected()) {
            if (max > 4 || min < 0)
                throw new IllegalArgumentException("Error: Only expecting combinations of 0, 1, 2, 3, and 4, but " +
                        "provided " + Arrays.toString(tmp));

            return results;
        }

        return null;
    }

}
