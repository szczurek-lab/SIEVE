package beast.app.geneannotator;

import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GeneAnnotatorDialog {

    private final static Map<String, Integer> SUBST_MODELS;
    private final static String DEFAULT_SUBST_MODEL;

    static {
        SUBST_MODELS = new HashMap<>();
        SUBST_MODELS.put("ScsFiniteMuExtendedModel", 0);
        SUBST_MODELS.put("ScsFiniteMuModel", 1);
        SUBST_MODELS.put("ScsFiniteMuDelModel", 2);
        SUBST_MODELS.put("ScsFiniteMuDelInsModel", 3);

        DEFAULT_SUBST_MODEL = "ScsFiniteMuExtendedModel";
    }


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private final JComboBox<String> substModelCombo = new JComboBox<>(SUBST_MODELS.keySet().toArray(new String[0]));

    private final JCheckBox fromAnnovarCheckBox = new JCheckBox("From Annovar", true);
    private final JRadioButton tabSepRadioButton1 = new JRadioButton("Tab-separated", true);
    private final JRadioButton commaSepRadioButton1 = new JRadioButton("Comma-separated", false);

    private final JRadioButton tabSepRadioButton2 = new JRadioButton("Tab-separated", true);
    private final JRadioButton commaSepRadioButton2 = new JRadioButton("Comma-separated", false);

    private final JTextField filteringGenesColText = new JTextField();

    private File treeFile = null;
    private File snvFile = null;
    private File mutationMapFile = null;
    private File filteringGenesFile = null;
    private File outputFile = null;

    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public GeneAnnotatorDialog(final JFrame frame) {
        this.frame = frame;
        this.optionPanel = new OptionsPanel(10, 10);

        // Choose substitution model
        substModelCombo.setSelectedItem(DEFAULT_SUBST_MODEL);
        this.optionPanel.addComponentWithLabel("Estimates from MCMC samples: ", substModelCombo);
        substModelCombo.setEnabled(true);

        // Choose tree file
        final JButton treeFileButton = new JButton("Choose file...");
        final JTextField treeFileNameText = new JTextField("not selected", 30);

        treeFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select tree file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            treeFile = new File(dialog.getDirectory(), dialog.getFile());
            treeFileNameText.setText(treeFile.getName());
        });
        treeFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(treeFileNameText, BorderLayout.CENTER);
        panel1.add(treeFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Tree file: ", panel1);

        // Choose variant sites (SNV) file
        final JButton snvFileButton = new JButton("Choose file...");
        final JTextField snvFileNameText = new JTextField("not selected", 30);

        snvFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select SNV file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            snvFile = new File(dialog.getDirectory(), dialog.getFile());
            snvFileNameText.setText(snvFile.getName());
        });
        snvFileNameText.setEditable(false);

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.add(snvFileNameText, BorderLayout.CENTER);
        panel2.add(snvFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("SNV file: ", panel2);

        // Choose mutation map file
        final JButton mutationMapFileButton = new JButton("Choose file...");
        final JTextField mutationMapFileNameText = new JTextField("not selected", 30);

        mutationMapFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select mutation map file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            mutationMapFile = new File(dialog.getDirectory(), dialog.getFile());
            mutationMapFileNameText.setText(mutationMapFile.getName());
        });
        mutationMapFileNameText.setEditable(false);

        JPanel panel3 = new JPanel(new BorderLayout(0, 0));
        panel3.add(mutationMapFileNameText, BorderLayout.CENTER);
        panel3.add(mutationMapFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Mutation map file: ", panel3);

        // Whether the mutation map file from Annovar or not?
        tabSepRadioButton1.setActionCommand("Tab");
        commaSepRadioButton1.setActionCommand("Comma");

        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(tabSepRadioButton1);
        buttonGroup1.add(commaSepRadioButton1);

        fromAnnovarCheckBox.addActionListener(ae -> {
            if (fromAnnovarCheckBox.isSelected()) {
                tabSepRadioButton1.setEnabled(true);
                commaSepRadioButton1.setEnabled(true);
            } else {
                tabSepRadioButton1.setEnabled(false);
                commaSepRadioButton1.setEnabled(false);
            }
        });

        JPanel panel4 = new JPanel(new BorderLayout(0, 0));
        panel4.add(fromAnnovarCheckBox, BorderLayout.WEST);
        panel4.add(tabSepRadioButton1, BorderLayout.CENTER);
        panel4.add(commaSepRadioButton1, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Mutation map file option", panel4);

        this.optionPanel.addSeparator();

        JLabel title = new JLabel("Optional");
        title.setFont(new Font("Default", Font.ITALIC, 14));
        this.optionPanel.addSpanningComponent(title);

        // Choose filtering genes file (can be empty)
        final JButton filteringGenesButton = new JButton("Choose file...");
        final JTextField filteringGenesFileNameText = new JTextField("not selected", 30);

        filteringGenesButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select filtering genes file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                return;
            }

            filteringGenesFile = new File(dialog.getDirectory(), dialog.getFile());
            filteringGenesFileNameText.setText(filteringGenesFile.getName());
        });
        filteringGenesFileNameText.setEditable(false);

        JPanel panel5 = new JPanel(new BorderLayout(0, 0));
        panel5.add(filteringGenesFileNameText, BorderLayout.CENTER);
        panel5.add(filteringGenesButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Filtering genes file (headers required): ", panel5);

        // Filtering genes file options (can be empty)
        ButtonGroup buttonGroup2 = new ButtonGroup();
        buttonGroup2.add(tabSepRadioButton2);
        buttonGroup2.add(commaSepRadioButton2);
        JPanel panel6 = new JPanel(new BorderLayout(0, 0));
        panel6.add(tabSepRadioButton2, BorderLayout.CENTER);
        panel6.add(commaSepRadioButton2, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Filtering genes file option: ", panel6);

        // Column index of filtering genes in the corresponding file
        filteringGenesColText.setColumns(12);
        filteringGenesColText.setText("0");
        filteringGenesColText.setVisible(true);
        filteringGenesColText.setEnabled(true);
        this.optionPanel.addComponentWithLabel("Column index of filtering genes (starting from 0): ",
                filteringGenesColText);

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

        JPanel panel7 = new JPanel(new BorderLayout(0, 0));
        panel7.add(outputFileNameText, BorderLayout.CENTER);
        panel7.add(outputFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Output file: ", panel7);
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

    public String getTreeFileName() {
        if (treeFile == null) return null;
        return treeFile.getPath();
    }

    public String getSNVFileName() {
        if (snvFile == null) return null;
        return snvFile.getPath();
    }

    public int getSubstModelLabel() {
        return SUBST_MODELS.get(substModelCombo.getSelectedItem());
    }

    public String getMutationMapFileName() {
        if (mutationMapFile == null) return null;
        return mutationMapFile.getPath();
    }

    public boolean isMapFromAnnovar() {
        return fromAnnovarCheckBox.isSelected();
    }

    public int getMapSeparator() {
        if (isMapFromAnnovar()) {
            if (tabSepRadioButton1.isSelected())
                return 0;
            else if (commaSepRadioButton1.isSelected())
                return 1;
        }

        return -1;
    }

    public String getFilteringFileName() {
        if (filteringGenesFile == null) return null;
        return filteringGenesFile.getPath();
    }

    public int getFilteringFileSeparator() {
        if (tabSepRadioButton2.isSelected())
            return 0;
        else if (commaSepRadioButton2.isSelected())
            return 1;

        return -1;
    }

    public int getColIndexFilteringGenes() {
        String col = filteringGenesColText.getText().trim();
        final Matcher matcher = Pattern.compile("\\d+").matcher(col);

        if (!matcher.matches())
            throw new IllegalArgumentException("Error, column index of filtering genes should be a number.");

        return Integer.parseInt(col);
    }

    public String getOutputFileName() {
        if (outputFile == null) return null;
        return outputFile.getPath();
    }

}
