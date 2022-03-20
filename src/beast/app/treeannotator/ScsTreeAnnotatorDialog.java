package beast.app.treeannotator;

import beast.app.util.WholeNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;

public class ScsTreeAnnotatorDialog {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private final WholeNumberField burninText = new WholeNumberField(0, Integer.MAX_VALUE);
    private final RealNumberField limitText = new RealNumberField(0.0, 1.0);

    private final JComboBox<ScsTreeAnnotator.Target> summaryTreeCombo = new JComboBox<>(ScsTreeAnnotator.Target.values());

    private final JComboBox<ScsTreeAnnotator.HeightsSummary> nodeHeightsCombo = new JComboBox<>(ScsTreeAnnotator.HeightsSummary.values());

    private final JCheckBox lowMemCheckbox = new JCheckBox();
    private final JCheckBox simpleTreeCheckBox = new JCheckBox();
    private final JCheckBox hasTrunkCheckBox = new JCheckBox("", true);

    private File targetFile = null;
    private File inputFile = null;
    private File outputFile = null;


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
    public ScsTreeAnnotatorDialog(boolean combined, final JFrame frame, final OptionsPanel optionPanel) {
        this.frame = frame;
        this.optionPanel = new OptionsPanel(12, 12);

        if (combined) {
            JLabel title = new JLabel("Tree Annotation Section");
            title.setFont(new Font("Default", Font.ITALIC, 14));
            this.optionPanel.addSpanningComponent(title);
        }

        burninText.setColumns(12);
        burninText.setValue(0);
        if (combined) {
            burninText.setVisible(false);
        } else {
            burninText.setVisible(true);
            this.optionPanel.addComponentWithLabel("Burnin percentage: ", burninText);
        }

        limitText.setColumns(12);
        limitText.setValue(0.0);
        this.optionPanel.addComponentWithLabel("Posterior probability limit: ", limitText);

        this.optionPanel.addComponentWithLabel("Target tree type: ", summaryTreeCombo);
        this.optionPanel.addComponentWithLabel("Node heights: ", nodeHeightsCombo);

        this.optionPanel.addSeparator();

        final JButton targetFileButton = new JButton("Choose File...");
        final JTextField targetFileNameText = new JTextField("not selected", 16);

        targetFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select target file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                // the dialog was cancelled...
                return;
            }

            targetFile = new File(dialog.getDirectory(), dialog.getFile());
            targetFileNameText.setText(targetFile.getName());

        });
        targetFileNameText.setEditable(false);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.add(targetFileNameText, BorderLayout.CENTER);
        panel1.add(targetFileButton, BorderLayout.EAST);
        final JLabel label1 = this.optionPanel.addComponentWithLabel("Target tree file: ", panel1);

        JButton inputFileButton = new JButton("Choose File...");
        final JTextField inputFileNameText = new JTextField("not selected", 16);

        inputFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select input tree file...",
                    FileDialog.LOAD);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                // the dialog was cancelled...
                return;
            }

            inputFile = new File(dialog.getDirectory(), dialog.getFile());
            inputFileNameText.setText(inputFile.getName());

        });
        inputFileNameText.setEditable(false);

        label1.setEnabled(false);
        targetFileNameText.setEnabled(false);
        targetFileButton.setEnabled(false);

        summaryTreeCombo.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent itemEvent) {
                boolean selected = summaryTreeCombo.getSelectedItem().toString().equals("User target tree");
                label1.setEnabled(selected);
                targetFileNameText.setEnabled(selected);
                targetFileButton.setEnabled(selected);
            }
        });

        JPanel panel2 = new JPanel(new BorderLayout(0, 0));
        panel2.add(inputFileNameText, BorderLayout.CENTER);
        panel2.add(inputFileButton, BorderLayout.EAST);

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        new FileDrop(null, inputFileNameText, focusBorder, new FileDrop.Listener() {
            @Override
            public void filesDropped(java.io.File[] files) {
                inputFile = files[0];
                inputFileNameText.setText(inputFile.getName());
            }   // end filesDropped
        }); // end FileDrop.Listener

        this.optionPanel.addComponentWithLabel("Input trees file: ", panel2);

        JButton outputFileButton = new JButton("Choose File...");
        final JTextField outputFileNameText = new JTextField("not selected", 16);

        outputFileButton.addActionListener(ae -> {
            FileDialog dialog = new FileDialog(frame,
                    "Select output file...",
                    FileDialog.SAVE);

            dialog.setVisible(true);
            if (dialog.getFile() == null) {
                // the dialog was cancelled...
                return;
            }

            outputFile = new File(dialog.getDirectory(), dialog.getFile());
            outputFileNameText.setText(outputFile.getName());

        });
        outputFileNameText.setEditable(false);

        JPanel panel3 = new JPanel(new BorderLayout(0, 0));
        panel3.add(outputFileNameText, BorderLayout.CENTER);
        panel3.add(outputFileButton, BorderLayout.EAST);
        this.optionPanel.addComponentWithLabel("Output tree file: ", panel3);

        this.optionPanel.addComponentWithLabel("Low memory: ", lowMemCheckbox);
        this.optionPanel.addComponentWithLabel("Save simple tree: ", simpleTreeCheckBox);
        this.optionPanel.addComponentWithLabel("Has trunk: ", hasTrunkCheckBox);
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

    public int getBurninPercentage() {
        return burninText.getValue();
    }

    public double getPosteriorLimit() {
        return limitText.getValue();
    }

    public ScsTreeAnnotator.Target getTargetOption() {
        return (ScsTreeAnnotator.Target) summaryTreeCombo.getSelectedItem();
    }

    public ScsTreeAnnotator.HeightsSummary getHeightsOption() {
        return (ScsTreeAnnotator.HeightsSummary) nodeHeightsCombo.getSelectedItem();
    }

    public String getTargetFileName() {
        if (targetFile == null) return null;
        return targetFile.getPath();
    }

    public String getInputFileName() {
        if (inputFile == null) return null;
        return inputFile.getPath();
    }

    public String getOutputFileName() {
        if (outputFile == null) return null;
        return outputFile.getPath();
    }

    public boolean useLowMem() {
        return lowMemCheckbox.isSelected();
    }

    public boolean saveSimpleTree() {
        return simpleTreeCheckBox.isSelected();
    }

    public boolean hasTrunk() {
        return hasTrunkCheckBox.isSelected();
    }

}
