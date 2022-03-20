package beast.app.postprocessor;

import beast.app.treeannotator.ScsTreeAnnotatorDialog;
import beast.app.util.WholeNumberField;
import beast.app.variantcaller.VariantCallerDialog;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

public class PostProcessorDialog {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    private JFrame frame;

    private final OptionsPanel optionPanel;

    private final WholeNumberField burninText = new WholeNumberField(0, Integer.MAX_VALUE);

    private ScsTreeAnnotatorDialog TADialog;

    private VariantCallerDialog VCDialog;

    //***********************************************
    //*                 Constructor                 *
    //***********************************************

    public PostProcessorDialog(final JFrame frame) {
        this.frame = frame;
        this.optionPanel = new OptionsPanel(100, 100);

        burninText.setColumns(12);
        burninText.setValue(0);
        this.optionPanel.addComponentWithLabel("Global burnin percentage: ", burninText);

        OptionsPanel TAPanel = getPanel(Color.GRAY);
        TADialog = new ScsTreeAnnotatorDialog(true, frame, TAPanel);
        this.optionPanel.addSpanningComponent(TAPanel);

        OptionsPanel VCPanel = getPanel(Color.GRAY);
        VCDialog = new VariantCallerDialog(true, frame, VCPanel);
        this.optionPanel.addSpanningComponent(VCPanel);
    }


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    private OptionsPanel getPanel(Color c) {
        OptionsPanel panel = new OptionsPanel();
        panel.setBorder(BorderFactory.createLineBorder(c));
        return panel;
    } // getPanel

    public boolean showDialog(String title) {

        JOptionPane optionPane = new JOptionPane(
                optionPanel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                new String[]{"Run", "Quit"},
                null
        );
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, title);
        dialog.setResizable(true);
        dialog.pack();

        dialog.setVisible(true);

        return optionPane.getValue().equals("Run");
    }

    public int getBurninText() {
        return burninText.getValue();
    }

    public ScsTreeAnnotatorDialog getTADialog() {
        return TADialog;
    }

    public VariantCallerDialog getVCDialog() {
        return VCDialog;
    }

}
