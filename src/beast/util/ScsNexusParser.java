package beast.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

public class ScsNexusParser extends NexusParser {

    protected String nodeType;

    public ScsNexusParser(final String nodeType) {
        this.nodeType = nodeType;
    }

    @Override
    protected void parseTreesBlock(final BufferedReader fin) throws IOException {
        trees = new ArrayList<>();
        // read to first command within trees block
        NexusCommand nextCommand = readNextCommand(fin);

        int origin = -1;

        // if first non-empty line is "translate" then parse translate block
        if (nextCommand.isCommand("translate")) {
            translationMap = parseTranslateCommand(nextCommand.arguments);
            origin = getIndexedTranslationMapOrigin(translationMap);
            if (origin != -1) {
                taxa = getIndexedTranslationMap(translationMap, origin);
            }
        }

        // read trees
        while (nextCommand != null && !nextCommand.isEndOfBlock()) {
            if (nextCommand.isCommand("tree")) {
                String treeString = nextCommand.arguments;
                final int i = treeString.indexOf('(');
                if (i > 0) {
                    treeString = treeString.substring(i);
                }
                TreeParser treeParser;

                if (origin != -1) {
                    treeParser = new ScsTreeParser(taxa, treeString, origin, false, nodeType);
                } else {
                    try {
                        treeParser = new ScsTreeParser(taxa, treeString, 0, false, nodeType);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        treeParser = new ScsTreeParser(taxa, treeString, 1, false, nodeType);
                    }
                }

                // this needs to go after translation map or listeners have an incomplete tree!
                for (final NexusParserListener listener : listeners) {
                    listener.treeParsed(trees.size(), treeParser);
                }

                // this must come after listener or trees.size() gives the wrong index to treeParsed
                trees.add(treeParser);

            }
            nextCommand = readNextCommand(fin);
        }
    }

}
