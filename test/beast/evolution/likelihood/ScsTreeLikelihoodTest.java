package beast.evolution.likelihood;

import beast.evolution.tree.Node;
import beast.evolution.tree.ScsNode;
import beast.evolution.tree.ScsTree;
import beast.evolution.tree.Tree;
import beast.util.ElementComparator;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class ScsTreeLikelihoodTest {

    private Tree tree;
    private Node targetNode;
    private List<int[]> MLGenotypes;
    private ScsGenericTreeLikelihood treeLikelihood;

    private Node newNode() {
        try {
            return ScsNode.class.newInstance();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }

    @SuppressWarnings("deprecation")
    private void makeTree() {
        tree = new ScsTree();

        Node n17 = newNode();
        Node n16 = newNode();
        assert n17 != null;
        n17.setNr(17);
        n17.setLeft(n16);
        assert n16 != null;
        n16.setNr(16);
        n16.setParent(n17);

        tree.setRoot(n17);

        Node n15 = newNode();
        Node n8 = newNode();
        n16.setLeft(n15);
        n16.setRight(n8);
        assert n15 != null;
        n15.setNr(15);
        n15.setParent(n16);
        assert n8 != null;
        n8.setNr(8);
        n8.setParent(n16);

        Node n9 = newNode();
        Node n6 = newNode();
        n15.setLeft(n9);
        n15.setRight(n6);
        assert n9 != null;
        n9.setNr(9);
        n9.setParent(n15);
        assert n6 != null;
        n6.setNr(6);
        n6.setParent(n15);

        targetNode = n9;

        Node n10 = newNode();
        Node n12 = newNode();
        n9.setLeft(n10);
        n9.setRight(n12);
        assert n10 != null;
        n10.setNr(10);
        n10.setParent(n9);
        assert n12 != null;
        n12.setNr(12);
        n12.setParent(n9);

        Node n1 = newNode();
        Node n13 = newNode();
        n10.setLeft(n1);
        n10.setRight(n13);
        assert n1 != null;
        n1.setNr(1);
        n1.setParent(n10);
        assert n13 != null;
        n13.setParent(n10);
        n13.setNr(13);

        Node n3 = newNode();
        Node n7 = newNode();
        n13.setLeft(n3);
        n13.setRight(n7);
        assert n3 != null;
        n3.setNr(3);
        n3.setParent(n13);
        assert n7 != null;
        n7.setNr(7);
        n7.setParent(n13);

        Node n5 = newNode();
        Node n14 = newNode();
        n12.setLeft(n5);
        n12.setRight(n14);
        assert n5 != null;
        n5.setNr(5);
        n5.setParent(n12);
        assert n14 != null;
        n14.setNr(14);
        n14.setParent(n12);

        Node n4 = newNode();
        Node n11 = newNode();
        n14.setLeft(n4);
        n14.setRight(n11);
        assert n4 != null;
        n4.setNr(4);
        n4.setParent(n14);
        assert n11 != null;
        n11.setNr(11);
        n11.setParent(n14);

        Node n2 = newNode();
        Node n0 = newNode();
        n11.setLeft(n2);
        n11.setRight(n0);
        assert n2 != null;
        n2.setNr(2);
        n2.setParent(n11);
        assert n0 != null;
        n0.setNr(0);
        n0.setParent(n11);
    }

    @BeforeMethod
    private void setUp() {
        makeTree();
        MLGenotypes = new ArrayList<>();
        final int nrOfNodes = 18;
        treeLikelihood = new ScsTreeLikelihood(nrOfNodes);
    }

    @Test
    public void testUpdateMaxLikelihoodGenotypes1() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {

        // arrange
        final int[] e1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] e2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        final int[] e3 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0};
        final int[] e4 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};

        // act
        final int[] g1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] g2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        MLGenotypes.add(g1);
        MLGenotypes.add(g2);

        Method testedMethod = ScsTreeLikelihood.class.getDeclaredMethod("updateMaxLikelihoodGenotypes", List.class, Node.class, int.class);
        testedMethod.setAccessible(true);
        testedMethod.invoke(treeLikelihood, MLGenotypes, targetNode, 0);

        // assert
        assertEquals(4, MLGenotypes.size());
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e1));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e2));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e3));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e4));
    }

    @Test
    public void testUpdateMaxLikelihoodGenotypes2() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {

        // arrange
        final int[] e1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] e2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        final int[] e3 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0};
        final int[] e4 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
        final int[] e5 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0};
        final int[] e6 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};
        final int[] e7 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0};
        final int[] e8 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0};

        // act
        final int[] g1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] g2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        final int[] g3 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0};
        final int[] g4 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};
        MLGenotypes.add(g1);
        MLGenotypes.add(g2);
        MLGenotypes.add(g3);
        MLGenotypes.add(g4);

        Method testedMethod = ScsTreeLikelihood.class.getDeclaredMethod("updateMaxLikelihoodGenotypes", List.class, Node.class, int.class);
        testedMethod.setAccessible(true);
        testedMethod.invoke(treeLikelihood, MLGenotypes, targetNode, 0);

        // assert
        assertEquals(8, MLGenotypes.size());
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e1));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e2));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e3));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e4));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e5));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e6));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e7));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e8));
    }

    @Test
    public void testUpdateMaxLikelihoodGenotypes3() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {

        // arrange
        final int[] e1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] e2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        final int[] e3 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0};
        final int[] e4 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
        final int[] e5 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0};
        final int[] e6 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};
        final int[] e7 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0};
        final int[] e8 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0};

        // act
        final int[] g1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] g2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0};
        MLGenotypes.add(g1);
        MLGenotypes.add(g2);

        Method testedMethod = ScsTreeLikelihood.class.getDeclaredMethod("updateMaxLikelihoodGenotypes", List.class, Node.class, int.class);
        testedMethod.setAccessible(true);
        testedMethod.invoke(treeLikelihood, MLGenotypes, targetNode, 0);

        // assert
        assertEquals(8, MLGenotypes.size());
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e1));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e2));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e3));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e4));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e5));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e6));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e7));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e8));
    }

    @Test
    public void testUpdateMaxLikelihoodGenotypes4() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {

        // arrange
        final int[] e1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] e2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        final int[] e3 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0};
        final int[] e4 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0};

        final int[] e5 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] e6 = {1, 1, 1, 0, 1, 1, 1, 2, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0};

        // act
        final int[] g1 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] g2 = {0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
        final int[] g3 = {1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0};
        final int[] g4 = {1, 1, 1, 0, 1, 1, 1, 2, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0};
        MLGenotypes.add(g1);
        MLGenotypes.add(g2);
        MLGenotypes.add(g3);
        MLGenotypes.add(g4);

        Set<Integer> genotypes = new HashSet<>();
        genotypes.add(0);
        genotypes.add(1);

        Method testedMethod = ScsTreeLikelihood.class.getDeclaredMethod("updateMaxLikelihoodGenotypes", List.class, Node.class, Set.class);
        testedMethod.setAccessible(true);
        testedMethod.invoke(treeLikelihood, MLGenotypes, targetNode, genotypes);

        // assert
        assertEquals(6, MLGenotypes.size());
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e1));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e2));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e3));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e4));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e5));
        assertTrue(ElementComparator.compIntArrays(MLGenotypes, e6));
    }

    @AfterMethod
    private void tearDown() {
        tree = null;
        MLGenotypes = null;
    }

}