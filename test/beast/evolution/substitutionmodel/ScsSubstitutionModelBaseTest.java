package beast.evolution.substitutionmodel;

import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class ScsSubstitutionModelBaseTest {

    @Test
    public void testGetEvolutionaryEventsMu() {
        ScsSubstitutionModelBase substModel = new ScsFiniteMuModel();


        // arrange
        ScsSubstitutionModelBase.EvolutionaryEventType[] z2oExpected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] z2tExpected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] o2zExpected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] o2tExpected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] t2zExpected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] t2oExpected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};


        // act
        ScsSubstitutionModelBase.EvolutionaryEventType[] z2zActual = substModel.getEvolutionaryEvents(0, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] z2oActual = substModel.getEvolutionaryEvents(0, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] z2tActual = substModel.getEvolutionaryEvents(0, 2);

        ScsSubstitutionModelBase.EvolutionaryEventType[] o2zActual = substModel.getEvolutionaryEvents(1, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] o2oActual = substModel.getEvolutionaryEvents(1, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] o2tActual = substModel.getEvolutionaryEvents(1, 2);

        ScsSubstitutionModelBase.EvolutionaryEventType[] t2zActual = substModel.getEvolutionaryEvents(2, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] t2oActual = substModel.getEvolutionaryEvents(2, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] t2tActual = substModel.getEvolutionaryEvents(2, 2);


        // assert
        assertNull(z2zActual);
        testHelper(z2oActual, z2oExpected);
        testHelper(z2tActual, z2tExpected);

        testHelper(o2zActual, o2zExpected);
        assertNull(o2oActual);
        testHelper(o2tActual, o2tExpected);

        testHelper(t2zActual, t2zExpected);
        testHelper(t2oActual, t2oExpected);
        assertNull(t2tActual);
    }

    @Test
    public void testGetEvolutionaryEventsMuDel() {
        ScsSubstitutionModelBase substModel = new ScsFiniteMuDelModel();


        // arrange
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};


        // act
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g0Actual = substModel.getEvolutionaryEvents(0, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g1Actual = substModel.getEvolutionaryEvents(0, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g2Actual = substModel.getEvolutionaryEvents(0, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g3Actual = substModel.getEvolutionaryEvents(0, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g4Actual = substModel.getEvolutionaryEvents(0, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g5Actual = substModel.getEvolutionaryEvents(0, 5);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g0Actual = substModel.getEvolutionaryEvents(1, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g1Actual = substModel.getEvolutionaryEvents(1, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g2Actual = substModel.getEvolutionaryEvents(1, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g3Actual = substModel.getEvolutionaryEvents(1, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g4Actual = substModel.getEvolutionaryEvents(1, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g5Actual = substModel.getEvolutionaryEvents(1, 5);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g0Actual = substModel.getEvolutionaryEvents(2, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g1Actual = substModel.getEvolutionaryEvents(2, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g2Actual = substModel.getEvolutionaryEvents(2, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g3Actual = substModel.getEvolutionaryEvents(2, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g4Actual = substModel.getEvolutionaryEvents(2, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g5Actual = substModel.getEvolutionaryEvents(2, 5);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g3Actual = substModel.getEvolutionaryEvents(3, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g4Actual = substModel.getEvolutionaryEvents(3, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g5Actual = substModel.getEvolutionaryEvents(3, 5);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g3Actual = substModel.getEvolutionaryEvents(4, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g4Actual = substModel.getEvolutionaryEvents(4, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g5Actual = substModel.getEvolutionaryEvents(4, 5);


        // assert
        assertNull(g02g0Actual);
        testHelper(g02g1Actual, g02g1Expected);
        testHelper(g02g2Actual, g02g2Expected);
        testHelper(g02g3Actual, g02g3Expected);
        testHelper(g02g4Actual, g02g4Expected);
        testHelper(g02g5Actual, g02g5Expected);

        testHelper(g12g0Actual, g12g0Expected);
        assertNull(g12g1Actual);
        testHelper(g12g2Actual, g12g2Expected);
        testHelper(g12g3Actual, g12g3Expected);
        testHelper(g12g4Actual, g12g4Expected);
        testHelper(g12g5Actual, g12g5Expected);

        testHelper(g22g0Actual, g22g0Expected);
        testHelper(g22g1Actual, g22g1Expected);
        assertNull(g22g2Actual);
        testHelper(g22g3Actual, g22g3Expected);
        testHelper(g22g4Actual, g22g4Expected);
        testHelper(g22g5Actual, g22g5Expected);

        try {
            substModel.getEvolutionaryEvents(3, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(3, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(3, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        assertNull(g32g3Actual);
        testHelper(g32g4Actual, g32g4Expected);
        testHelper(g32g5Actual, g32g5Expected);

        try {
            substModel.getEvolutionaryEvents(4, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(4, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(4, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g42g3Actual, g42g3Expected);
        assertNull(g42g4Actual);
        testHelper(g42g5Actual, g42g5Expected);

        try {
            substModel.getEvolutionaryEvents(5, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(5, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(5, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(5, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(5, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        assertNull(substModel.getEvolutionaryEvents(5, 5));
    }

    @Test
    public void testGetEvolutionaryEventsMuDelIns() {
        ScsSubstitutionModelBase substModel = new ScsFiniteMuDelInsModel();


        // arrange
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g10Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g0Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g1Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g2Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g3Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g4Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g5Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g6Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g8Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION, ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g13Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};

        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g7Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g9Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION, ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g11Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.INSERTION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g12Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.SINGLE_BACK_MUTATION};
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g14Expected = new ScsSubstitutionModelBase.EvolutionaryEventType[]{ScsSubstitutionModelBase.EvolutionaryEventType.DELETION};


        // act
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g0Actual = substModel.getEvolutionaryEvents(0, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g1Actual = substModel.getEvolutionaryEvents(0, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g2Actual = substModel.getEvolutionaryEvents(0, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g3Actual = substModel.getEvolutionaryEvents(0, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g4Actual = substModel.getEvolutionaryEvents(0, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g5Actual = substModel.getEvolutionaryEvents(0, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g6Actual = substModel.getEvolutionaryEvents(0, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g7Actual = substModel.getEvolutionaryEvents(0, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g8Actual = substModel.getEvolutionaryEvents(0, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g9Actual = substModel.getEvolutionaryEvents(0, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g10Actual = substModel.getEvolutionaryEvents(0, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g11Actual = substModel.getEvolutionaryEvents(0, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g12Actual = substModel.getEvolutionaryEvents(0, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g13Actual = substModel.getEvolutionaryEvents(0, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g02g14Actual = substModel.getEvolutionaryEvents(0, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g0Actual = substModel.getEvolutionaryEvents(1, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g1Actual = substModel.getEvolutionaryEvents(1, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g2Actual = substModel.getEvolutionaryEvents(1, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g3Actual = substModel.getEvolutionaryEvents(1, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g4Actual = substModel.getEvolutionaryEvents(1, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g5Actual = substModel.getEvolutionaryEvents(1, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g6Actual = substModel.getEvolutionaryEvents(1, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g7Actual = substModel.getEvolutionaryEvents(1, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g8Actual = substModel.getEvolutionaryEvents(1, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g9Actual = substModel.getEvolutionaryEvents(1, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g10Actual = substModel.getEvolutionaryEvents(1, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g11Actual = substModel.getEvolutionaryEvents(1, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g12Actual = substModel.getEvolutionaryEvents(1, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g13Actual = substModel.getEvolutionaryEvents(1, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g12g14Actual = substModel.getEvolutionaryEvents(1, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g0Actual = substModel.getEvolutionaryEvents(2, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g1Actual = substModel.getEvolutionaryEvents(2, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g2Actual = substModel.getEvolutionaryEvents(2, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g3Actual = substModel.getEvolutionaryEvents(2, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g4Actual = substModel.getEvolutionaryEvents(2, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g5Actual = substModel.getEvolutionaryEvents(2, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g6Actual = substModel.getEvolutionaryEvents(2, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g7Actual = substModel.getEvolutionaryEvents(2, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g8Actual = substModel.getEvolutionaryEvents(2, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g9Actual = substModel.getEvolutionaryEvents(2, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g10Actual = substModel.getEvolutionaryEvents(2, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g11Actual = substModel.getEvolutionaryEvents(2, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g12Actual = substModel.getEvolutionaryEvents(2, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g13Actual = substModel.getEvolutionaryEvents(2, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g22g14Actual = substModel.getEvolutionaryEvents(2, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g0Actual = substModel.getEvolutionaryEvents(3, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g1Actual = substModel.getEvolutionaryEvents(3, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g2Actual = substModel.getEvolutionaryEvents(3, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g3Actual = substModel.getEvolutionaryEvents(3, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g4Actual = substModel.getEvolutionaryEvents(3, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g5Actual = substModel.getEvolutionaryEvents(3, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g6Actual = substModel.getEvolutionaryEvents(3, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g7Actual = substModel.getEvolutionaryEvents(3, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g8Actual = substModel.getEvolutionaryEvents(3, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g9Actual = substModel.getEvolutionaryEvents(3, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g10Actual = substModel.getEvolutionaryEvents(3, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g11Actual = substModel.getEvolutionaryEvents(3, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g12Actual = substModel.getEvolutionaryEvents(3, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g13Actual = substModel.getEvolutionaryEvents(3, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g32g14Actual = substModel.getEvolutionaryEvents(3, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g0Actual = substModel.getEvolutionaryEvents(4, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g1Actual = substModel.getEvolutionaryEvents(4, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g2Actual = substModel.getEvolutionaryEvents(4, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g3Actual = substModel.getEvolutionaryEvents(4, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g4Actual = substModel.getEvolutionaryEvents(4, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g5Actual = substModel.getEvolutionaryEvents(4, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g6Actual = substModel.getEvolutionaryEvents(4, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g7Actual = substModel.getEvolutionaryEvents(4, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g8Actual = substModel.getEvolutionaryEvents(4, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g9Actual = substModel.getEvolutionaryEvents(4, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g10Actual = substModel.getEvolutionaryEvents(4, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g11Actual = substModel.getEvolutionaryEvents(4, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g12Actual = substModel.getEvolutionaryEvents(4, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g13Actual = substModel.getEvolutionaryEvents(4, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g42g14Actual = substModel.getEvolutionaryEvents(4, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g0Actual = substModel.getEvolutionaryEvents(5, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g1Actual = substModel.getEvolutionaryEvents(5, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g2Actual = substModel.getEvolutionaryEvents(5, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g3Actual = substModel.getEvolutionaryEvents(5, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g4Actual = substModel.getEvolutionaryEvents(5, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g5Actual = substModel.getEvolutionaryEvents(5, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g6Actual = substModel.getEvolutionaryEvents(5, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g7Actual = substModel.getEvolutionaryEvents(5, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g8Actual = substModel.getEvolutionaryEvents(5, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g9Actual = substModel.getEvolutionaryEvents(5, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g10Actual = substModel.getEvolutionaryEvents(5, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g11Actual = substModel.getEvolutionaryEvents(5, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g12Actual = substModel.getEvolutionaryEvents(5, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g13Actual = substModel.getEvolutionaryEvents(5, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g52g14Actual = substModel.getEvolutionaryEvents(5, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g0Actual = substModel.getEvolutionaryEvents(6, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g1Actual = substModel.getEvolutionaryEvents(6, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g2Actual = substModel.getEvolutionaryEvents(6, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g3Actual = substModel.getEvolutionaryEvents(6, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g4Actual = substModel.getEvolutionaryEvents(6, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g5Actual = substModel.getEvolutionaryEvents(6, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g6Actual = substModel.getEvolutionaryEvents(6, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g7Actual = substModel.getEvolutionaryEvents(6, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g8Actual = substModel.getEvolutionaryEvents(6, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g9Actual = substModel.getEvolutionaryEvents(6, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g10Actual = substModel.getEvolutionaryEvents(6, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g11Actual = substModel.getEvolutionaryEvents(6, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g12Actual = substModel.getEvolutionaryEvents(6, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g13Actual = substModel.getEvolutionaryEvents(6, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g62g14Actual = substModel.getEvolutionaryEvents(6, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g7Actual = substModel.getEvolutionaryEvents(7, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g9Actual = substModel.getEvolutionaryEvents(7, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g11Actual = substModel.getEvolutionaryEvents(7, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g12Actual = substModel.getEvolutionaryEvents(7, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g13Actual = substModel.getEvolutionaryEvents(7, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g72g14Actual = substModel.getEvolutionaryEvents(7, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g0Actual = substModel.getEvolutionaryEvents(8, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g1Actual = substModel.getEvolutionaryEvents(8, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g2Actual = substModel.getEvolutionaryEvents(8, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g3Actual = substModel.getEvolutionaryEvents(8, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g4Actual = substModel.getEvolutionaryEvents(8, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g5Actual = substModel.getEvolutionaryEvents(8, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g6Actual = substModel.getEvolutionaryEvents(8, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g7Actual = substModel.getEvolutionaryEvents(8, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g8Actual = substModel.getEvolutionaryEvents(8, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g9Actual = substModel.getEvolutionaryEvents(8, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g10Actual = substModel.getEvolutionaryEvents(8, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g11Actual = substModel.getEvolutionaryEvents(8, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g12Actual = substModel.getEvolutionaryEvents(8, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g13Actual = substModel.getEvolutionaryEvents(8, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g82g14Actual = substModel.getEvolutionaryEvents(8, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g7Actual = substModel.getEvolutionaryEvents(9, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g9Actual = substModel.getEvolutionaryEvents(9, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g11Actual = substModel.getEvolutionaryEvents(9, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g12Actual = substModel.getEvolutionaryEvents(9, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g13Actual = substModel.getEvolutionaryEvents(9, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g92g14Actual = substModel.getEvolutionaryEvents(9, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g0Actual = substModel.getEvolutionaryEvents(10, 0);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g1Actual = substModel.getEvolutionaryEvents(10, 1);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g2Actual = substModel.getEvolutionaryEvents(10, 2);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g3Actual = substModel.getEvolutionaryEvents(10, 3);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g4Actual = substModel.getEvolutionaryEvents(10, 4);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g5Actual = substModel.getEvolutionaryEvents(10, 5);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g6Actual = substModel.getEvolutionaryEvents(10, 6);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g7Actual = substModel.getEvolutionaryEvents(10, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g8Actual = substModel.getEvolutionaryEvents(10, 8);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g9Actual = substModel.getEvolutionaryEvents(10, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g10Actual = substModel.getEvolutionaryEvents(10, 10);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g11Actual = substModel.getEvolutionaryEvents(10, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g12Actual = substModel.getEvolutionaryEvents(10, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g13Actual = substModel.getEvolutionaryEvents(10, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g102g14Actual = substModel.getEvolutionaryEvents(10, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g7Actual = substModel.getEvolutionaryEvents(11, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g9Actual = substModel.getEvolutionaryEvents(11, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g11Actual = substModel.getEvolutionaryEvents(11, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g12Actual = substModel.getEvolutionaryEvents(11, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g13Actual = substModel.getEvolutionaryEvents(11, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g112g14Actual = substModel.getEvolutionaryEvents(11, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g7Actual = substModel.getEvolutionaryEvents(12, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g9Actual = substModel.getEvolutionaryEvents(12, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g11Actual = substModel.getEvolutionaryEvents(12, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g12Actual = substModel.getEvolutionaryEvents(12, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g13Actual = substModel.getEvolutionaryEvents(12, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g122g14Actual = substModel.getEvolutionaryEvents(12, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g7Actual = substModel.getEvolutionaryEvents(13, 7);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g9Actual = substModel.getEvolutionaryEvents(13, 9);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g11Actual = substModel.getEvolutionaryEvents(13, 11);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g12Actual = substModel.getEvolutionaryEvents(13, 12);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g13Actual = substModel.getEvolutionaryEvents(13, 13);
        ScsSubstitutionModelBase.EvolutionaryEventType[] g132g14Actual = substModel.getEvolutionaryEvents(13, 14);

        ScsSubstitutionModelBase.EvolutionaryEventType[] g142g14Actual = substModel.getEvolutionaryEvents(14, 14);


        // assert
        assertNull(g02g0Actual);
        testHelper(g02g1Actual, g02g1Expected);
        testHelper(g02g2Actual, g02g2Expected);
        testHelper(g02g3Actual, g02g3Expected);
        testHelper(g02g4Actual, g02g4Expected);
        testHelper(g02g5Actual, g02g5Expected);
        testHelper(g02g6Actual, g02g6Expected);
        testHelper(g02g7Actual, g02g7Expected);
        testHelper(g02g8Actual, g02g8Expected);
        testHelper(g02g9Actual, g02g9Expected);
        testHelper(g02g10Actual, g02g10Expected);
        testHelper(g02g11Actual, g02g11Expected);
        testHelper(g02g12Actual, g02g12Expected);
        testHelper(g02g13Actual, g02g13Expected);
        testHelper(g02g14Actual, g02g14Expected);

        testHelper(g12g0Actual, g12g0Expected);
        assertNull(g12g1Actual);
        testHelper(g12g2Actual, g12g2Expected);
        testHelper(g12g3Actual, g12g3Expected);
        testHelper(g12g4Actual, g12g4Expected);
        testHelper(g12g5Actual, g12g5Expected);
        testHelper(g12g6Actual, g12g6Expected);
        testHelper(g12g7Actual, g12g7Expected);
        testHelper(g12g8Actual, g12g8Expected);
        testHelper(g12g9Actual, g12g9Expected);
        testHelper(g12g10Actual, g12g10Expected);
        testHelper(g12g11Actual, g12g11Expected);
        testHelper(g12g12Actual, g12g12Expected);
        testHelper(g12g13Actual, g12g13Expected);
        testHelper(g12g14Actual, g12g14Expected);

        testHelper(g22g0Actual, g22g0Expected);
        testHelper(g22g1Actual, g22g1Expected);
        assertNull(g22g2Actual);
        testHelper(g22g3Actual, g22g3Expected);
        testHelper(g22g4Actual, g22g4Expected);
        testHelper(g22g5Actual, g22g5Expected);
        testHelper(g22g6Actual, g22g6Expected);
        testHelper(g22g7Actual, g22g7Expected);
        testHelper(g22g8Actual, g22g8Expected);
        testHelper(g22g9Actual, g22g9Expected);
        testHelper(g22g10Actual, g22g10Expected);
        testHelper(g22g11Actual, g22g11Expected);
        testHelper(g22g12Actual, g22g12Expected);
        testHelper(g22g13Actual, g22g13Expected);
        testHelper(g22g14Actual, g22g14Expected);

        testHelper(g32g0Actual, g32g0Expected);
        testHelper(g32g1Actual, g32g1Expected);
        testHelper(g32g2Actual, g32g2Expected);
        assertNull(g32g3Actual);
        testHelper(g32g4Actual, g32g4Expected);
        testHelper(g32g5Actual, g32g5Expected);
        testHelper(g32g6Actual, g32g6Expected);
        testHelper(g32g7Actual, g32g7Expected);
        testHelper(g32g8Actual, g32g8Expected);
        testHelper(g32g9Actual, g32g9Expected);
        testHelper(g32g10Actual, g32g10Expected);
        testHelper(g32g11Actual, g32g11Expected);
        testHelper(g32g12Actual, g32g12Expected);
        testHelper(g32g13Actual, g32g13Expected);
        testHelper(g32g14Actual, g32g14Expected);

        testHelper(g42g0Actual, g42g0Expected);
        testHelper(g42g1Actual, g42g1Expected);
        testHelper(g42g2Actual, g42g2Expected);
        testHelper(g42g3Actual, g42g3Expected);
        assertNull(g42g4Actual);
        testHelper(g42g5Actual, g42g5Expected);
        testHelper(g42g6Actual, g42g6Expected);
        testHelper(g42g7Actual, g42g7Expected);
        testHelper(g42g8Actual, g42g8Expected);
        testHelper(g42g9Actual, g42g9Expected);
        testHelper(g42g10Actual, g42g10Expected);
        testHelper(g42g11Actual, g42g11Expected);
        testHelper(g42g12Actual, g42g12Expected);
        testHelper(g42g13Actual, g42g13Expected);
        testHelper(g42g14Actual, g42g14Expected);

        testHelper(g52g0Actual, g52g0Expected);
        testHelper(g52g1Actual, g52g1Expected);
        testHelper(g52g2Actual, g52g2Expected);
        testHelper(g52g3Actual, g52g3Expected);
        testHelper(g52g4Actual, g52g4Expected);
        assertNull(g52g5Actual);
        testHelper(g52g6Actual, g52g6Expected);
        testHelper(g52g7Actual, g52g7Expected);
        testHelper(g52g8Actual, g52g8Expected);
        testHelper(g52g9Actual, g52g9Expected);
        testHelper(g52g10Actual, g52g10Expected);
        testHelper(g52g11Actual, g52g11Expected);
        testHelper(g52g12Actual, g52g12Expected);
        testHelper(g52g13Actual, g52g13Expected);
        testHelper(g52g14Actual, g52g14Expected);

        testHelper(g62g0Actual, g62g0Expected);
        testHelper(g62g1Actual, g62g1Expected);
        testHelper(g62g2Actual, g62g2Expected);
        testHelper(g62g3Actual, g62g3Expected);
        testHelper(g62g4Actual, g62g4Expected);
        testHelper(g62g5Actual, g62g5Expected);
        assertNull(g62g6Actual);
        testHelper(g62g7Actual, g62g7Expected);
        testHelper(g62g8Actual, g62g8Expected);
        testHelper(g62g9Actual, g62g9Expected);
        testHelper(g62g10Actual, g62g10Expected);
        testHelper(g62g11Actual, g62g11Expected);
        testHelper(g62g12Actual, g62g12Expected);
        testHelper(g62g13Actual, g62g13Expected);
        testHelper(g62g14Actual, g62g14Expected);

        try {
            substModel.getEvolutionaryEvents(7, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(7, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(7, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(7, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(7, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(7, 5);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(7, 6);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        assertNull(g72g7Actual);
        try {
            substModel.getEvolutionaryEvents(7, 8);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g72g9Actual, g72g9Expected);
        try {
            substModel.getEvolutionaryEvents(7, 10);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g72g11Actual, g72g11Expected);
        testHelper(g72g12Actual, g72g12Expected);
        testHelper(g72g13Actual, g72g13Expected);
        testHelper(g72g14Actual, g72g14Expected);

        testHelper(g82g0Actual, g82g0Expected);
        testHelper(g82g1Actual, g82g1Expected);
        testHelper(g82g2Actual, g82g2Expected);
        testHelper(g82g3Actual, g82g3Expected);
        testHelper(g82g4Actual, g82g4Expected);
        testHelper(g82g5Actual, g82g5Expected);
        testHelper(g82g6Actual, g82g6Expected);
        testHelper(g82g7Actual, g82g7Expected);
        assertNull(g82g8Actual);
        testHelper(g82g9Actual, g82g9Expected);
        testHelper(g82g10Actual, g82g10Expected);
        testHelper(g82g11Actual, g82g11Expected);
        testHelper(g82g12Actual, g82g12Expected);
        testHelper(g82g13Actual, g82g13Expected);
        testHelper(g82g14Actual, g82g14Expected);

        try {
            substModel.getEvolutionaryEvents(9, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(9, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(9, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(9, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(9, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(9, 5);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(9, 6);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g92g7Actual, g92g7Expected);
        try {
            substModel.getEvolutionaryEvents(9, 8);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        assertNull(g92g9Actual);
        try {
            substModel.getEvolutionaryEvents(9, 10);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g92g11Actual, g92g11Expected);
        testHelper(g92g12Actual, g92g12Expected);
        testHelper(g92g13Actual, g92g13Expected);
        testHelper(g92g14Actual, g92g14Expected);

        testHelper(g102g0Actual, g102g0Expected);
        testHelper(g102g1Actual, g102g1Expected);
        testHelper(g102g2Actual, g102g2Expected);
        testHelper(g102g3Actual, g102g3Expected);
        testHelper(g102g4Actual, g102g4Expected);
        testHelper(g102g5Actual, g102g5Expected);
        testHelper(g102g6Actual, g102g6Expected);
        testHelper(g102g7Actual, g102g7Expected);
        testHelper(g102g8Actual, g102g8Expected);
        testHelper(g102g9Actual, g102g9Expected);
        assertNull(g102g10Actual);
        testHelper(g102g11Actual, g102g11Expected);
        testHelper(g102g12Actual, g102g12Expected);
        testHelper(g102g13Actual, g102g13Expected);
        testHelper(g102g14Actual, g102g14Expected);

        try {
            substModel.getEvolutionaryEvents(11, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(11, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(11, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(11, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(11, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(11, 5);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(11, 6);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g112g7Actual, g112g7Expected);
        try {
            substModel.getEvolutionaryEvents(11, 8);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g112g9Actual, g112g9Expected);
        try {
            substModel.getEvolutionaryEvents(11, 10);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        assertNull(g112g11Actual);
        testHelper(g112g12Actual, g112g12Expected);
        testHelper(g112g13Actual, g112g13Expected);
        testHelper(g112g14Actual, g112g14Expected);

        try {
            substModel.getEvolutionaryEvents(12, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(12, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(12, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(12, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(12, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(12, 5);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(12, 6);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g122g7Actual, g122g7Expected);
        try {
            substModel.getEvolutionaryEvents(12, 8);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g122g9Actual, g122g9Expected);
        try {
            substModel.getEvolutionaryEvents(12, 10);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g122g11Actual, g122g11Expected);
        assertNull(g122g12Actual);
        testHelper(g122g13Actual, g122g13Expected);
        testHelper(g122g14Actual, g122g14Expected);

        try {
            substModel.getEvolutionaryEvents(13, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(13, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(13, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(13, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(13, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(13, 5);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(13, 6);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g132g7Actual, g132g7Expected);
        try {
            substModel.getEvolutionaryEvents(13, 8);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g132g9Actual, g132g9Expected);
        try {
            substModel.getEvolutionaryEvents(13, 10);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        testHelper(g132g11Actual, g132g11Expected);
        testHelper(g132g12Actual, g132g12Expected);
        assertNull(g132g13Actual);
        testHelper(g132g14Actual, g132g14Expected);

        try {
            substModel.getEvolutionaryEvents(14, 0);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 1);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 2);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 3);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 4);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 5);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 6);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 7);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 8);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 9);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 10);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 11);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 12);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        try {
            substModel.getEvolutionaryEvents(14, 13);
            fail();
        } catch (IllegalStateException e) {
            //
        }
        assertNull(g142g14Actual);
    }

    public void testHelper(
            ScsSubstitutionModelBase.EvolutionaryEventType[] actual,
            ScsSubstitutionModelBase.EvolutionaryEventType[] expected
    ) {
        assertEquals(actual.length, expected.length);

        for (int i = 0; i < expected.length; i++)
            assertEquals(actual[i], expected[i]);
    }

}