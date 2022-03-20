package beast.evolution.operators;

import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ScsSwapOperator extends SwapOperator {

    private List<Integer> masterList = null;

    @Override
    public void initAndValidate() {
        if (parameterInput.get() != null) {
            parameter = parameterInput.get();
        } else {
            parameter = intparameterInput.get();
        }

        howMany = howManyInput.get();
        if (howMany * 2 > parameter.getDimension()) {
            throw new IllegalArgumentException("howMany too large: must be less than half the parameter dimension");
        }

        filter = parameterFilterInput.get();
        if (filter != null) {
            filter.initAndValidate();
            if (filter.getDimension() != parameter.getDimension())
                throw new IllegalArgumentException("Filter vector should have the same length as parameter");
        }

        List<Integer> list = new ArrayList<>();
        for (int i = 0; i < parameter.getDimension(); i++) {
            if (filter == null) {
                list.add(i);
            } else if (filter.getValue(i)) {
                list.add(i);
            }
        }
        masterList = Collections.unmodifiableList(list);
    }

    @Override
    public double proposal() {
        List<Integer> allIndices = new ArrayList<>(masterList);
        int left, right;

        if (Arrays.stream(parameter.getValues()).distinct().toArray().length <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        for (int i = 0; i < howMany; i++) {
            do {
                left = allIndices.remove(Randomizer.nextInt(allIndices.size()));
                right = allIndices.remove(Randomizer.nextInt(allIndices.size()));
            } while (parameter.getValue(left) == parameter.getValue(right));

            parameter.swap(left, right);
        }

        return 0.0;
    }

}
