package grok_connect.providers.utils;

import serialization.Column;
import serialization.DataFrame;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class DataFrameComparator {
    public boolean isDataFramesEqual(DataFrame that, DataFrame other) {
        return that.getClass().equals(other.getClass())
                && that.rowCount.equals(other.rowCount) && isColumnsEqual(that.getColumns(), other.getColumns());
    }

    public boolean isDataFramesEqualUnOrdered(DataFrame that, DataFrame other) {
        return that.getClass().equals(other.getClass())
                && that.rowCount.equals(other.rowCount) && isColumnsEqualUnOrdered(that.getColumns(), other.getColumns());
    }

    public boolean isColumnsEqualUnOrdered(List<Column<?>> columns, List<Column<?>> columns1) {
        if (columns.size() != columns1.size())
            return false;
        for (Column<?> column: columns) {
            List<Column<?>> filtered = columns1.stream()
                    .filter(col -> col.getName().equals(column.getName()))
                    .collect(Collectors.toList());
            if (filtered.size() == 0)
                return false;
            if (!isTwoColumnsEqual(column, filtered.get(0)))
                return false;
        }
        return true;
    }

    public boolean isColumnsEqual(List<Column<?>> columns, List<Column<?>> columns1) {
        if (columns.size() != columns1.size())
            return false;
        for (int i = 0; i < columns.size(); i++)
            if (!isTwoColumnsEqual(columns.get(i), columns1.get(i)))
                return false;
        return true;
    }

    public boolean isTwoColumnsEqual(Column<?> column, Column<?> column1) {
        if (!column.getType().equals(column1.getType()))
            return false;
        if (!column.getName().equals(column1.getName()))
            return false;
        return arraysEqual(column.toArray(), column1.toArray());
    }

    private boolean arraysEqual(Object a, Object b) {
        if (a instanceof int[] && b instanceof int[])
            return Arrays.equals((int[]) a, (int[]) b);
        if (a instanceof float[] && b instanceof float[])
            return Arrays.equals((float[]) a, (float[]) b);
        if (a instanceof double[] && b instanceof double[])
            return Arrays.equals((double[]) a, (double[]) b);
        if (a instanceof Object[] && b instanceof Object[])
            return Arrays.deepEquals((Object[]) a, (Object[]) b);
        return false;
    }
}
