package grok_connect.providers.utils;

import serialization.Column;
import serialization.DataFrame;

/**
 * DataFrame builder for reducing amount of code when generating data for tests
 */
public class DataFrameBuilder {
    private static DataFrameBuilder builder;
    private DataFrame dataFrame;

    private DataFrameBuilder() {
    }

    public static DataFrameBuilder getBuilder() {
        if (builder == null) {
            builder = new DataFrameBuilder();
            builder.restore();
        }
        return builder;
    }

    public DataFrameBuilder setRowCount(int rowCount) {
        dataFrame.rowCount = rowCount;
        return this;
    }

    public DataFrameBuilder setColumn(Column<?> column, String columnName) {
        column.name = columnName;
        dataFrame.addColumn(column);
        return this;
    }

    public <T> DataFrameBuilder setColumn(Column<T> column, String columnName, T[] data) {
        column.name = columnName;
        for (T t: data) {
            column.add(t);
        }
        dataFrame.addColumn(column);
        return this;
    }

    public DataFrame build() {
        DataFrame frame = dataFrame;
        restore();
        return frame;
    }

    public void restore() {
        dataFrame = new DataFrame();
    }
}
