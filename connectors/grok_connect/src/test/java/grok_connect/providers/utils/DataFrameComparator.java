package grok_connect.providers.utils;

import serialization.BoolColumn;
import serialization.Column;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;
import java.util.Arrays;
import java.util.List;

public class DataFrameComparator {
    public boolean isDataFramesEqual(DataFrame that, DataFrame other) {
        return that.getClass().equals(other.getClass())
                && that.rowCount.equals(other.rowCount) && isColumnsEqual(that.columns, other.columns);

    }

    private boolean isColumnsEqual(List<Column> columns, List<Column> columns1) {
        if (columns.size() != columns1.size()) {
            return false;
        }
        for (int i = 0; i < columns.size(); i++) {
            Column column = columns.get(i);
            Column column1 = columns1.get(i);
            if (!isTwoColumnsEqual(column, column1)) {
                return false;
            }
        }
        return true;
    }

    private boolean isTwoColumnsEqual(Column<?> column, Column<?> column1) {
        if (!column.getType().equals(column1.getType())) {
            return false;
        }
        switch (column.getType()) {
            case Types.STRING:
            case Types.BIG_INT:
                StringColumn sc = (StringColumn) column;
                StringColumn sc1 = (StringColumn) column1;
                return sc.name.equals(sc1.name) && Arrays.equals(sc.getData(), sc1.getData());
            case Types.INT:
                IntColumn ic = (IntColumn) column;
                IntColumn ic1 = (IntColumn) column1;
                return ic.name.equals(ic1.name) && Arrays.equals(ic.getData(), ic1.getData());
            case Types.FLOAT:
                FloatColumn fc = (FloatColumn) column;
                FloatColumn fc1 = (FloatColumn) column1;
                return fc.name.equals(fc1.name) && Arrays.equals(fc.getData(), fc1.getData());
            case Types.DATE_TIME:
                DateTimeColumn dtc = (DateTimeColumn) column;
                DateTimeColumn dtc1 = (DateTimeColumn) column1;
                return dtc.name.equals(dtc1.name) && Arrays.equals(dtc.getData(), dtc1.getData());
            case Types.BOOL:
                BoolColumn bc = (BoolColumn) column;
                BoolColumn bc1 = (BoolColumn) column1;
                return bc.name.equals(bc1.name) && Arrays.equals(bc.getData(), bc1.getData());
            default:
                return false;
        }
    }
}
