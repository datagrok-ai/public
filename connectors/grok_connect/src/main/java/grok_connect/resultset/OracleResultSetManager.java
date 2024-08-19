package grok_connect.resultset;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.bigint_column.OracleBigIntColumnManager;
import serialization.BigIntColumn;
import serialization.Column;
import serialization.IntColumn;
import serialization.Types;
import java.math.BigInteger;
import java.sql.ResultSetMetaData;
import java.util.Collection;

public class OracleResultSetManager extends DefaultResultSetManager {
    public OracleResultSetManager(Collection<ColumnManager<?>> columnManagers) {
        super(columnManagers);
    }

    @Override
    public void init(ResultSetMetaData meta, int initColumnSize) {
        super.init(meta, initColumnSize);
        for (int i = 0; i < columns.length; i++) {
            Column column = columns[i];
            if (column.getType().equals(Types.BIG_INT)) {
                columns[i] = new IntColumn(initColumnSize);
                columns[i].name = column.name;
            }
        }
    }

    @Override
    protected void setValue(Object o, int index) {
        Column column = columns[index - 1];
        if (column.getType().equals(Types.INT) && currentManagers[index - 1].getClass().equals(OracleBigIntColumnManager.class))
            setBigIntValue(o, index, column);
        else
            super.setValue(o, index);
    }

    private void setBigIntValue(Object o, int index, Column column) {
        if (o == null)
            column.add(null);
        else {
            String str = o.toString();
            BigInteger bigIntValue = new BigInteger(str);
            if (bigIntValue.compareTo(BigInteger.valueOf(Integer.MAX_VALUE)) <= 0 &&
                    bigIntValue.compareTo(BigInteger.valueOf(Integer.MIN_VALUE)) >= 0) {
                column.add(bigIntValue.intValue());
            } else {
                Column bigIntColumn = new BigIntColumn();
                for (int i = 0; i < column.length; i++) {
                    Object currentValue = column.get(i);
                    bigIntColumn.add(currentValue.equals(IntColumn.None) ? "" : currentValue.toString());
                }
                bigIntColumn.add(str);
                bigIntColumn.name = column.name;
                columns[index - 1] = bigIntColumn;
            }
        }
    }
}
