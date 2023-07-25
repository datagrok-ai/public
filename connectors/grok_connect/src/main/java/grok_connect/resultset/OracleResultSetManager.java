package grok_connect.resultset;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.bigint_column.OracleBigIntColumnManager;
import grok_connect.managers.string_column.DefaultStringColumnManager;
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
    public void processValue(Object o, int index) {
        if (isInit) {
            Column currentColumn = columns[index - 1];
            ColumnManager<?> currentManager = currentManagers[index - 1];
            if (currentColumn.getType().equals(Types.INT) && currentManager.getClass().equals(OracleBigIntColumnManager.class)) {
                setBigIntValue(o, index, currentColumn);
            } else
                currentColumn.add(currentManager
                        .convert(o, columnsMeta[index - 1]));
        } else
            throw new RuntimeException("ResultSetManager should be init");
    }


    private void setBigIntValue(Object o, int index, Column column) {
        if (o == null)
            column.add(o);
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
