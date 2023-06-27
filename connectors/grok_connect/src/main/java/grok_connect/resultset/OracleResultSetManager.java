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
    public void processValue(Object o, int index, ResultSetMetaData meta) {
        boolean inBounds = (index - 1 >= 0) && (index - 1 < columns.size());
        if (!inBounds) {
            processHeaders(o, index, meta);
        } else {
            Column currentColumn = columns.get(index - 1);
            ColumnManager<?> currentManager = currentManagers.get(index - 1);
            if (currentManager.getClass().equals(OracleBigIntColumnManager.class) && currentColumn.getType().equals(Types.INT)) {
                if (o == null)
                    currentColumn.add(o);
                else {
                    String str = o.toString();
                    BigInteger bigIntValue = new BigInteger(str);
                    if (bigIntValue.compareTo(BigInteger.valueOf(Integer.MAX_VALUE)) <= 0 &&
                            bigIntValue.compareTo(BigInteger.valueOf(Integer.MIN_VALUE)) >= 0) {
                        currentColumn.add(bigIntValue.intValue());
                    } else {
                        Column bigIntColumn = new BigIntColumn();
                        for (int i = 0; i < currentColumn.length; i++) {
                            Object currentValue = currentColumn.get(i);
                            bigIntColumn.add(currentValue.equals(IntColumn.None) ? "" : currentValue.toString());
                        }
                        bigIntColumn.add(str);
                        bigIntColumn.name = currentColumn.name;
                        columns.set(index - 1, bigIntColumn);
                    }
                }
            } else
                currentColumn.add(currentManager
                        .convert(o, columnsMeta.get(index - 1).getColumnLabel()));
        }
    }
}
