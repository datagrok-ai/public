package grok_connect.managers.float_column;

import grok_connect.resultset.ColumnMeta;

public class CassandraFloatColumnManager extends DefaultFloatColumnManager {
    // The driver resolves a codec per CQL type and refuses the widening the typed getters ask for
    // ("Codec not found ... [FLOAT <-> java.lang.Double]"), so this provider reads through
    // getObject and the converter chain instead.
    @Override
    public boolean canReadFast(ColumnMeta columnMeta) {
        return false;
    }
}
