package grok_connect.managers.integer_column;

import grok_connect.resultset.ColumnMeta;

public class CassandraIntColumnManager extends DefaultIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        return !columnMeta.getTypeName().equalsIgnoreCase("varint")
                && ((type == java.sql.Types.INTEGER) || (type == java.sql.Types.TINYINT) || (type == java.sql.Types.SMALLINT));
    }

    // The driver resolves a codec per CQL type and refuses the widening the typed getters ask for
    // ("Codec not found ... [TINYINT <-> java.lang.Integer]"), so this provider reads through
    // getObject and the converter chain instead.
    @Override
    public boolean canReadFast(ColumnMeta columnMeta) {
        return false;
    }
}
