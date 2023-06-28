package grok_connect.resultset;

public class ColumnMeta {
    private final int type;
    private final String typeName;
    private final int precision;
    private final int scale;
    private final String columnLabel;

    public ColumnMeta(int type, String typeName, int precision, int scale, String columnLabel) {
        this.type = type;
        this.typeName = typeName;
        this.precision = precision;
        this.scale = scale;
        this.columnLabel = columnLabel;
    }

    public int getType() {
        return type;
    }

    public String getTypeName() {
        return typeName;
    }

    public int getPrecision() {
        return precision;
    }

    public int getScale() {
        return scale;
    }

    public String getColumnLabel() {
        return columnLabel;
    }
}
