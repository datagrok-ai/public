package grok_connect.type;

import java.sql.Types;

public class DefaultTypeCheckers {
    public static final TypeChecker DEFAULT_BIGINT_TYPECHECKER = (type, typeName, precision, scale) ->
            type == java.sql.Types.BIGINT
                    || typeName.equalsIgnoreCase("int8")
                    || typeName.equalsIgnoreCase("serial8");
    public static final TypeChecker DEFAULT_BITSTRING_TYPECHECKER = (type, typeName, precision, scale) ->
            type == java.sql.Types.BIT && precision > 1 || (scale > 1 && type == java.sql.Types.BIT)
                    || typeName.equalsIgnoreCase("varbit");
    public static final TypeChecker DEFAULT_BOOL_TYPECHECKER = (type, typeName, precision, scale) ->
            (type == java.sql.Types.BOOLEAN) ||
                    typeName.equalsIgnoreCase("bool") || (type == java.sql.Types.BIT
                    && precision == 1 && scale == 0);
    public static final TypeChecker DEFAULT_COMPLEX_TYPECHECKER = (type, typeName, precision, scale) -> false;
    public static final TypeChecker DEFAULT_DATETIME_TYPECHECKER = (type, typeName, precision, scale) ->
            (type == java.sql.Types.DATE) || (type == java.sql.Types.TIME) || (type == java.sql.Types.TIMESTAMP)
                    || type == java.sql.Types.TIMESTAMP_WITH_TIMEZONE
                    || type == java.sql.Types.TIME_WITH_TIMEZONE
                    || typeName.equalsIgnoreCase("timetz")
                    || typeName.equalsIgnoreCase("timestamptz")
                    || (typeName.equalsIgnoreCase("TIMESTAMP WITH TIME ZONE"))
                    || (typeName.equalsIgnoreCase("datetimeoffset"));
    public static final TypeChecker DEFAULT_FLOAT_TYPECHECKER = (type, typeName, precision, scale) ->
            type == Types.FLOAT || type == java.sql.Types.DOUBLE || type == java.sql.Types.REAL ||
                    type == Types.DECIMAL ||
                    typeName.equalsIgnoreCase("float8") ||
                    typeName.equalsIgnoreCase("float4") ||
                    typeName.equalsIgnoreCase("money") ||
                    typeName.equalsIgnoreCase("binary_float") ||
                    typeName.equalsIgnoreCase("binary_double") ||
                    typeName.equalsIgnoreCase("numeric") ||
                    typeName.equalsIgnoreCase("DECFLOAT") ||
                    (typeName.equalsIgnoreCase("number") && scale > 0);
    public static final TypeChecker DEFAULT_INT_TYPECHECKER = (type, typeName, precision, scale) ->
            (type == java.sql.Types.INTEGER) || (type == java.sql.Types.TINYINT) || (type == java.sql.Types.SMALLINT)
                    || typeName.equalsIgnoreCase("int4") || typeName.equalsIgnoreCase("int2")
                    || typeName.equalsIgnoreCase("int") || typeName.equalsIgnoreCase("serial2")
                    || typeName.equalsIgnoreCase("serial4") || typeName.equalsIgnoreCase("UInt16")
                    || typeName.equalsIgnoreCase("UInt8") || (typeName.equalsIgnoreCase("NUMBER")
                    && precision < 10 && scale == 0);
    public static final TypeChecker DEFAULT_ARRAY_TYPECHECKER = (type, typeName, precision, scale) ->
            type == java.sql.Types.ARRAY || typeName.equalsIgnoreCase("ARRAY");
    public static final TypeChecker DEFAULT_STRING_TYPECHECKER = (type, typeName, precision, scale) ->
            ((type == java.sql.Types.VARCHAR)|| (type == java.sql.Types.CHAR) ||
                    (type == java.sql.Types.LONGVARCHAR) || (type == java.sql.Types.CLOB)
                    || (type == java.sql.Types.NCLOB) ||
                    typeName.equalsIgnoreCase("varchar") ||
                    typeName.equalsIgnoreCase("nvarchar") ||
                    typeName.equalsIgnoreCase("nchar") ||
                    typeName.equalsIgnoreCase("ntext")) &&
                    !typeName.equalsIgnoreCase("uuid") &&
                    !typeName.equalsIgnoreCase("set");
    public static final TypeChecker DEFAULT_XML_TYPECHECKER = (type, typeName, precision, scale) ->
            (type == java.sql.Types.SQLXML || typeName.equalsIgnoreCase("xml"));
}
