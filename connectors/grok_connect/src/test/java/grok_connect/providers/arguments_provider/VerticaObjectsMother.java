package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.*;

import java.util.stream.Stream;

@SuppressWarnings("unused")
public class VerticaObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("table_schema", new String[] {"v_internal", "v_func", "store",
                                "v_catalog", "online_sales", "v_txtindex", "v_monitor", "public"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "public";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn(firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema}),
                new StringColumn(secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table}),
                new StringColumn(thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "date", "some_number"}),
                new StringColumn(fourthColumnName, new String[] {"int", "varchar(50)",
                        "varchar(50)", "varchar(50)", "varchar(50)", "varchar(50)",
                        "boolean", "varchar(50)", "date", "numeric(5,2)"}),
                new IntColumn(fifthColumnName, new Integer[] {0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkCharacterTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("char_type", new String[]{"datagrok"}),
                new StringColumn("varchar_type", new String[]{"Hello World"}),
                new StringColumn("longvarchar_type", new String[]{"Datagrok"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM CHARACTER_TYPES;");
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkDateTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new DateTimeColumn("date", new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "1999-01-08")}),
                new DateTimeColumn("time", new Double[]{parser.parseDateToDouble("HH:mm:ss.SSS",
                        "04:05:06.789")}),
                new DateTimeColumn("stamp", new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                        "1999-01-08 04:05:06")}),
                new DateTimeColumn("zoned_time", new Double[]{parser.parseDateToDouble("HH:mm:ss X",
                        "04:05:06 -08:00")}),
                new DateTimeColumn("zoned_stamp", new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss X",
                        "1999-01-08 04:05:06 -08:00")}),
                new StringColumn("interval1", new String[]{"1-0"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM DATES;");
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkNumericTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("int_type", new String[] {"922337", "-92230937", "0"}),
                new BigIntColumn("big_int", new String[] {"922337203685477588", "-922337203685477580", "0"}),
                new FloatColumn("double_precision", new Float[]{Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY,
                        1.24412412E11f}),
                new FloatColumn("float_type", new Float[]{Float.NaN, 0.00124412f,
                        1.0E-5f}),
                new FloatColumn("decimal_type", new Float[]{123456.78f, 1.22223456E8f,
                        1.0E14f}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM NUMERIC_TYPE;");
        return Stream.of(Arguments.of(Named.of("NUMERIC TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkSpatialTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("geography_type", new String[]{"POLYGON ((1 2, 3 4, 2 3, 1 2))"}),
                new StringColumn("geometry_type", new String[]{"POINT (3.14 -1.34)"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT ST_AsText(geography_type) AS geography_type, "
                + "ST_AsText(geometry_type) AS geometry_type FROM SPATIAL;");
        return Stream.of(Arguments.of(Named.of("SPATIAL TYPE SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkUuidTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("uuid_data", new String[]{"9fb01de0-1d63-4d09-9415-90e0b4e93b9a"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM UUID_TYPE;");
        return Stream.of(Arguments.of(Named.of("UUID TYPE SUPPORT", funcCall), expected));
    }
}
