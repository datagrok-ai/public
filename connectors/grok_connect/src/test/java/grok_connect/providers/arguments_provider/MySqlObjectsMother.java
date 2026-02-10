package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BigIntColumn;
import serialization.BoolColumn;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class MySqlObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"datagrok", "information_schema", "performance_schema"}),
                        "table_schema")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "datagrok";
        String table = "mock_data";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "date", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"bigint", "varchar",
                        "varchar", "varchar", "varchar", "varchar",
                        "tinyint", "varchar", "date", "decimal"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkListParameterSupport_ok() {
        // list<string>
        List<String> values = new ArrayList<>();
        values.add("Poland");
        values.add("Brazil");
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new BigIntColumn(new String[]{"2", "5", "20"}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Nicholle", "Mitchell", "Lucius",}), "first_name")
                .setColumn(new StringColumn(new String[]{"Karoly", "Haglington", "Edelmann"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"nkaroly1@alexa.com", "mhaglington4@indiegogo.com",
                        "ledelmannj@bravesites.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Female", "Male", "Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"255.233.247.118/32", "209.93.181.190/32", "66.174.30.225/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false, true, false}), "bool")
                .setColumn(new StringColumn(new String[]{"Poland", "Poland", "Brazil"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd", "2014-02-27",
                                "2020-10-09","1999-06-22")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{864.09f, 15.22f, 378.73f}), "some_number")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values\n" +
                        "SELECT * FROM mock_data WHERE country IN (@values)")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values = ['Poland','Brazil']\n" +
                        "SELECT * FROM mock_data WHERE country IN (@values)")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall1), expected),
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall2), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_jsonType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"{\"key1\": \"value1\", \"key2\": \"value2\"}",
                                "{\"phones\": [{\"type\": \"mobile\", \"phone\": \"001001\"},"
                                        + " {\"type\": \"fix\", \"phone\": \"002002\"}]}", "{\"reading\": 0.0000123}"}),
                        "json_type")
                .build();
        return Stream.of(Arguments.of(Named.of("JSON TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM JSON_TYPE")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_bitType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[]{"10000000", "1"}), "bit_type1")
                .setColumn(new BoolColumn(new Boolean[]{true, false}), "bit_type2")
                .build();
        return Stream.of(Arguments.of(Named.of("BIT TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT BIN(bit_type1) AS bit_type1, bit_type2 FROM BIT_TYPE")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "9999-12-31")}), "date_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("hh:mm:ss.SSS" ,
                        "18:59:59.000000")}), "time_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS",
                        "1970-01-01 00:00:01.000000")}), "timestamp_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS",
                                "1000-01-01 00:00:00.000000")}),
                        "datetime_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "2023-01-01")}), "year_type")
                .build();
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_spatialTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[]{"POLYGON((0 0,10 0,10 10,0 10,0 0),(5 5,7 5,7 7,5 7,5 5))",
                "GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 1,2 2,3 3,4 4))"}), "geometry_type")
                .setColumn(new StringColumn(new String[]{"POINT(1 1)", "POINT(1 0)"}), "point_type")
                .build();
        return Stream.of(Arguments.of(Named.of("SPATIAL TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT ST_AsText(geometry_type) AS geometry_type, " +
                        "ST_AsText(point_type) AS point_type FROM GEOMETRY_TYPE")),
                expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_integerTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{12, 0}), "tinyint_type")
                .setColumn(new IntColumn(new Integer[]{32000, 1212}), "smallint_type")
                .setColumn(new IntColumn(new Integer[]{167772, -1000}), "mediumint_type")
                .setColumn(new IntColumn(new Integer[]{2147483647, -2147483648}), "int_type")
                .setColumn(new BigIntColumn(new String[]{"9223372036854775807", "-9223372036854775808"}),
                        "bigint_type")
                .build();
        return Stream.of(Arguments.of(Named.of("INTEGER TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM INTEGER_TYPES")),
                expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_floatTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new FloatColumn(new Float[]{-3.40282E+38f, 3.40282E+38f, -4.54654f,
                        -0.000124f}), "float_type")
                .setColumn(new FloatColumn(new Float[]{Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY,
                        4.457745745745457f, 0.002f}), "double_type")
                .setColumn(new FloatColumn(new Float[]{999.9999f, 0.9999f, 23542363246234234234.46456456f, 0.00001f}),
                        "decimal_type")
                .build();
        return Stream.of(Arguments.of(Named.of("FLOAT TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES")),
                expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_characterTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"datagrok"}), "char_type1")
                .setColumn(new StringColumn(new String[]{"A"}), "char_type2")
                .setColumn(new StringColumn(new String[]{"Hello World"}), "varchar_type")
                .setColumn(new StringColumn(new String[]{"Hello Datagrok"}), "text_type")
                .setColumn(new StringColumn(new String[]{"HelloolleH"}), "mediumtext_type")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "longtext_type")
                .build();
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM CHARACTER_TYPES")),
                expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_binaryTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"Hello"}), "binary_type")
                .setColumn(new StringColumn(new String[]{"datagrok"}), "varbinary_type")
                .build();
        return Stream.of(Arguments.of(Named.of("BINARY TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT CAST(binary_type as CHAR(5)) as binary_type, CAST(varbinary_type AS CHAR(8)) "
                                + "as varbinary_type FROM BINARY_TYPES")),
                expected));
    }
}
