package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
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

@SuppressWarnings("unused")
public class ClickHouseObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("table_schema", new String[] {"INFORMATION_SCHEMA", "default",
                                "information_schema", "system"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "default";
        String table = "mock_data";
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn(firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema}),
                new StringColumn(secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table}),
                new StringColumn(thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "date", "some_number"}),
                new StringColumn(fourthColumnName, new String[] {"UInt64", "String",
                        "String", "String", "String", "String",
                        "Bool", "String", "Date", "Float64"}),
                new IntColumn(fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkListParameterSupport_ok() {
        // list<string>
        List<String> values = new ArrayList<>();
        values.add("Poland");
        values.add("Brazil");
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"2", "5", "20"}),
                new StringColumn("first_name", new String[]{"Nicholle", "Mitchell", "Lucius",}),
                new StringColumn("last_name", new String[]{"Karoly", "Haglington", "Edelmann"}),
                new StringColumn("email", new String[]{"nkaroly1@alexa.com", "mhaglington4@indiegogo.com",
                        "ledelmannj@bravesites.com"}),
                new StringColumn("gender", new String[]{"Female", "Male", "Male"}),
                new StringColumn("ip_address", new String[]{"255.233.247.118/32", "209.93.181.190/32", "66.174.30.225/32"}),
                new BoolColumn("bool", new Boolean[]{false, true, false}),
                new StringColumn("country", new String[]{"Poland", "Poland", "Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles("yyyy-MM-dd", "2014-02-27",
                                "2020-10-09","1999-06-22")),
                new FloatColumn("some_number", new Float[]{864.09f, 15.22f, 378.73f}));
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

    public static Stream<Arguments> checkOutputDataFrame_uuidType_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("uuid", new String[] {"417ddc5d-e556-4d27-95dd-a34d84e46a50"}));
        return Stream.of(Arguments.of(Named.of("UUID TYPE",
                        FuncCallBuilder.fromQuery("SELECT * FROM uuid_type;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericType_ok() {
        DataFrame expected1 = DataFrame.fromColumns(
                new IntColumn("int8_type", new Integer[]{127, -128}),
                new IntColumn("int16_type", new Integer[]{32767, -32768}),
                new IntColumn("int32_type", new Integer[]{2147483647, -2147483648}),
                new BigIntColumn("int64_type", new String[] {"9223372036854775807", "-9223372036854775808"}),
                new BigIntColumn("int128_type", new String[] {"170141183460469231731687303715884105727",
                        "-170141183460469231731687303715884105728"}),
                new BigIntColumn("int256_type", new String[] {"57896044618658097711785492504343953926634992332820282019728792003956564819967",
                        "-57896044618658097711785492504343953926634992332820282019728792003956564819968"}));
        DataFrame expected2 = DataFrame.fromColumns(
                new IntColumn("uint8_type", new Integer[]{255}),
                new IntColumn("uint16_type", new Integer[]{65535}),
                new BigIntColumn("uint32_type", new String[]{"4294967295"}),
                new BigIntColumn("uint64_type", new String[] {"18446744073709551615"}),
                new BigIntColumn("uint128_type", new String[] {"340282366920938463463374607431768211455"}),
                new BigIntColumn("uint256_type", new String[] {"11579208923731619542357098500868790785"
                        + "3269984665640564039457584007913129639935"}));
        DataFrame expected3 = DataFrame.fromColumns(
                new FloatColumn("float32_type", new Float[]{Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY}),
                new FloatColumn("float64_type", new Float[]{Float.POSITIVE_INFINITY, Float.NaN}),
                new FloatColumn("decimal32_type", new Float[]{999.9999f, -999.9999f}),
                new FloatColumn("decimal64_type", new Float[]{1.0E9f, -1.0E9f}),
                new FloatColumn("decimal128_type", new Float[]{2.11111117E11f, -2.11111117E10f}),
                new FloatColumn("decimal256_type", new Float[]{9.9999998E17f, -1.0E21f}));
        return Stream.of(Arguments.of(Named.of("SIGNED INTEGER TYPES",
                FuncCallBuilder.fromQuery("SELECT * FROM SIGNED_INTEGER_TYPES ORDER BY int8_type DESC;")), expected1),
                Arguments.of(Named.of("UNSIGNED INTEGER TYPES",
                        FuncCallBuilder.fromQuery("SELECT * FROM UNSIGNED_INTEGER_TYPES;")), expected2),
                Arguments.of(Named.of("FLOAT TYPES",
                        FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES ORDER BY float32_type DESC;")), expected3));
    }

    public static Stream<Arguments> checkOutputDataFrame_arrayType_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("array_type", new String[] {"[1, 2, 3, 4]"}),
                new StringColumn("array_array_type", new String[] {"[[24], [421, 12, 4], []]"}));
        return Stream.of(Arguments.of(Named.of("ARRAY TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM ARRAY_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_tupleType_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("tuple", new String[] {"[1244, Datagrok]"}));
        return Stream.of(Arguments.of(Named.of("TUPLE TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM TUPLE_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_mapType_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("a", new String[] {"{key1=1, key2=10}",
                        "{key1=2, key2=20}", "{key1=3, key2=30}"}));
        return Stream.of(Arguments.of(Named.of("MAP TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM MAP_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new DateTimeColumn("date_type", new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd", "2016-06-15")}),
                new DateTimeColumn("date32_type", new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd", "2016-06-15")}),
                new DateTimeColumn("datetime_type", new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                                        "2016-06-15 23:00:00")}),
                new DateTimeColumn("datetime64_type", new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd'T'HH:mm:ssX",
                                        "2019-01-01T02:00:00+0200")}));
        return Stream.of(Arguments.of(Named.of("DATE TYPES",
                FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_nestedType_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new DateTimeColumn("EventDate", new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                                "2016-01-01")}),
                new BigIntColumn("UserID", new String[]{"123"}),
                new StringColumn("Attrs.Key", new String[]{"[price, color]"}),
                new StringColumn("Attrs.Value", new String[]{"[high, red]"}));
        return Stream.of(Arguments.of(Named.of("NESTED TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM NESTED_TYPE;")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_geoTypes_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("p", new String[]{"[10.0, 10.0]"}),
                new StringColumn("r", new String[]{"[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]"}),
                new StringColumn("pg", new String[]{"[[[20.0, 20.0], [50.0, 20.0], [50.0, 50.0], [20.0, 50.0]], "
                        + "[[30.0, 30.0], [50.0, 50.0], [50.0, 30.0]]]"}),
                new StringColumn("mpg", new String[]{"[[[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]], "
                                + "[[[20.0, 20.0], [50.0, 20.0], [50.0, 50.0], [20.0, 50.0]], [[30.0, 30.0], [50.0, 50.0], [50.0, 30.0]]]]"}));
        return Stream.of(Arguments.of(Named.of("GEO TYPE",
                FuncCallBuilder.fromQuery("SELECT * FROM GEO_TYPES;")), expected));
    }
}
