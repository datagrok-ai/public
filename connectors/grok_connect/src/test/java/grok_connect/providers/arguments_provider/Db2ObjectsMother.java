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
public class Db2ObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("TABLE_SCHEMA", new String[] {"DATAGROK", "SYSCAT", "SYSIBM", "SYSIBMADM",
                        "SYSSTAT"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "TABLE_SCHEMA";
        String secondColumnName = "TABLE_NAME";
        String thirdColumnName = "COLUMN_NAME";
        String fourthColumnName = "DATA_TYPE";
        String fifthColumnName = "IS_VIEW";
        String schema = "DATAGROK";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn(firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema}),
                new StringColumn(secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table}),
                new StringColumn(thirdColumnName, new String[] {"BOOL", "COUNTRY", "DATE", "EMAIL",
                        "FIRST_NAME", "GENDER", "ID", "IP_ADDRESS", "LAST_NAME", "SOME_NUMBER"}),
                new StringColumn(fourthColumnName, new String[] {"BOOLEAN", "VARCHAR",
                        "DATE", "VARCHAR", "VARCHAR", "VARCHAR",
                        "BIGINT", "VARCHAR", "VARCHAR", "DECIMAL"}),
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

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected1 = DataFrame.fromColumns(
                new IntColumn("SMALLINT_TYPE", new Integer[]{32767, 1}),
                new IntColumn("INTEGER_TYPE", new Integer[]{2147483647, -2147483648}),
                new BigIntColumn("BIGINT_TYPE", new String[]{"9223372036854775807", "12415"}));
        FuncCall funcCall1 = FuncCallBuilder.fromQuery("SELECT * FROM INTEGER_TYPES;");
        DataFrame expected2 = DataFrame.fromColumns(
                new FloatColumn("DECIMAL_TYPE", new Float[]{0.1243124f, 0.0004124f}),
                new FloatColumn("DECFLOAT_TYPE", new Float[]{0.0f, -0.0f}),
                new FloatColumn("REAL_TYPE", new Float[]{214412.244f, -5.4E-9f}),
                new FloatColumn("DOUBLE_TYPE", new Float[]{-0.0f, 0.0f}));
        FuncCall funcCall2 = FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("INTEGER TYPES SUPPORT", funcCall1), expected1),
                Arguments.of(Named.of("FLOAT TYPES SUPPORT", funcCall2), expected2)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_xmlTypes_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("XML_TYPE", new String[]{"<foo>Hello World!</foo>",
                        "<book><title>Manual</title><chapter>...</chapter></book>"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM XML_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("XML TYPE SUPPORT", funcCall), expected)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new DateTimeColumn("DATE_TYPE", new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd", "1999-01-08")}),
                new DateTimeColumn("TIME_TYPE", new Double[]{parser.parseDateToDouble("HH:mm:ss", "04:05:06")}),
                new DateTimeColumn("TIMESTAMP_TYPE", new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                                "1999-01-08 04:05:06")}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_characterTypes_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("CHARACTER_TYPE", new String[]{"Datagrok"}),
                new StringColumn("VARCHAR_TYPE", new String[]{"Datagrok"}),
                new StringColumn("CLOB_TYPE", new String[]{"Datagrok"}),
                new StringColumn("GRAPHIC_TYPE", new String[]{"Datagrok"}),
                new StringColumn("VARGRAPHIC_TYPE", new String[]{"Datagrok"}),
                new StringColumn("DBCLOB_TYPE", new String[]{"Datagrok"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM CHARACTER_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("CHARACTER TYPES SUPPORT", funcCall), expected)
        );
    }
}
