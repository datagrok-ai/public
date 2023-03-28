package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.*;

import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class Db2ObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(5)
                .setColumn(new StringColumn(new String[] {"DATAGROK", "SYSCAT", "SYSIBM", "SYSIBMADM",
                        "SYSSTAT"}), "TABLE_SCHEMA")
                .build();
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
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"BOOL", "COUNTRY", "DATE", "EMAIL",
                        "FIRST_NAME", "GENDER", "ID", "IP_ADDRESS", "LAST_NAME", "SOME_NUMBER"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"BOOLEAN", "VARCHAR",
                        "DATE", "VARCHAR", "VARCHAR", "VARCHAR",
                        "BIGINT", "VARCHAR", "VARCHAR", "DECIMAL"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
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

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{32767, 1}), "SMALLINT_TYPE")
                .setColumn(new IntColumn(new Integer[]{2147483647, -2147483648}), "INTEGER_TYPE")
                .setColumn(new BigIntColumn(new String[]{"9223372036854775807", "12415"}), "BIGINT_TYPE")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.fromQuery("SELECT * FROM INTEGER_TYPES;");
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[]{0.1243124f, 0.0004124f}), "DECIMAL_TYPE")
                .setColumn(new FloatColumn(new Float[]{0.0f, -0.0f}), "DECFLOAT_TYPE")
                .setColumn(new FloatColumn(new Float[]{214412.244f, -5.4E-9f}), "REAL_TYPE")
                .setColumn(new FloatColumn(new Float[]{-0.0f, 0.0f}), "DOUBLE_TYPE")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("INTEGER TYPES SUPPORT", funcCall1), expected1),
                Arguments.of(Named.of("FLOAT TYPES SUPPORT", funcCall2), expected2)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_xmlTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[]{"<foo>Hello World!</foo>",
                        "<book><title>Manual</title><chapter>...</chapter></book>"}), "XML_TYPE")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM XML_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("XML TYPE SUPPORT", funcCall), expected)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{
                        parser.parseDateToDouble("yyyy-MM-dd", "1999-01-08")}),
                        "DATE_TYPE")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("HH:mm:ss", "04:05:06")}),
                        "TIME_TYPE")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                                "1999-01-08 04:05:06")}),
                        "TIMESTAMP_TYPE")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_characterTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "CHARACTER_TYPE")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "VARCHAR_TYPE")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "CLOB_TYPE")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "GRAPHIC_TYPE")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "VARGRAPHIC_TYPE")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "DBCLOB_TYPE")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM CHARACTER_TYPES;");
        return Stream.of(
                Arguments.of(Named.of("CHARACTER TYPES SUPPORT", funcCall), expected)
        );
    }
}
