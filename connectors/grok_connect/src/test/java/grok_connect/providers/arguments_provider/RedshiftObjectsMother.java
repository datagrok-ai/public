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
import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class RedshiftObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "public";
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
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"bigint", "character varying",
                        "character varying", "character varying", "character varying", "character varying",
                        "boolean", "character varying", "date", "numeric"})
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

    public static Stream<Arguments> checkDatesParameterSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        LocalDate now = LocalDate.now();
        int dayOfWeek = now.getDayOfWeek().getValue();
        int dayOfMonth = now.getDayOfMonth();
        int dayOfYear = now.getDayOfYear();
        LocalDate firstDayOfWeek = now.minusDays(dayOfWeek - 1);
        LocalDate lastDayOfWeek = now.plusDays(7 - dayOfWeek);
        LocalDate firstDayOfMonth = now.with(TemporalAdjusters.firstDayOfMonth());
        LocalDate lastDayOfMonth = now.with(TemporalAdjusters.lastDayOfMonth());
        LocalDate firstDayOfYear = now.with(TemporalAdjusters.firstDayOfYear());
        LocalDate lastDayOfYear = now.with(TemporalAdjusters.lastDayOfYear());
        LocalDate yesterday = now.minusDays(1);
        LocalDate dayOfLastYear = now.minusDays(150);
        // --input: string date = "today" {pattern: datetime}
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble(datePattern, now.toString())}),
                        "date")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"today\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date)\n"
                        + "-- end")
                .addFuncParam("string", "","date", "today", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfWeek == 1 ? 2 : 3)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString())),
                        "date")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"this week\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string", "","date", "this week", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfMonth > 1 && dayOfMonth < 31 - 6 ? 3 : 2)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())),
                        "date")
                .build();
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"this month\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string", "","date", "this month", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfYear > 1 && dayOfYear < Year.of(now.getYear()).length() - 6 ? 3 : 2)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(),
                                dayOfLastYear.getYear() == now.getYear() ? dayOfLastYear.toString() : null)),
                        "date")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"this year\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string","", "date", "this year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, yesterday.toString())),
                        "date")
                .build();
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date)\n"
                        + "-- end")
                .addFuncParam("string","", "date", "yesterday", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                yesterday.getYear() < now.getYear() ? yesterday.toString() : null,
                                dayOfLastYear.getYear() < now.getYear() ? dayOfLastYear.toString() : null)),
                        "date")
                .build();
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"last year\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string","","date", "last year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(5)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2021-04-09",
                                dayOfLastYear.toString(), yesterday.toString(),
                                now.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())),
                        "date")
                .build();
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery( "-- input: string date = \"anytime\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string", "","date", "anytime", "datetime")
                .addFuncCallOptionsPattern("date", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2021-04-09")),
                        "date")
                .build();

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "-- input: string date = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date)\n"
                        + "-- end")
                .addFuncParam("string", "","date", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string", "","date", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, dayOfLastYear.toString(), yesterday.toString(),
                                now.toString(), lastDayOfWeek.getMonthValue() >
                                        lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())),
                        "date")
                .build();
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date;\n"
                        + "-- end")
                .addFuncParam("string", "","date", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT * FROM dates_patterns WHERE @date(date) ORDER BY date\n"
                        + "-- end")
                .addFuncParam("string", "","date", "April 2021", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        Year.of(2021).atMonth(4).atDay(1).toString(),
                        Year.of(2021).atMonth(5).atDay(1).toString())
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: string; operator: today; pattern: datetime", funcCall1), expected1),
                Arguments.of(Named.of("type: string; operator: this week; pattern: datetime", funcCall2), expected2),
                Arguments.of(Named.of("type: string; operator: this month; pattern: datetime", funcCall3), expected3),
                Arguments.of(Named.of("type: string; operator: this year; pattern: datetime", funcCall4), expected4),
                Arguments.of(Named.of("type: string; operator: yesterday; pattern: datetime", funcCall5), expected5),
                Arguments.of(Named.of("type: string; operator: last year; pattern: datetime", funcCall6), expected6),
                Arguments.of(Named.of("type: string; operator: anytime; pattern: datetime", funcCall7), expected7),
                Arguments.of(Named.of("type: string; operator: range -; pattern: datetime", funcCall8), expected8),
                Arguments.of(Named.of("type: string; operator: before; pattern: datetime", funcCall9), expected8),
                Arguments.of(Named.of("type: string; operator: after; pattern: datetime", funcCall10), expected9),
                Arguments.of(Named.of("type: string; operator: April 2021; pattern: datetime", funcCall11), expected8)
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_binaryTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"Grok"}), "varbyte_type")
                .setColumn(new StringColumn(new String[] {"Grok Connect"}), "varbinary_type")
                .build();
        return Stream.of(Arguments.of(Named.of("BINARY TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT from_varbyte(varbyte_type, 'utf-8') varbyte_type, "
                        + "from_varbyte(varbinary_type, 'utf-8') varbinary_type FROM BINARY_TYPES")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_characterTypes_ok() {
        String bpchar = new String(new char[249]).replace("\0", " ");
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"Hello"}), "char_type")
                .setColumn(new StringColumn(new String[] {"Datagrok"}), "varchar_type")
                .setColumn(new StringColumn(new String[] {"Groks"}), "nchar_type")
                .setColumn(new StringColumn(new String[] {"Datagrok"}), "nvarchar_type")
                .setColumn(new StringColumn(new String[] {"Hello World!"}), "text_type")
                .setColumn(new StringColumn(new String[] {"Groking" + bpchar}), "bpchar_type")
                .build();
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM CHARACTER_TYPES")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        Parser parser = new DateParser();
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "2008-06-01")}), "date_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("hh:mm:ss",
                        "00:00:00")}), "time_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("hh:mm:ss.SSS X",
                        "04:05:06.789 +00:00")}), "timetz_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss",
                        "2008-06-01 09:59:59"),}), "timestamp_type")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss X",
                        "2001-02-16 18:38:40 +00:00"),}), "timestamptz_type")
                .build();
        return Stream.of(Arguments.of(Named.of("DATE, TIME, TIMETZ, TIMESTAMP, TIMESTAMPTZ TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM DATE_TYPES")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{0, -32768}), "smallint_type")
                .setColumn(new IntColumn(new Integer[]{2147483647, -2147483648}), "int_type")
                .setColumn(new BigIntColumn(new String[]{"9223372036854775807", "-9223372036854775808"}),
                        "bigint_type")
                .build();
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new FloatColumn(new Float[]{5.23553f, 0.9999f}), "decimal_type")
                .setColumn(new FloatColumn(new Float[]{-2E-10f, 214144412.2f}), "real_type")
                .setColumn(new FloatColumn(new Float[]{3.14f, Float.NEGATIVE_INFINITY}), "double_precision_type")
                .build();
        return Stream.of(Arguments.of(Named.of("INTEGER TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM INTEGER_TYPES")), expected1),
                Arguments.of(Named.of("FLOAT TYPES SUPPORT",
                        FuncCallBuilder.fromQuery("SELECT * FROM FLOAT_TYPES")), expected2));
    }

    public static Stream<Arguments> checkOutputDataFrame_superType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[] {"{\"bar\":\"baz\",\"balance\":7.77,\"active\":false}",
                "[10,100,1000]"}), "super_type")
                .build();
        return Stream.of(Arguments.of(Named.of("SUPER TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM SUPER_TYPE")), expected));
    }
}
