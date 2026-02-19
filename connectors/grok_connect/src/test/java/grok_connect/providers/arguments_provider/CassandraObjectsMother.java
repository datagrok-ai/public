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
import java.time.LocalDate;
import java.time.Year;
import java.time.temporal.TemporalAdjusters;
import java.util.stream.Stream;

@SuppressWarnings("unused")
public class CassandraObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn("TABLE_SCHEMA", new String[] {"datagrok"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String schema = "datagrok";
        String table = "mock_data";
        DataFrame expected = DataFrame.fromColumns(
                new StringColumn(firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema}),
                new StringColumn(secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table}),
                new StringColumn(thirdColumnName, new String[] {"bool", "country", "date",
                        "email", "first_name", "gender", "id", "ip_address", "last_name", "some_number"}),
                new StringColumn(fourthColumnName, new String[] {"boolean", "text", "date",
                        "text", "text", "text", "int", "text", "text", "double"}));
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkParameterSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        // --input: int id = 20
        DataFrame expected1 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{20}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "1999-06-22")),
                new StringColumn("email", new String[]{"ledelmannj@bravesites.com"}),
                new StringColumn("first_name", new String[]{"Lucius"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"66.174.30.225/32"}),
                new StringColumn("last_name", new String[]{"Edelmann"}),
                new FloatColumn("some_number", new Float[]{378.73f}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: int id = 20\n"
                        + "SELECT * FROM datagrok.mock_data WHERE id = @id ALLOW FILTERING")
                .addFuncParam("int", "", "id", 20, "")
                .build();
        // --input: string id = ">28" {pattern: int}
        DataFrame expected2 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{30, 29}),
                new BoolColumn("bool", new Boolean[]{false, false}),
                new StringColumn("country", new String[]{"France", "Sweden"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                        "2016-07-10", "2009-10-02")),
                new StringColumn("email", new String[]{"blonglandst@tripod.com", "gfayters@desdev.cn"}),
                new StringColumn("first_name", new String[]{"Bran", "Grantham"}),
                new StringColumn("gender", new String[]{"Genderqueer", "Male"}),
                new StringColumn("ip_address", new String[]{"14.92.3.30/32", "26.120.76.78/32"}),
                new StringColumn("last_name", new String[]{"Longlands", "Fayter"}),
                new FloatColumn("some_number", new Float[]{879.94f, 595.22f}));

        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \">28\" {pattern: int}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @id(id) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "id", ">28", "int")
                .addFuncCallOptionsPattern("id", ">28", ">",
                        null, null, 28)
                .build();
        // input: string id = ">=29" {pattern: int}
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \">=29\" {pattern: int}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @id(id) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "id", ">=29", "int")
                .addFuncCallOptionsPattern("id", ">=29", ">=",
                        null, null, 29)
                .build();
        // --input: string id = "<=1" {pattern: int}
        DataFrame expected3 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{1}),
                new BoolColumn("bool", new Boolean[]{true}),
                new StringColumn("country", new String[]{"China"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2017-09-20")),
                new StringColumn("email", new String[]{"bkemery0@businesswire.com"}),
                new StringColumn("first_name", new String[]{"Burk"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"249.64.22.121/32"}),
                new StringColumn("last_name", new String[]{"Kemery"}),
                new FloatColumn("some_number", new Float[]{510.32f}));
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"<=1\" {pattern: int}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @id(id) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "id", "<=1", "int")
                .addFuncCallOptionsPattern("id", "<=1", "<=",
                        null, null, 1)
                .build();
        // --input: string id = "<2" {pattern: int}
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"<2\" {pattern: int}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @id(id) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "id", "<2", "int")
                .addFuncCallOptionsPattern("id", "<2", "<",
                        null, null, 2)
                .build();
        // --input: string id = "in(29, 30)" {pattern: int}
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"in(29, 30)\" {pattern: int}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "","id", "in(29, 30)", "int")
                .addFuncCallOptionsPattern("id", "in(29, 30)", "in",
                        null, null, 29, 30)
                .build();
        // --input: string id = "min-max 29-30" {pattern: int}
        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"min-max 29-30\" {pattern: int}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @id(id) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "id", "min-max 29-30",
                        "int")
                .addFuncCallOptionsPattern("id", "29-30",
                        "-", null, null, 29, 30)
                .build();
        //--input: double some_number = 510.32
        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: double some_number = 510.32\n"
                        + "SELECT * FROM datagrok.mock_data WHERE some_number = @some_number ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("double", "", "some_number", 510.32, "double")
                .build();
        // --input: string some_number = ">975" {pattern: double}
        DataFrame expected5 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{10, 26}),
                new BoolColumn("bool", new Boolean[]{false, false}),
                new StringColumn("country", new String[]{"Vietnam", "Honduras"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2003-01-04",
                        "2010-05-04")),
                new StringColumn("email", new String[]{"sformilli9@aol.com", "doshaughnessyp@com.com"}),
                new StringColumn("first_name", new String[]{"Scottie", "Daryle"}),
                new StringColumn("gender", new String[]{"Male", "Male"}),
                new StringColumn("ip_address", new String[]{"101.241.191.228/32", "204.107.16.207/32"}),
                new StringColumn("last_name", new String[]{"Formilli", "O'Shaughnessy"}),
                new FloatColumn("some_number", new Float[]{978.01f, 983.03f}));
        FuncCall funcCall10 =  FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \">975.0\" {pattern: double}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @some_number(some_number) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string",  "", "some_number", ">975", "double")
                .addFuncCallOptionsPattern("some_number", ">975", ">",
                        null, null, 975.0)
                .build();
        // --input: string some_number = ">=975" {pattern: double}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \">=975.0\" {pattern: double}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @some_number(some_number) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "some_number", ">=975", "double")
                .addFuncCallOptionsPattern("some_number", ">=975", ">=",
                        null, null, 975.0)
                .build();
        //--input: string some_number = "<20" {pattern: double}
        DataFrame expected6 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{5}),
                new BoolColumn("bool", new Boolean[]{true}),
                new StringColumn("country", new String[]{"Poland"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2020-10-09")),
                new StringColumn("email", new String[]{"mhaglington4@indiegogo.com"}),
                new StringColumn("first_name", new String[]{"Mitchell"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"209.93.181.190/32"}),
                new StringColumn("last_name", new String[]{"Haglington"}),
                new FloatColumn("some_number", new Float[]{15.22f}));
        FuncCall funcCall12 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \"<20.0\" {pattern: double}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @some_number(some_number) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "some_number", "<20", "double")
                .addFuncCallOptionsPattern("some_number", "<20", "<",
                        null, null, 20.0)
                .build();
        // --input: string some_number = "<=20" {pattern: double}
        FuncCall funcCall13 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \"<=20.0\" {pattern: double}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @some_number(some_number) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "", "some_number", "<=20", "double")
                .addFuncCallOptionsPattern("some_number", "<=20", "<=",
                        null, null, 20.0)
                .build();
        // --input: string first_name = 'contains Z' {pattern: string}
        DataFrame expected7 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{25}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Bosnia and Herzegovina"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2003-02-12")),
                new StringColumn("email", new String[]{"zwimmerso@hatena.ne.jp"}),
                new StringColumn("first_name", new String[]{"Zolly"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"123.12.225.114/32"}),
                new StringColumn("last_name", new String[]{"Wimmers"}),
                new FloatColumn("some_number", new Float[]{217.18f}));
        FuncCall funcCall14 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = 'contains Z' {pattern: string}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @first_name(first_name) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "contains Z", "string")
                .addFuncCallOptionsPattern("first_name", "contains Z", "contains",
                        null, null, "Z")
                .build();
        // --input: string first_name = 'starts with W' {pattern: string}
        DataFrame expected8 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{23}),
                new BoolColumn("bool", new Boolean[]{true}),
                new StringColumn("country", new String[]{"Sweden"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2011-12-18")),
                new StringColumn("email", new String[]{"wroglierom@berkeley.edu"}),
                new StringColumn("first_name", new String[]{"Waly"}),
                new StringColumn("gender", new String[]{"Female"}),
                new StringColumn("ip_address", new String[]{"122.90.196.231/32"}),
                new StringColumn("last_name", new String[]{"Rogliero"}),
                new FloatColumn("some_number", new Float[]{147.69f}));
        FuncCall funcCall15 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = 'starts with W' {pattern: string}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @first_name(first_name) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "starts with W", "string")
                .addFuncCallOptionsPattern("first_name", "starts with W", "starts with",
                        null, null, "W")
                .build();
        // --input: string first_name = 'ends with y' {pattern: string}
        DataFrame expected9 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{20}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "1999-06-22")),
                new StringColumn("email", new String[]{"ledelmannj@bravesites.com"}),
                new StringColumn("first_name", new String[]{"Lucius"}),
                new StringColumn("gender", new String[]{"Male"}),
                new StringColumn("ip_address", new String[]{"66.174.30.225/32"}),
                new StringColumn("last_name", new String[]{"Edelmann"}),
                new FloatColumn("some_number", new Float[]{378.73f}));
        FuncCall funcCall16 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = 'ends with s' {pattern: string}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @first_name(first_name) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "ends with s", "string")
                .addFuncCallOptionsPattern("first_name", "ends with s", "ends with",
                        null, null, "s")
                .build();
        // --input: string country = 'in (Poland, Brazil)' {pattern: string}
        DataFrame expected10 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{5, 2, 20}),
                new BoolColumn("bool", new Boolean[]{true, false, false}),
                new StringColumn("country", new String[]{"Poland", "Poland", "Brazil"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2020-10-09", "2014-02-27", "1999-06-22")),
                new StringColumn("email", new String[]{"mhaglington4@indiegogo.com", "nkaroly1@alexa.com", "ledelmannj@bravesites.com"}),
                new StringColumn("first_name", new String[]{"Mitchell", "Nicholle", "Lucius",}),
                new StringColumn("gender", new String[]{"Male", "Female", "Male"}),
                new StringColumn("ip_address", new String[]{"209.93.181.190/32", "255.233.247.118/32", "66.174.30.225/32"}),
                new StringColumn("last_name", new String[]{"Haglington", "Karoly", "Edelmann"}),
                new FloatColumn("some_number", new Float[]{15.22f, 864.09f, 378.73f}));
        FuncCall funcCall17 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string country = 'in (Poland, Brazil)' {pattern: string}\n" +
                        "SELECT * FROM datagrok.mock_data WHERE @country(country) ALLOW FILTERING\n" +
                        "--end")
                .addFuncParam("string", "", "country", "in (Poland, Brazil)", "string")
                .addFuncCallOptionsPattern("country", "in (Poland, Brazil)", "in",
                        null, null, "Poland", "Brazil")
                .build();
        DataFrame expected11 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{29, 30}),
                new BoolColumn("bool", new Boolean[]{false, false}),
                new StringColumn("country", new String[]{"Sweden", "France"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                        "2009-10-02", "2016-07-10")),
                new StringColumn("email", new String[]{"gfayters@desdev.cn", "blonglandst@tripod.com"}),
                new StringColumn("first_name", new String[]{"Grantham", "Bran"}),
                new StringColumn("gender", new String[]{"Male", "Genderqueer"}),
                new StringColumn("ip_address", new String[]{"26.120.76.78/32", "14.92.3.30/32"}),
                new StringColumn("last_name", new String[]{"Fayter", "Longlands"}),
                new FloatColumn("some_number", new Float[]{595.22f, 879.94f}));
        return Stream.of(
                Arguments.of(Named.of("type: int; operator: =; pattern: none", funcCall1), expected1),
                Arguments.of(Named.of("type: string; operator: >; pattern: int", funcCall2), expected2),
                Arguments.of(Named.of("type: string; operator: >=; pattern: int", funcCall3), expected2),
                Arguments.of(Named.of("type: string; operator: <=; pattern: int", funcCall4), expected3),
                Arguments.of(Named.of("type: string; operator: <; pattern: int", funcCall5), expected3),
                Arguments.of(Named.of("type: string; operator: in; pattern: int", funcCall6), expected11),
                Arguments.of(Named.of("type: string; operator: min-max; pattern: int", funcCall8), expected2),
                Arguments.of(Named.of("type: double; operator: =; pattern: none", funcCall9), expected3),
                Arguments.of(Named.of("type: string; operator: >; pattern: double", funcCall10), expected5),
                Arguments.of(Named.of("type: string; operator: >=; pattern: double", funcCall11), expected5),
                Arguments.of(Named.of("type: string; operator: <; pattern: double", funcCall12), expected6),
                Arguments.of(Named.of("type: string; operator: <=; pattern: double", funcCall13), expected6),
                Arguments.of(Named.of("type: string; operator: contains; pattern: string", funcCall14), expected7),
                Arguments.of(Named.of("type: string; operator: starts with; pattern: string", funcCall15), expected8),
                Arguments.of(Named.of("type: string; operator: ends with; pattern: string", funcCall16), expected9),
                Arguments.of(Named.of("type: string; operator: in; pattern: string", funcCall17), expected10)
        );
    }

    public static Stream<Arguments> checkMultipleParametersSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        // --input: string first_name = "starts with p" {pattern: string}
        //--input: string id = ">1" {pattern :int}
        //--input: bool bool = false
        //--input: string email = "contains com" {pattern: string}
        //--input: string some_number = ">20" {pattern: double}
        //--input: string country = "in (Indonesia)" {pattern: string}
        //--input: string date = "before 1/1/2022" {pattern: datetime}
        DataFrame expected1 = DataFrame.fromColumns(
                new IntColumn("id", new Integer[]{13}),
                new BoolColumn("bool", new Boolean[]{false}),
                new StringColumn("country", new String[]{"Indonesia"}),
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2012-01-14")),
                new StringColumn("email", new String[]{"pboxellc@moonfruit.com"}),
                new StringColumn("first_name", new String[]{"Pail"}),
                new StringColumn("gender", new String[]{"Genderqueer"}),
                new StringColumn("ip_address", new String[]{"2.37.160.155/32"}),
                new StringColumn("last_name", new String[]{"Boxell"}),
                new FloatColumn("some_number", new Float[]{73.47f}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = \"starts with p\" {pattern: string}\n"
                        + "--input: string id = \">1\" {pattern :int}\n"
                        + "--input: bool bool = false\n"
                        + "--input: string email = \"contains com\" {pattern: string}\n"
                        + "--input: string some_number = \">20.0\" {pattern: double}\n"
                        + "--input: string country = \"in (Indonesia)\" {pattern: string}\n"
                        + "--input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM datagrok.mock_data WHERE @first_name(first_name) AND @id(id) AND bool = @bool "
                        + "AND @email(email) AND @some_number(some_number) "
                        + "AND @country(country) AND @date(date) ALLOW FILTERING\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "starts with p", "string")
                .addFuncParam("string", "","id", ">1", "int")
                .addFuncParam("bool", "","bool", false, "")
                .addFuncParam("string", "","email", "contains com", "string")
                .addFuncParam("string", "","some_number", ">20", "double")
                .addFuncParam("string", "","country", "in (Indonesia)", "string")
                .addFuncParam("string", "","date", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("first_name", "starts with p",
                        "starts with", null, null, "p")
                .addFuncCallOptionsPattern("id", ">1", ">", null,
                        null, "1")
                .addFuncCallOptionsPattern("email", "contains com",
                        "contains", null, null, "com")
                .addFuncCallOptionsPattern("some_number", ">20", ">", null,
                        null, 20.0)
                .addFuncCallOptionsPattern("country", "in (Indonesia)", "in",
                        null, null, "Indonesia")
                .addFuncCallOptionsPattern("date", "before 1/1/2022", "before",
                        true, true, Year.of(2022).atMonth(1).atDay(1).toString())
                .build();
        return Stream.of(Arguments.of(Named.of("type: multiple; operator: multiple; pattern: multiple", funcCall1),
                expected1));
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
        DataFrame expected1 = DataFrame.fromColumns(
                new DateTimeColumn("date", new Double[]{parser.parseDateToDouble(datePattern, now.toString())}));
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"today\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string", "","date", "today", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString(),
                                now.toString())));
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"this week\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string", "","date", "this week", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(), now.toString())));
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"this month\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string", "","date", "this month", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                dayOfYear == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(),
                                dayOfLastYear.getYear() == now.getYear() ? dayOfLastYear.toString() : null,
                                now.toString())));
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"this year\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string","", "date", "this year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, yesterday.toString())));
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string","", "date", "yesterday", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern,
                                yesterday.getYear() < now.getYear() ? yesterday.toString() : null,
                                dayOfLastYear.getYear() < now.getYear() ? dayOfLastYear.toString() : null)));
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"last year\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string","","date", "last year", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, "2021-04-09")));

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "-- input: string date = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string", "","date", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
                        + "-- end")
                .addFuncParam("string", "","date", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrame.fromColumns(
                new DateTimeColumn("date", parser.parseDatesToDoubles(datePattern, yesterday.toString(),
                                lastDayOfWeek.toString(), dayOfLastYear.toString(), now.toString())));
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING;\n"
                        + "-- end")
                .addFuncParam("string", "","date", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("date", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string date = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT date FROM datagrok.dates_patterns WHERE @date(date) ALLOW FILTERING\n"
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
                Arguments.of(Named.of("type: string; operator: range -; pattern: datetime", funcCall8), expected8),
                Arguments.of(Named.of("type: string; operator: before; pattern: datetime", funcCall9), expected8),
                Arguments.of(Named.of("type: string; operator: after; pattern: datetime", funcCall10), expected9),
                Arguments.of(Named.of("type: string; operator: April 2021; pattern: datetime", funcCall11), expected8)
        );
    }

    public static Stream<Arguments> checkIntegerTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new IntColumn("tinyint_type", new Integer[]{-128, 127}),
                new BigIntColumn("bigint_type", new String[]{"-999999999999999999", "999999999999999999"}),
                new IntColumn("int_type", new Integer[]{-2147483648, 2147483647}),
                new IntColumn("smallint_type", new Integer[]{-32768, 32767}),
                new BigIntColumn("varint_type", new String[]{"-9999999999999999999", "999999999999999999"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM datagrok.integers;");
        return Stream.of(Arguments.of(Named.of("INTEGER TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkFloatTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new FloatColumn("decimal_type", new Float[]{-20.01001f, 123003.336f}),
                new FloatColumn("double_type", new Float[]{-3.134576f, 3.134576f}),
                new FloatColumn("float_type", new Float[]{-9990.992f, 1.23E-4f}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM datagrok.float_types;");
        return Stream.of(Arguments.of(Named.of("FLOAT TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkDateTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"2", "1"}),
                new DateTimeColumn("date_type", parser.parseDatesToDoubles("yyyy-MM-dd",
                        "2023-04-27", "2011-02-03")),
                new StringColumn("duration_type", new String[]{"9m51s", "89h4m48s"}),
                new DateTimeColumn("stamp_type", parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss X",
                                "2011-02-03 04:05:08 +03:00", "2011-02-03 04:05:01 +01:00")),
                new DateTimeColumn("time_type", parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS",
                        "1970-01-01 12:50:54.333", "1970-01-01 08:12:54.111")));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM datagrok.date_types;");
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkComplexTypesSupport_ok() {
        DataFrame expected = DataFrame.fromColumns(
                new BigIntColumn("id", new String[]{"1"}),
                new StringColumn("list_type", new String[]{"[hello, world]"}),
                new StringColumn("map_type", new String[]{"{band=Beatles, fruit=Apple}"}),
                new StringColumn("set_type", new String[]{"[3, 15, 16]"}),
                new StringColumn("tuple_type", new String[]{"{3, hours}"}),
                new StringColumn("udt_type", new String[]{"{country_code=380, number=+3800900909}"}));
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM datagrok.complex_types;");
        return Stream.of(Arguments.of(Named.of("COMPLEX TYPES SUPPORT", funcCall), expected));
    }
}
