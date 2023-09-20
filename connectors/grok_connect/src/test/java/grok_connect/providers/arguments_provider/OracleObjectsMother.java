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

/**
 * Provides data for Oracle provider tests
 */
public class OracleObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"DATAGROK", "GSMADMIN_INTERNAL",
                        "LBACSYS"}), "TABLE_SCHEMA")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "TABLE_SCHEMA";
        String secondColumnName = "TABLE_NAME";
        String thirdColumnName = "COLUMN_NAME";
        String fourthColumnName = "DATA_TYPE";
        String fifthColumnname = "IS_VIEW";
        String schema = "DATAGROK";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(9)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"ID", "FIRST_NAME", "LAST_NAME", "EMAIL",
                        "GENDER", "IP_ADDRESS", "COUNTRY", "DAT", "SOME_NUMBER"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"NUMBER(18, 0)", "VARCHAR2",
                        "VARCHAR2", "VARCHAR2", "VARCHAR2", "VARCHAR2",
                        "VARCHAR2", "DATE", "NUMBER(5, 2)"})
                .setColumn(new IntColumn(), fifthColumnname, new Integer[]{0, 0, 0, 0, 0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> checkDatesParameterSupport_ok() {
        Parser parser = new DateParser();
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
                        "dat")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"today\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "today", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        now.toString(), now.plusDays(1).toString())
                .build();
        // --input: string date = "this week" {pattern: datetime}
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfWeek == 1 ? 2 : 3)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfWeek == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.equals(now) ? null : lastDayOfWeek.toString())),
                        "dat")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this week\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string","", "dat", "this week", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfWeek.toString(),
                        lastDayOfWeek.plusDays(1).toString())
                .build();
        // --input: string date = "this month" {pattern: datetime}
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfMonth > 1 && dayOfMonth < 31 - 6 ? 3 : 2)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfMonth == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getMonthValue() >  lastDayOfMonth.getMonthValue() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString())),
                        "dat")
                .build();
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this month\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "","dat", "this month", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfMonth.toString(),
                        lastDayOfMonth.plusDays(1).toString())
                .build();
        // --input: string date = "this year" {pattern: datetime}
        DataFrame expected4 = DataFrameBuilder.getBuilder()
                .setRowCount(dayOfYear > 1 && dayOfYear < Year.of(now.getYear()).length() - 6 ? 3 : 2)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                now.toString(),
                                dayOfYear == 1 ? null : yesterday.toString(),
                                lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(),
                                dayOfLastYear.getYear() == now.getYear() ? dayOfLastYear.toString() : null)),
                        "dat")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"this year\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "","dat", "this year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.toString(),
                        lastDayOfYear.plusDays(1).toString())
                .build();
        // --input: string date = "yesterday" {pattern: datetime}
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, yesterday.toString())),
                        "dat")
                .build();
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string date = \"yesterday\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "","dat", "yesterday", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        yesterday.toString(),
                        now.toString())
                .build();
        // --input: string date = "last year" {pattern: datetime}
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,
                                yesterday.getYear() < now.getYear() ? yesterday.toString() : null,
                                dayOfLastYear.getYear() < now.getYear() ? dayOfLastYear.toString() : null)),
                        "dat")
                .build();
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"last year\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "last year", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
                        firstDayOfYear.minusYears(1).toString(), firstDayOfYear.toString())
                .build();
        // --input: string date = "anytime" {pattern: datetime}
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(5)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.getYear() >  now.getYear() || lastDayOfWeek.equals(now)?
                                        null : lastDayOfWeek.toString(), dayOfLastYear.toString(),
                                "2021-04-09")),
                        "dat")
                .build();
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"anytime\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "anytime", "datetime")
                .addFuncCallOptionsPattern("dat", "", "none", true, true)
                .build();
        // --input: string date = "2021-2022" {pattern: datetime}
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2021-04-09")),
                        "dat")
                .build();

        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string dat = \"2021-2021\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "", "dat", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "before 1/1/2022" {pattern: datetime}

        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "","dat", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "before", true, true,
                        Year.of(2022).atDay(1).toString())
                .build();
        // --input: string date = "after 1/1/2022" {pattern: datetime}
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, now.toString(),
                                yesterday.toString(), lastDayOfWeek.toString(), dayOfLastYear.toString())),
                        "dat")
                .build();
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"after 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "","dat", "after 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("dat", "", "after", true, true,
                        LocalDate.parse("2022-01-01").toString())
                .build();
        // --input: string date = "April 2021" {pattern: datetime}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string dat = \"April 2021\" {pattern: datetime}\n"
                        + "SELECT TO_DATE(TO_CHAR (dat, 'YYYY-MM-DD'), 'YYYY-MM-DD') AS DAT FROM dates_patterns WHERE @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "","dat", "April 2021", "datetime")
                .addFuncCallOptionsPattern("dat", "", "range", true, false,
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

    public static Stream<Arguments> checkParameterSupport_ok() {
        String datePattern = "yyyy-MM-dd";
        // --input: int id = 20
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{20}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Lucius"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Edelmann"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"ledelmannj@bravesites.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"66.174.30.225/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Brazil"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "1999-06-22")), "date")
                .setColumn(new FloatColumn(new Float[]{378.73f}), "some_number")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: int id = 20\n"
                        + "SELECT * FROM mock_data WHERE id = @id")
                .addFuncParam("int", "", "id", 20, "")
                .build();
        // --input: string id = ">28" {pattern: int}
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{29, 30}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Grantham", "Bran"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Fayter", "Longlands"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"gfayters@desdev.cn", "blonglandst@tripod.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male", "Genderqueer"}), "gender")
                .setColumn(new StringColumn(new String[]{"26.120.76.78/32", "14.92.3.30/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false, false}), "bool")
                .setColumn(new StringColumn(new String[]{"Sweden", "France"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2009-10-02",
                        "2016-07-10")), "date")
                .setColumn(new FloatColumn(new Float[]{595.22f, 879.94f}), "some_number")
                .build();

        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \">28\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "id", ">28", "int")
                .addFuncCallOptionsPattern("id", ">28", ">",
                        null, null, 28)
                .build();
        // input: string id = ">=29" {pattern: int}
        FuncCall funcCall3 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \">=29\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "id", ">=29", "int")
                .addFuncCallOptionsPattern("id", ">=29", ">=",
                        null, null, 29)
                .build();
        // --input: string id = "<=1" {pattern: int}
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{1}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Burk"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Kemery"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"bkemery0@businesswire.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"249.64.22.121/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bool")
                .setColumn(new StringColumn(new String[]{"China"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2017-09-20")), "date")
                .setColumn(new FloatColumn(new Float[]{510.32f}), "some_number")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"<=1\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "", "id", "<=1", "int")
                .addFuncCallOptionsPattern("id", "<=1", "<=",
                        null, null, 1)
                .build();
        // --input: string id = "<2" {pattern: int}
        FuncCall funcCall5 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"<2\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id)\n"
                        + "--end")
                .addFuncParam("string", "", "id", "<2", "int")
                .addFuncCallOptionsPattern("id", "<2", "<",
                        null, null, 2)
                .build();
        // --input: string id = "in(29, 30)" {pattern: int}
        FuncCall funcCall6 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"in(29, 30)\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","id", "in(29, 30)", "int")
                .addFuncCallOptionsPattern("id", "in(29, 30)", "in",
                        null, null, 29, 30)
                .build();
        // --input: string id = "not in(11, 12, 13, 14, 15, "
        //                                + "16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)" {pattern: int}
        DataFrame expected4 = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new IntColumn(
                                Stream.iterate(1, x -> x + 1).limit(10)
                                        .toArray(Integer[]::new)),
                        "id")
                .setColumn(new StringColumn(new String[]{"Burk", "Nicholle", "Orlando", "Gothart",
                        "Mitchell", "Jeromy", "Joela", "Darren", "Marlie", "Scottie"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Kemery", "Karoly", "Westgate",
                                "Cokayne", "Haglington", "Twinn", "Cornau", "Juares", "Mayze", "Formilli"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"bkemery0@businesswire.com", "nkaroly1@alexa.com",
                        "owestgate2@dedecms.com", "gcokayne3@plala.or.jp", "mhaglington4@indiegogo.com",
                        "jtwinn5@globo.com", "jcornau6@imgur.com", "djuares7@hexun.com",
                        "mmayze8@google.com.au", "sformilli9@aol.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male", "Female", "Polygender", "Male", "Male", "Male",
                        "Female", "Male", "Female", "Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"249.64.22.121/32", "255.233.247.118/32",
                                "75.0.252.254/32", "196.83.12.163/32", "209.93.181.190/32", "25.13.2.132/32",
                                "195.47.88.236/32", "94.170.16.96/32", "68.41.25.65/32", "101.241.191.228/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true, false, false, true, true, true,
                        false, false, false, false}), "bool")
                .setColumn(new StringColumn(new String[]{"China", "Poland", "Netherlands",
                        "Philippines", "Poland", "Serbia", "Indonesia", "China",
                        "France", "Vietnam"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern,"2017-09-20",
                        "2014-02-27", "2020-09-03", "2001-01-31", "2020-10-09",
                        "2014-10-04", "2020-03-19", "2011-04-09", "2011-11-10", "2003-01-04")), "date")
                .setColumn(new FloatColumn(new Float[]{510.32f, 864.09f, 822.7f, 251.05f, 15.22f,
                        378.4f, 349.11f, 631.89f, 561.72f, 978.01f}), "some_number")
                .build();
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"not in(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, "
                        + "21, 22, 23, 24, 25, 26, 27, 28, 29, 30)\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "id", "not in(11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, " +
                                "22, 23, 24, 25, 26, 27, 28, 29, 30)",
                        "int")
                .addFuncCallOptionsPattern("id", "not in(11, 12, 13, 14, 15, "
                                + "16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)",
                        "not in", null, null, 11, 12, 13, 14, 15, 16, 17,
                        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
                .build();
        // --input: string id = "min-max 29-30" {pattern: int}
        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string id = \"min-max 29-30\" {pattern: int}\n"
                        + "SELECT * FROM mock_data WHERE @id(id) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "id", "min-max 29-30",
                        "int")
                .addFuncCallOptionsPattern("id", "29-30",
                        "-", null, null, 29, 30)
                .build();
        //--input: double some_number = 510.32
        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: double some_number = 510.32\n"
                        + "SELECT * FROM mock_data WHERE some_number = @some_number\n"
                        + "--end")
                .addFuncParam("double", "", "some_number", 510.32, "double")
                .build();
        // --input: string some_number = ">975" {pattern: double}
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new IntColumn(new Integer[]{10, 26}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Scottie", "Daryle"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Formilli", "O'Shaughnessy"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"sformilli9@aol.com", "doshaughnessyp@com.com"}),
                        "email")
                .setColumn(new StringColumn(new String[]{"Male", "Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"101.241.191.228/32", "204.107.16.207/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false, false}), "bool")
                .setColumn(new StringColumn(new String[]{"Vietnam", "Honduras"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2003-01-04",
                        "2010-05-04")), "date")
                .setColumn(new FloatColumn(new Float[]{978.01f, 983.03f}), "some_number")
                .build();
        FuncCall funcCall10 =  FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \">975\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string",  "", "some_number", ">975", "double")
                .addFuncCallOptionsPattern("some_number", ">975", ">",
                        null, null, 975)
                .build();
        // --input: string some_number = ">=975" {pattern: double}
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \">=975\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "some_number", ">=975", "double")
                .addFuncCallOptionsPattern("some_number", ">=975", ">=",
                        null, null, 975)
                .build();
        //--input: string some_number = "<20" {pattern: double}
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{5}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Mitchell"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Haglington"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"mhaglington4@indiegogo.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"209.93.181.190/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bool")
                .setColumn(new StringColumn(new String[]{"Poland"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2020-10-09")), "date")
                .setColumn(new FloatColumn(new Float[]{15.22f}), "some_number")
                .build();
        FuncCall funcCall12 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \"<20\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "some_number", "<20", "double")
                .addFuncCallOptionsPattern("some_number", "<20", "<",
                        null, null, 20)
                .build();
        // --input: string some_number = "<=20" {pattern: double}
        FuncCall funcCall13 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string some_number = \"<=20\" {pattern: double}\n"
                        + "SELECT * FROM mock_data WHERE @some_number(some_number) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "", "some_number", "<=20", "double")
                .addFuncCallOptionsPattern("some_number", "<=20", "<=",
                        null, null, 20)
                .build();
        // --input: string first_name = 'contains Z' {pattern: string}
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{25}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Zolly"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Wimmers"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"zwimmerso@hatena.ne.jp"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"123.12.225.114/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Bosnia and Herzegovina"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2003-02-12")), "date")
                .setColumn(new FloatColumn(new Float[]{217.18f}), "some_number")
                .build();
        FuncCall funcCall14 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = 'contains Z' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "contains Z", "string")
                .addFuncCallOptionsPattern("first_name", "contains Z", "contains",
                        null, null, "Z")
                .build();
        // --input: string first_name = 'starts with W' {pattern: string}
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{23}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Waly"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Rogliero"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"wroglierom@berkeley.edu"}), "email")
                .setColumn(new StringColumn(new String[]{"Female"}), "gender")
                .setColumn(new StringColumn(new String[]{"122.90.196.231/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{true}), "bool")
                .setColumn(new StringColumn(new String[]{"Sweden"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2011-12-18")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{147.69f}), "some_number")
                .build();
        FuncCall funcCall15 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = 'starts with W' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "starts with W", "string")
                .addFuncCallOptionsPattern("first_name", "starts with W", "starts with",
                        null, null, "W")
                .build();
        // --input: string first_name = 'ends with y' {pattern: string}
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{20}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Lucius"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Edelmann"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"ledelmannj@bravesites.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Male"}), "gender")
                .setColumn(new StringColumn(new String[]{"66.174.30.225/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Brazil"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "1999-06-22")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{378.73f}), "some_number")
                .build();
        FuncCall funcCall16 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = 'ends with s' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name) ORDER BY id\n"
                        + "--end")
                .addFuncParam("string", "","first_name", "ends with s", "string")
                .addFuncCallOptionsPattern("first_name", "ends with s", "ends with",
                        null, null, "s")
                .build();
        // --input: string country = 'in (Poland, Brazil)' {pattern: string}
        DataFrame expected10 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{2, 5, 20}),
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
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2014-02-27",
                                "2020-10-09","1999-06-22")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{864.09f, 15.22f, 378.73f}), "some_number")
                .build();
        FuncCall funcCall17 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string country = 'in (Poland, Brazil)' {pattern: string}\n" +
                        "SELECT * FROM mock_data WHERE @country(country) ORDER BY id\n" +
                        "--end")
                .addFuncParam("string", "", "country", "in (Poland, Brazil)", "string")
                .addFuncCallOptionsPattern("country", "in (Poland, Brazil)", "in",
                        null, null, "Poland", "Brazil")
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: int; operator: =; pattern: none", funcCall1), expected1),
                Arguments.of(Named.of("type: string; operator: >; pattern: int", funcCall2), expected2),
                Arguments.of(Named.of("type: string; operator: >=; pattern: int", funcCall3), expected2),
                Arguments.of(Named.of("type: string; operator: <=; pattern: int", funcCall4), expected3),
                Arguments.of(Named.of("type: string; operator: <; pattern: int", funcCall5), expected3),
                Arguments.of(Named.of("type: string; operator: in; pattern: int", funcCall6), expected2),
                Arguments.of(Named.of("type: string; operator: not in; pattern: int", funcCall7), expected4),
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
        Parser parser = new DateParser();
        String datePattern = "yyyy-MM-dd";
        // --input: string first_name = "starts with p" {pattern: string}
        //--input: string id = ">1" {pattern :int}
        //--input: bool bool = false
        //--input: string email = "contains com" {pattern: string}
        //--input: string some_number = ">20" {pattern: double}
        //--input: string country = "in (Indonesia)" {pattern: string}
        //--input: string date = "before 1/1/2022" {pattern: datetime}
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{13}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Pail"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Boxell"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"pboxellc@moonfruit.com"}), "email")
                .setColumn(new StringColumn(new String[]{"Genderqueer"}), "gender")
                .setColumn(new StringColumn(new String[]{"2.37.160.155/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"Indonesia"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles(datePattern, "2012-01-14")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{73.47f}), "some_number")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.getBuilder()
                .addQuery("--input: string first_name = \"starts with p\" {pattern: string}\n"
                        + "--input: string id = \">1\" {pattern: int}\n"
                        + "--input: string email = \"contains com\" {pattern: string}\n"
                        + "--input: string some_number = \">20\" {pattern: double}\n"
                        + "--input: string country = \"in (Indonesia)\" {pattern: string}\n"
                        + "--input: string dat = \"before 1/1/2022\" {pattern: datetime}\n"
                        + "SELECT * FROM mock_data WHERE @first_name(first_name) AND @id(id) "
                        + "AND @email(email) AND @some_number(some_number) "
                        + "AND @country(country) AND @dat(dat)\n"
                        + "--end")
                .addFuncParam("string", "", "first_name", "starts with p", "string")
                .addFuncParam("string", "","id", ">1", "string")
                .addFuncParam("string", "","email", "contains com", "string")
                .addFuncParam("string", "","some_number", ">20", "double")
                .addFuncParam("string", "","country", "in (Indonesia)", "string")
                .addFuncParam("string", "","dat", "before 1/1/2022", "datetime")
                .addFuncCallOptionsPattern("first_name", "starts with p",
                        "starts with", null, null, "p")
                .addFuncCallOptionsPattern("id", ">1", ">", null,
                        null, "1")
                .addFuncCallOptionsPattern("email", "contains com",
                        "contains", null, null, "com")
                .addFuncCallOptionsPattern("some_number", ">20", ">", null,
                        null, 20)
                .addFuncCallOptionsPattern("country", "in (Indonesia)", "in",
                        null, null, "Indonesia")
                .addFuncCallOptionsPattern("dat", "before 1/1/2022", "before",
                        true, true, Year.of(2022).atMonth(1).atDay(1).toString())
                .build();
        return Stream.of(Arguments.of(Named.of("type: multiple; operator: multiple; pattern: multiple", funcCall1),
                expected1));
    }

    public static Stream<Arguments> checkListParameterSupport_ok() {
        // list<string>
        List<String> values = new ArrayList<>();
        values.add("Poland");
        values.add("Brazil");
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{2, 5, 20}),
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
                        "SELECT * FROM mock_data WHERE country = ANY(@values) ORDER BY id")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> values = ['Poland','Brazil']\n" +
                        "SELECT * FROM mock_data WHERE country = ANY(@values) ORDER BY id")
                .addFuncParam("list", "string", "values", values, "")
                .addFuncCallOptionsPattern("country", "", "",
                        null, null, "Poland", "Brazil")
                .build();
        return Stream.of(
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall1), expected),
                Arguments.of(Named.of("type: list<string>; operator: none; pattern: none", funcCall2), expected));
    }

    public static Stream<Arguments> checkRegexSupport_ok() {
        // --input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new IntColumn(new Integer[]{9}),
                        "id")
                .setColumn(new StringColumn(new String[]{"Marlie"}), "first_name")
                .setColumn(new StringColumn(new String[]{"Mayze"}),
                        "last_name")
                .setColumn(new StringColumn(new String[]{"mmayze8@google.com.au"}), "email")
                .setColumn(new StringColumn(new String[]{"Female"}), "gender")
                .setColumn(new StringColumn(new String[]{"68.41.25.65/32"}),
                        "ip_address")
                .setColumn(new BoolColumn(new Boolean[]{false}), "bool")
                .setColumn(new StringColumn(new String[]{"France"}), "country")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd", "2011-11-10")),
                        "date")
                .setColumn(new FloatColumn(new Float[]{561.72f}), "some_number")
                .build();
        FuncCall funcCall = FuncCallBuilder.getBuilder()
                .addQuery("-- input: string email = 'regex ^([A-Za-z0-9_]+@google.com.au)$' {pattern: string}\n"
                        + "SELECT * FROM mock_data WHERE @email(email)\n"
                        + "-- end")
                .addFuncParam("string", "", "email", "regex ^([A-Za-z0-9_]+@google.com.au)$", "string")
                .addFuncCallOptionsPattern("email", "regex ^([A-Za-z0-9_]+@google.com.au)$",
                        "regex", null, null, "^([A-Za-z0-9_]+@google.com.au)$")
                .build();
        return Stream.of(Arguments.of(Named.of("type: string; operator: regex; pattern: string",
                funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_xmlType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[]{"<foo>Hello World!</foo>\n",
                        "<book>\n  <title>Manual</title>\n  <chapter>...</chapter>\n</book>\n"}), "DATA")
                .build();
        return Stream.of(Arguments.of(
                        Named.of("XML TYPE SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM xml_data")), expected
                )
        );
    }

    public static Stream<Arguments> checkOutputDataFrame_characterTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"Hello"}), "CH")
                .setColumn(new StringColumn(new String[]{"World"}), "VARCH")
                .setColumn(new StringColumn(new String[]{"Datagrok"}), "NCH")
                .setColumn(new StringColumn(new String[]{"Groking"}), "NVARCH")
                .build();
        return Stream.of(Arguments.of(Named.of("CHARACTER TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT TRIM(CH) AS CH, VARCH, TRIM(NCH) AS NCH, "
                        + "NVARCH FROM character_type")), // Use trim function, because of how values stored in char and varchar columns
                expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        Parser parser = new DateParser();
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("dd-MM-yyyy",
                                "01-01-2023")}),
                        "DAT")
                .setColumn(new DateTimeColumn(new Double[]{
                        parser.parseDateToDouble(
                                "dd-MMM-yy HH:mm:ss.SS", "03-AUG-17 11:20:30.450 AM")}), "STAMP")
                .setColumn(new DateTimeColumn(new Double[]{
                        parser.parseDateToDouble("dd-MMM-yyyy HH:mm:ss X", "21-FEB-2009 18:00:00 -05:00")}),
                        "ZONED_STAMP")
                .setColumn(new StringColumn(new String[]{"10-2"}), "INTERVAL1")
                .setColumn(new StringColumn(new String[]{"4 5:12:10.222"}), "INTERVAL2")
                .build();
        return Stream.of(Arguments.of(
                Named.of("DATE TYPES SUPPORT", FuncCallBuilder.fromQuery("SELECT * FROM dates_type")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_jsonType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"{\"phones\":[{\"type\":\"mobile\",\"phone\":\"001001\"},"
                        + "{\"type\":\"fix\",\"phone\":\"002002\"}]}", "{\"bar\":\"baz\",\"balance\":7.77,\"active\":false}",
                        "{\"reading\":0.0000123}"}), "DATA")
                .build();
        return Stream.of(Arguments.of(Named.of("JSON TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM JSON_DATA")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_uriType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"/home/oe/doc1.xml",
                        "/HR/EMPLOYEES/ROW[EMPLOYEE_ID=205]/SALARY", "https://datagrok.ai"}), "URI")
                .build();
        // should be used only with .getURI() otherwise unreadable
        return Stream.of(Arguments.of(Named.of("URI TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT u.uri.getURL() AS URI FROM uri_types u")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_varrayType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"[Brenda, Richard]"}), "MEMBERS")
                .build();
        return Stream.of(Arguments.of(Named.of("VARRAY TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM varrays")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new FloatColumn(new Float[]{1.99f, 1.9f, 13.0f}), "NUMBER_VALUE")
                .setColumn(new IntColumn(new Integer[]{1123, 1, 1244124}), "SMALL_VALUE")
                .setColumn(new FloatColumn(new Float[]{1.2e-4f, 0.55f, 12.123f}), "FLOAT_VALUE")
                .setColumn(new FloatColumn(new Float[]{1.17549E-38f, Float.POSITIVE_INFINITY,
                        Float.NEGATIVE_INFINITY}), "BINARY_FLOAT_VALUE")
                .setColumn(new FloatColumn(new Float[]{0.0f, 0.2222f, 2.2222f}), "BINARY_DOUBLE_VALUE")
                .build();
        return Stream.of(Arguments.of(Named.of("NUMERIC TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM numeric_type")), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_lobsTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"Grok"}), "BLOB_TYPE")
                .setColumn(new StringColumn(new String[]{"Grok"}), "CLOB_TYPE")
                .setColumn(new StringColumn(new String[]{"Grok"}), "NCLOB_TYPE")
                .build();
        return Stream.of(Arguments.of(Named.of("LOB TYPES SUPPORT",
                FuncCallBuilder.fromQuery("SELECT UTL_RAW.CAST_TO_VARCHAR2(BLOB_TYPE) AS BLOB_TYPE, TO_CHAR(CLOB_TYPE) AS CLOB_TYPE, "
                        + "TO_CHAR(NCLOB_TYPE) AS NCLOB_TYPE FROM lobs")), expected));
    }
}
