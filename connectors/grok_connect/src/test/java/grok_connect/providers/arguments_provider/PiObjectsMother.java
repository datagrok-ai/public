package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.StringColumn;
import java.time.Year;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class PiObjectsMother {
    public static final Parser parser = new DateParser();

    public static Stream<Arguments> testQueries() {
        String in1 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Shape (k)";
        String in2 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Scale (Theta)";
        String in3 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Initial Volume";
        String in4 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Sum of Squares";
        String in5 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Mean";
        String in6 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-HETP";
        String in7 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Asymmetric Factor";
        String in8 = "033_PI-OPCHDA-UNICORN_ColumnVolumeRetr";
        String in9 = "033_PI-OPCHDA-UNICORN_ColumnDiameter";
        String in10 = "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Bed Height";
        DataFrame expected1 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"033_PI-OPCHDA-UNICORN_ColumnDiameter",
                        "033_PI-OPCHDA-UNICORN_ColumnVolumeRetr", "033_PI-OPCHDA-UNICORN_ColumnVolumeRetr"}), "tag")
                .setColumn(new DateTimeColumn(new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS X",
                                        "2022-01-19 13:59:35.000 +00:00"), parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS X",
                                "2021-12-07 20:56:23.106 +00:00"), parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS X",
                                "2021-12-09 17:44:45.270 +00:00")}),
                        "time")
                .setColumn(new FloatColumn(new Float[]{null, 0.01005f, 0.01005f}), "value")
                .build();
        FuncCall funcCall1 = FuncCallBuilder.fromQuery("select top 3 tag,time,value from piarchive..picomp calc where tag in (\n" +
                String.format("'%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s'\n", in1, in2, in3, in4,
                        in5, in6, in7, in8, in9, in10) +
                ") order by tag,time");
        DataFrame expected2 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 21:02:24.917 +00:00", "2021-12-07 21:11:19.816 +00:00",
                        "2021-12-07 21:20:14.116 +00:00")),"time")
                .setColumn(new StringColumn(new String[]{"Phase 25mM Bis-Tris, 80mM NaCl, pH 6 ended.",
                                "Phase 25mM Bis-Tris, 1M NaCl, pH 6 ended.", "Phase 25mM Bis-Tris, 80mM NaCl, pH 6 ended."}),
                        "svalue")
                .build();
        FuncCall funcCall2 = FuncCallBuilder.fromQuery("select top 3 time,svalue from piarchive..picomp "
                + "where tag = '033_PI-OPCHDA-UNICORN_PhaseEnd' order by time");

        FuncCall funcCall3 = FuncCallBuilder.getBuilder().addQuery("--input: string tag = \"033_PI-OPCHDA-UNICORN_PhaseEnd\"\n"
                        + "select top 3 time,svalue from piarchive..picomp "
                        + "where tag = @tag order by time")
                .addFuncParam("string", "", "tag", "033_PI-OPCHDA-UNICORN_PhaseEnd",
                        "")
                .build();
        FuncCall funcCall4 = FuncCallBuilder.getBuilder().addQuery("--input: string tag = \"033_PI-OPCHDA-UNICORN_PhaseEnd\" {pattern: string}\n"
                        + "select top 3 time,svalue from piarchive..picomp "
                        + "where @tag(tag) order by time")
                .addFuncParam("string", "", "tag", "033_PI-OPCHDA-UNICORN_PhaseEnd",
                        "string")
                .addFuncCallOptionsPattern("tag", "033_PI-OPCHDA-UNICORN_PhaseEnd", "equals",
                        null, null, "033_PI-OPCHDA-UNICORN_PhaseEnd")
                .build();
        List<String> tags = new ArrayList<>();
        tags.add(in1);
        tags.add(in2);
        tags.add(in3);
        tags.add(in4);
        tags.add(in5);
        tags.add(in6);
        tags.add(in7);
        tags.add(in8);
        tags.add(in9);
        tags.add(in10);
        FuncCall funcCall5 = FuncCallBuilder.getBuilder().addQuery("--input: list<string> tags\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc "
                        + "where tag in (@tags) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .build();
        DataFrame expected3 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Shape (k)"}), "tag")
                .setColumn(new DateTimeColumn(new Double[]{
                                parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS X",
                                        "2021-12-09 18:19:48.221 +00:00")}),
                        "time")
                .setColumn(new FloatColumn(new Float[]{43.135532f}), "value")
                .build();
        FuncCall funcCall6 = FuncCallBuilder.getBuilder().addQuery("--input: list<string> tags\n"
                        + "--input: string value = \">42\" {pattern: double}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc "
                        + "where tag in (@tags) and @value(value) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncParam("string", "", "value", 42,
                        "double")
                .addFuncCallOptionsPattern("value", ">42", ">",
                        null, null, 42)
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .build();
        DataFrame expected4 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"033_PI-OPCHDA-UNICORN_ColumnVolumeRetr",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Asymmetric Factor",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Asymmetric Factor"}), "tag")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 20:56:23.106 +00:00", "2021-12-07 21:00:30.016 +00:00",
                        "2021-12-07 21:09:24.916 +00:00")), "time")
                .setColumn(new FloatColumn(new Float[]{0.01005f, 1.3200065f, 1.5109893f}), "value")
                .build();
        FuncCall funcCall7 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string time = \"before 2021-12-08\" {pattern: datetime}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @time(time) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "","time", "before 2021-12-08", "datetime")
                .addFuncCallOptionsPattern("time", "", "before", true, true,
                        Year.of(2021).atMonth(12).atDay(8).toString())
                .build();
        DataFrame expected5 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"033_PI-OPCHDA-UNICORN_ColumnVolumeRetr",
                        "033_PI-OPCHDA-UNICORN_ColumnVolumeRetr",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Asymmetric Factor"}), "tag")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 20:56:23.106 +00:00", "2021-12-09 17:44:45.270 +00:00",
                        "2021-12-07 21:00:30.016 +00:00")), "time")
                .setColumn(new FloatColumn(new Float[]{0.01005f, 0.01005f, 1.3200065f}), "value")
                .build();
        FuncCall funcCall8 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string time = \"2021-2022\" {pattern: datetime}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @time(time) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "","time", "2021-2022", "datetime")
                .addFuncCallOptionsPattern("time", "", "range", true, true,
                        Year.of(2021).atDay(1).toString(),
                        Year.of(2022).atDay(1).toString())
                .build();
        DataFrame expected6 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"033_PI-OPCHDA-UNICORN_ColumnDiameter"}), "tag")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2022-01-19 13:59:35.000 +00:00")), "time")
                .setColumn(new FloatColumn(new Float[]{null}), "value")
                .build();
        FuncCall funcCall9 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string time = \"January 2022\" {pattern: datetime}\n"
                        + "select tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @time(time) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "","time", "January 2022", "datetime")
                .addFuncCallOptionsPattern("time", "", "range", true, false,
                        Year.of(2022).atMonth(1).atDay(1).toString(),
                        Year.of(2022).atMonth(2).atDay(1).toString())
                .build();
        FuncCall funcCall10 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string time = \"anytime\" {pattern: datetime}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @time(time) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "","time", "anytime", "datetime")
                .addFuncCallOptionsPattern("time", "", "none", true, true)
                .build();
        DataFrame expected7 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"033_PI-OPCHDA-UNICORN_ColumnVolumeRetr"}), "tag")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 20:56:23.106 +00:00")), "time")
                .setColumn(new FloatColumn(new Float[]{0.01005f}), "value")
                .build();
        FuncCall funcCall11 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string value = \"<0.1\" {pattern: double}\n"
                        + "select top 1 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @value(value) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "", "value", 0.1,
                        "double")
                .addFuncCallOptionsPattern("value", "<0.1", "<",
                        null, null, 0.1)
                .build();
        DataFrame expected8 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-HETP",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-HETP",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-HETP"}), "tag")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 21:00:30.016 +00:00", "2021-12-07 21:09:24.916 +00:00",
                        "2021-12-07 21:18:19.117 +00:00")), "time")
                .setColumn(new FloatColumn(new Float[]{0.20504594f, 0.4594084f, 0.14371714f}), "value")
                .build();
        FuncCall funcCall12 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string value = \"min-max 0.1-1\" {pattern: double}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @value(value) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "", "value", "min-max 29-30",
                        "double")
                .addFuncCallOptionsPattern("value", "0.1-1",
                        "-", null, null, 0.1, 1)
                .build();
        FuncCall funcCall13 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string value = \"not in(0.1)\" {pattern: double}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @value(value) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "", "value", "not in(0.1)",
                        "double")
                .addFuncCallOptionsPattern("value", "not in(0.1)",
                        "not in", null, null, 0.1)
                .build();
        DataFrame expected9 = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[]{"GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Scale (Theta)",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Scale (Theta)",
                        "GDTA-JRD_MVN_AKTA-JRD_MVN_AKTA_Processing-Processing-Optimized Scale (Theta)"}), "tag")
                .setColumn(new DateTimeColumn(parser.parseDatesToDoubles("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 21:09:24.916 +00:00", "2022-03-16 18:23:12.122 +00:00",
                        "2022-03-17 19:02:59.267 +00:00")), "time")
                .setColumn(new FloatColumn(new Float[]{0.1000f, 0.1000f, 0.1000f}), "value")
                .build();
        FuncCall funcCall14 = FuncCallBuilder.getBuilder()
                .addQuery("--input: list<string> tags\n"
                        + "--input: string value = \"in(0.1)\" {pattern: double}\n"
                        + "select top 3 tag,time,value from piarchive..picomp calc where tag in (@tags) "
                        + "and @value(value) order by tag,time")
                .addFuncParam("list", "string", "tags", tags,
                        "")
                .addFuncCallOptionsPattern("tag", "", "",
                        null, null,
                        in1, in2, in3, in4, in5, in6, in7, in8, in9, in10)
                .addFuncParam("string", "", "value", "in(0.1)",
                        "double")
                .addFuncCallOptionsPattern("value", "in(0.1)",
                        "in", null, null, 0.1)
                .build();
        DataFrame expected10 = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SSS X",
                        "2021-12-07 20:56:23.917 +00:00")}), "time")
                .setColumn(new StringColumn(new String[]{"/Malvern/GP140 (Mosaic & Clade C)/gp140 (Mosaic)/(Stage 05) "
                        + "Capto MMC/2020 Reactor Day Purification/GDTA Test/GDTA Test 001 "
                        + "User NA/radelber at computer WCNTUSMJ09MX3X is in control"}), "svalue")
                .build();
        FuncCall funcCall15 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string tag = \"contains _PI-OPCHDA-UNICORN_ResultName\" {pattern: string}\n"
                        + "select top 1 time,svalue from piarchive..picomp where @tag(tag) order by time")
                .addFuncParam("string", "","tag", "contains _PI-OPCHDA-UNICORN_ResultName", "string")
                .addFuncCallOptionsPattern("tag", "contains _PI-OPCHDA-UNICORN_ResultName", "contains",
                        null, null, "_PI-OPCHDA-UNICORN_ResultNam")
                .build();
        FuncCall funcCall16 = FuncCallBuilder.getBuilder()
                .addQuery( "--input: string tag = \"in (033_PI-OPCHDA-UNICORN_ResultName)\" {pattern: string}\n"
                        + "select top 1 time,svalue from piarchive..picomp where @tag(tag) order by time")
                .addFuncParam("string", "", "tag", "in (033_PI-OPCHDA-UNICORN_ResultName)", "string")
                .addFuncCallOptionsPattern("tag", "in (033_PI-OPCHDA-UNICORN_ResultName)", "in",
                        null, null, "033_PI-OPCHDA-UNICORN_ResultName")
                .build();
        return Stream.of(
                Arguments.of(Named.of("test#1", funcCall1), expected1),
                Arguments.of(Named.of("test#2", funcCall2), expected2),
                Arguments.of(Named.of("test#3", funcCall3), expected2),
                Arguments.of(Named.of("test#4", funcCall4), expected2),
                Arguments.of(Named.of("test#5", funcCall5), expected1),
                Arguments.of(Named.of("test#6", funcCall6), expected3),
                Arguments.of(Named.of("test#7", funcCall7), expected4),
                Arguments.of(Named.of("test#8", funcCall8), expected5),
                Arguments.of(Named.of("test#9", funcCall9), expected6),
                Arguments.of(Named.of("test#10", funcCall10), expected1),
                Arguments.of(Named.of("test#11", funcCall11), expected7),
                Arguments.of(Named.of("test#12", funcCall12), expected8),
                Arguments.of(Named.of("test#13", funcCall13), expected1),
                Arguments.of(Named.of("test#14", funcCall14), expected9),
                Arguments.of(Named.of("test#15", funcCall15), expected10),
                Arguments.of(Named.of("test#16", funcCall16), expected10));
    }
}
