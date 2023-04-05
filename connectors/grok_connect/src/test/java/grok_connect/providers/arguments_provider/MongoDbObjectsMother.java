package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.FuncCallBuilder;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BoolColumn;
import serialization.DataFrame;
import serialization.FloatColumn;
import serialization.StringColumn;
import java.util.stream.Stream;

public class MongoDbObjectsMother {
    public static Stream<Arguments> checkOutputAllTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(6)
                .setColumn(new StringColumn(new String[]{"642d3cf9675463d439e9d23d",
                        "642d3cf9c355b25fe3f3268a", "642d3cf936ac9251742be958",
                        "642d3cf975fe1b1ac7e7933a", "642d3cf9abd475f558e55197", "642d3cf906dcdfea7bb8bf34"}),
                        "_id")
                .setColumn(new FloatColumn(new Float[]{0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f}), "index")
                .setColumn(new StringColumn(new String[]{"1606b22d-f927-418c-a5aa-f03f14f9ce48",
                                "9a45da61-6f43-4512-9d2e-9dca26c4ff56", "b92ad6a7-4d2f-4646-a42d-538085e29935",
                                "218bac3c-ddff-44fe-8429-192dfb235982", "a8b57021-1203-4a1c-b907-1c2fcd06519e",
                                "72516057-8b78-47c9-83ab-c7079c6184e1"}),
                        "guid")
                .setColumn(new BoolColumn(new Boolean[]{false, true, false, false, true, false}), "isActive")
                .setColumn(new StringColumn(new String[]{"$2,679.87", "$3,609.65", "$2,155.94", "$1,108.98",
                        "$3,429.22", "$1,061.53"}), "balance")
                .setColumn(new StringColumn(new String[]{"http://placehold.it/32x32",
                        "http://placehold.it/32x32", "http://placehold.it/32x32", "http://placehold.it/32x32",
                        "http://placehold.it/32x32", "http://placehold.it/32x32"}), "picture")
                .setColumn(new FloatColumn(new Float[]{28.0f, 25.0f, 31.0f, 25.0f, 33.0f, 21.0f}), "age")
                .setColumn(new FloatColumn(new Float[]{47.66397f, 82.178986f, -65.5806f, -2.464041f, 82.45335f, 74.18999f}), "latitude")
                .setColumn(new FloatColumn(new Float[]{-86.55363f, -117.3157f, -31.680079f, -76.23148f, 74.44026f, 51.04572f}), "longitude")
                .setColumn(new StringColumn(new String[]{"[esse, dolor, ex, excepteur, consectetur, et, ad]",
                "[mollit, laborum, laborum, do, Lorem, dolore, minim]", "[irure, ad, irure, occaecat, esse, amet, excepteur]",
                "[do, anim, labore, dolore, consequat, occaecat, ad]", "[deserunt, dolore, ad, eu, minim, deserunt, sint]",
                "[duis, ex, aliqua, cillum, duis, aliqua, amet]"}), "tags")
                .setColumn(new StringColumn(new String[]{"[Document{{id=0.0, name=Mccoy Lindsey}}, Document{{id=1.0, name=Joanna Lynch}}, Document{{id=2.0, name=Evangeline Combs}}]",
                "[Document{{id=0.0, name=Anita Nielsen}}, Document{{id=1.0, name=Campbell Hatfield}}, Document{{id=2.0, name=Mccormick Simon}}]",
                        "[Document{{id=0.0, name=Denise Curry}}, Document{{id=1.0, name=Joyner Soto}}, Document{{id=2.0, name=Chapman Baldwin}}]",
                        "[Document{{id=0.0, name=Jerri Poole}}, Document{{id=1.0, name=Minnie Hood}}, Document{{id=2.0, name=Talley Moody}}]",
                "[Document{{id=0.0, name=Wilkins Vaughan}}, Document{{id=1.0, name=Bettie Carson}}, Document{{id=2.0, name=House Pope}}]",
                        "[Document{{id=0.0, name=Nieves Hobbs}}, Document{{id=1.0, name=Fisher Vasquez}}, Document{{id=2.0, name=Davidson Norman}}]"}), "friends")
                .build();
        return Stream.of(Arguments.of(Named.of("ALL TYPES",
                FuncCallBuilder.fromQuery("db.mocks.find()")), expected));
    }

    public static Stream<Arguments> checkStringReturn_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new StringColumn(new String[]{"admin", "config", "local", "test"}),
                        "DATABASE_NAME")
                .build();
        return Stream.of(Arguments.of(Named.of("STRING RETURN TYPE",
                FuncCallBuilder.fromQuery("show dbs")), expected));
    }

    public static Stream<Arguments> checkOneDocument_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new StringColumn(new String[]{"642d3cf9675463d439e9d44d"}), "_id")
                .setColumn(new StringColumn(new String[]{"Gomez"}),
                        "name")
                .setColumn(new BoolColumn(new Boolean[]{true}), "isActive")
                .setColumn(new StringColumn(new String[]{"$1,316.05"}),
                        "balance")
                .build();
        return Stream.of(Arguments.of(Named.of("DOCUMENT RETURN TYPE",
                FuncCallBuilder.fromQuery("db.one_line.find()")), expected));
    }

    public static Stream<Arguments> checkNullResult_ok() {
        return Stream.of(Arguments.of(Named.of("NULL RETURN",
                FuncCallBuilder.fromQuery("db.some_foo_bar.find()")), new DataFrame()));
    }
}
