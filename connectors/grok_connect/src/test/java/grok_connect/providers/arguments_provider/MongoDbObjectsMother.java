package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.FuncCallBuilder;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BoolColumn;
import serialization.DataFrame;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.util.stream.Stream;

public class MongoDbObjectsMother {
    public static Stream<Arguments> checkOutputAllTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(5)
                .setColumn(new StringColumn(new String[]{"642d3cf9675463d439e9d23d",
                        "642d3cf9c355b25fe3f3268a", "642d3cf936ac9251742be958",
                        "642d3cf975fe1b1ac7e7933a", "642d3cf9abd475f558e55197"}),
                        "map._id")
                .setColumn(new IntColumn(new Integer[]{0, 1, 2, 3, 4}), "map.index")
                .setColumn(new StringColumn(new String[]{"1606b22d-f927-418c-a5aa-f03f14f9ce48",
                                "9a45da61-6f43-4512-9d2e-9dca26c4ff56", "b92ad6a7-4d2f-4646-a42d-538085e29935",
                                "218bac3c-ddff-44fe-8429-192dfb235982", "a8b57021-1203-4a1c-b907-1c2fcd06519e"}),
                        "map.guid")
                .setColumn(new BoolColumn(new Boolean[]{false, true, false, false, true}), "map.isActive")
                .setColumn(new StringColumn(new String[]{"$2,679.87", "$3,609.65", "$2,155.94", "$1,108.98",
                        "$3,429.22"}), "map.balance")
                .setColumn(new StringColumn(new String[]{"http://placehold.it/32x32",
                        "http://placehold.it/32x32", "http://placehold.it/32x32", "http://placehold.it/32x32",
                        "http://placehold.it/32x32"}), "map.picture")
                .setColumn(new IntColumn(new Integer[]{28, 25, 31, 25, 33}), "map.age")
                .setColumn(new FloatColumn(new Float[]{47.66397f, 82.178986f, -65.5806f, -2.464041f, 82.45335f}), "map.latitude")
                .setColumn(new FloatColumn(new Float[]{-86.55363f, -117.3157f, -31.680079f, -76.23148f, 74.44026f}), "map.longitude")
                .setColumn(new StringColumn(new String[]{"[esse, dolor, ex, excepteur, consectetur, et, ad]",
                "[mollit, laborum, laborum, do, Lorem, dolore, minim]", "[irure, ad, irure, occaecat, esse, amet, excepteur]",
                "[do, anim, labore, dolore, consequat, occaecat, ad]", "[deserunt, dolore, ad, eu, minim, deserunt, sint]"}),
                        "map.tags")
                .setColumn(new IntColumn(new Integer[] {0, 1, 2, 0, 0}), "map.friends.is")
                .setColumn(new StringColumn(new String[]{"Mccoy Lindsey", "Joanna Lynch",
                        "Evangeline Combs", "Anita Nielsen", "Denise Curry", "Jerri Poole", "Wilkins Vaughan"}), "map.friends.name")
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
