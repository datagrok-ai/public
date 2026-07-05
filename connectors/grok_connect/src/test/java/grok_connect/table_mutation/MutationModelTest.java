package grok_connect.table_mutation;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataQuery;
import grok_connect.table_query.FieldPredicate;
import org.apache.commons.io.IOUtils;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Gson round-trip tests for the table_mutation wire contract (GROK-20325).
 * Fixtures in src/test/resources/table_mutation/ pin the exact JSON shape the
 * Dart side (WO-7) emits; no Docker required.
 */
public class MutationModelTest {
    private String readFixture(String name) throws IOException {
        try (InputStream is = getClass().getResourceAsStream("/table_mutation/" + name)) {
            return IOUtils.toString(is, StandardCharsets.UTF_8);
        }
    }

    /**
     * Round-trips a fixture both directions through the DataQuery Gson seam:
     * deserializes to the expected concrete type, asserts the fixture JSON survives
     * serialization unchanged (as a subset — serialization adds defaulted fields),
     * and asserts a second deserialize/serialize cycle is stable.
     */
    private <T extends TableMutation> T roundTrip(String fixtureName, Class<T> expectedClass) throws IOException {
        String json = readFixture(fixtureName);
        DataQuery d1 = GrokConnect.gson.fromJson(json, DataQuery.class);
        Assertions.assertEquals(expectedClass, d1.getClass());

        String serialized = GrokConnect.gson.toJson(d1);
        assertSubset(JsonParser.parseString(json), JsonParser.parseString(serialized), "$");

        DataQuery d2 = GrokConnect.gson.fromJson(serialized, DataQuery.class);
        Assertions.assertEquals(expectedClass, d2.getClass());
        Assertions.assertEquals(JsonParser.parseString(serialized), JsonParser.parseString(GrokConnect.gson.toJson(d2)));
        return expectedClass.cast(d1);
    }

    private void assertSubset(JsonElement expected, JsonElement actual, String path) {
        if (expected.isJsonObject()) {
            Assertions.assertTrue(actual.isJsonObject(), path + ": expected an object");
            JsonObject actualObject = actual.getAsJsonObject();
            for (Map.Entry<String, JsonElement> e : expected.getAsJsonObject().entrySet()) {
                Assertions.assertTrue(actualObject.has(e.getKey()), path + "." + e.getKey() + " lost on serialization");
                assertSubset(e.getValue(), actualObject.get(e.getKey()), path + "." + e.getKey());
            }
        }
        else if (expected.isJsonArray()) {
            Assertions.assertTrue(actual.isJsonArray(), path + ": expected an array");
            JsonArray expectedArray = expected.getAsJsonArray();
            JsonArray actualArray = actual.getAsJsonArray();
            Assertions.assertEquals(expectedArray.size(), actualArray.size(), path + ": array size");
            for (int i = 0; i < expectedArray.size(); i++)
                assertSubset(expectedArray.get(i), actualArray.get(i), path + "[" + i + "]");
        }
        else
            Assertions.assertEquals(expected, actual, path);
    }

    @DisplayName("InsertRows: int, bigint-as-string, bool, null, ISO datetime")
    @Test
    public void insertRows_roundTrip() throws IOException {
        InsertRows insert = roundTrip("insert_rows.json", InsertRows.class);
        Assertions.assertEquals("InsertRows", insert.type);
        Assertions.assertEquals("orders", insert.tableName);
        Assertions.assertEquals("public", insert.schema);
        Assertions.assertEquals(Arrays.asList("id", "big_id", "active", "note", "created"), insert.columns);
        Assertions.assertEquals(Arrays.asList("int", "bigint", "bool", "string", "datetime"), insert.columnTypes);
        Assertions.assertEquals(2, insert.rows.size());
        List<Object> row = insert.rows.get(0);
        // Gson maps JSON numbers in Object fields to Double — the WO-3 binder coerces via columnTypes
        Assertions.assertTrue(row.get(0) instanceof Double);
        Assertions.assertEquals(1, ((Number) row.get(0)).intValue());
        // bigint values beyond 2^53 travel as JSON strings (wire contract)
        Assertions.assertEquals("9007199254740995", row.get(1));
        Assertions.assertEquals(Boolean.TRUE, row.get(2));
        Assertions.assertNull(row.get(3));
        Assertions.assertEquals("2026-07-05T12:34:56.789Z", row.get(4));
        Assertions.assertEquals("insert", insert.mode);
        Assertions.assertEquals(Collections.singletonList("id"), insert.keyColumns);
        Assertions.assertTrue(insert.allOrNothing);
        Assertions.assertFalse(insert.errorOnDuplicate);
        Assertions.assertTrue(insert.returnGeneratedKeys);
        Assertions.assertFalse(insert.bulk);
    }

    @DisplayName("UpsertRows: extends InsertRows, matchKeys required")
    @Test
    public void upsertRows_roundTrip() throws IOException {
        UpsertRows upsert = roundTrip("upsert_rows.json", UpsertRows.class);
        Assertions.assertEquals("UpsertRows", upsert.type);
        Assertions.assertEquals("upsert", upsert.mode);
        Assertions.assertEquals(Collections.singletonList("id"), upsert.matchKeys);
        Assertions.assertEquals(Arrays.asList("id", "qty"), upsert.columns);
        Assertions.assertEquals(20, ((Number) upsert.rows.get(1).get(1)).intValue());
    }

    @DisplayName("UpdateRows: setColumns/setValues/setTypes parallel lists + FieldPredicate matcher passthrough")
    @Test
    public void updateRows_roundTrip() throws IOException {
        UpdateRows update = roundTrip("update_rows.json", UpdateRows.class);
        Assertions.assertEquals("UpdateRows", update.type);
        Assertions.assertEquals(Arrays.asList("status", "qty"), update.setColumns);
        Assertions.assertEquals("shipped", update.setValues.get(0));
        Assertions.assertEquals(5, ((Number) update.setValues.get(1)).intValue());
        Assertions.assertEquals(Arrays.asList("string", "int"), update.setTypes);
        Assertions.assertEquals(update.setColumns.size(), update.setValues.size());
        Assertions.assertEquals(update.setColumns.size(), update.setTypes.size());
        Assertions.assertEquals("and", update.whereOp);
        FieldPredicate predicate = update.whereClauses.get(0);
        Assertions.assertEquals("id", predicate.field);
        Assertions.assertEquals("int", predicate.dataType);
        Assertions.assertEquals(">3", predicate.pattern);
        Assertions.assertEquals(">", predicate.matcher.op);
        Assertions.assertEquals(">3", predicate.matcher.expression);
        Assertions.assertEquals(Collections.singletonList("3"), predicate.matcher.values);
    }

    @DisplayName("DeleteRows: whereOp, allowFullTable safety rail, matcher passthrough")
    @Test
    public void deleteRows_roundTrip() throws IOException {
        DeleteRows delete = roundTrip("delete_rows.json", DeleteRows.class);
        Assertions.assertEquals("DeleteRows", delete.type);
        Assertions.assertEquals("or", delete.whereOp);
        Assertions.assertFalse(delete.allowFullTable);
        FieldPredicate predicate = delete.whereClauses.get(0);
        Assertions.assertEquals("status", predicate.field);
        Assertions.assertEquals("equals", predicate.matcher.op);
        Assertions.assertEquals(Boolean.TRUE, predicate.matcher.include1);
        Assertions.assertEquals(Collections.singletonList("obsolete"), predicate.matcher.values);
    }

    @DisplayName("MutationBatch: nested heterogeneous operations, polymorphic elements")
    @Test
    public void mutationBatch_roundTrip() throws IOException {
        MutationBatch batch = roundTrip("mutation_batch.json", MutationBatch.class);
        Assertions.assertEquals("MutationBatch", batch.type);
        Assertions.assertEquals(4, batch.operations.size());
        Assertions.assertEquals(InsertRows.class, batch.operations.get(0).getClass());
        Assertions.assertEquals(UpdateRows.class, batch.operations.get(1).getClass());
        Assertions.assertEquals(DeleteRows.class, batch.operations.get(2).getClass());
        Assertions.assertEquals(MutationBatch.class, batch.operations.get(3).getClass());
        InsertRows insert = (InsertRows) batch.operations.get(0);
        Assertions.assertEquals("9007199254740997", insert.rows.get(0).get(1));
        MutationBatch nested = (MutationBatch) batch.operations.get(3);
        Assertions.assertEquals(UpsertRows.class, nested.operations.get(0).getClass());
        Assertions.assertEquals(Collections.singletonList("id"), ((UpsertRows) nested.operations.get(0)).matchKeys);
    }

    @DisplayName("TableMutation adapter: direct deserialization and unknown-type error")
    @Test
    public void tableMutationAdapter() throws IOException {
        TableMutation mutation = GrokConnect.gson.fromJson(readFixture("mutation_batch.json"), TableMutation.class);
        Assertions.assertEquals(MutationBatch.class, mutation.getClass());
        Assertions.assertThrows(JsonParseException.class,
                () -> GrokConnect.gson.fromJson("{\"#type\": \"DropTable\"}", TableMutation.class));
    }

    @DisplayName("DataQuery seam: missing #type is a structured error, unknown #type names the type, known types still work")
    @Test
    public void dataQuerySeam_typeGuards() {
        JsonParseException missing = Assertions.assertThrows(JsonParseException.class,
                () -> GrokConnect.gson.fromJson("{\"query\": \"select 1\"}", DataQuery.class));
        Assertions.assertEquals("Missing #type in DataQuery JSON", missing.getMessage());
        Assertions.assertThrows(JsonParseException.class,
                () -> GrokConnect.gson.fromJson("{\"#type\": null, \"query\": \"select 1\"}", DataQuery.class));
        JsonParseException unknown = Assertions.assertThrows(JsonParseException.class,
                () -> GrokConnect.gson.fromJson("{\"#type\": \"BogusQuery\", \"query\": \"select 1\"}", DataQuery.class));
        Assertions.assertEquals("Unknown DataQuery type: BogusQuery", unknown.getMessage());
        // plain queries keep working — the wire always carries #type (Dart PropMixin emits it)
        DataQuery plain = GrokConnect.gson.fromJson("{\"#type\": \"DataQuery\", \"query\": \"select 1\"}", DataQuery.class);
        Assertions.assertEquals(DataQuery.class, plain.getClass());
        Assertions.assertEquals("select 1", plain.query);
    }

    @DisplayName("TableMutation adapter: missing #type is a structured error, not an NPE")
    @Test
    public void tableMutationAdapter_missingType() {
        JsonParseException e = Assertions.assertThrows(JsonParseException.class,
                () -> GrokConnect.gson.fromJson("{\"tableName\": \"orders\"}", TableMutation.class));
        Assertions.assertEquals("Missing #type in TableMutation JSON", e.getMessage());
        Assertions.assertThrows(JsonParseException.class,
                () -> GrokConnect.gson.fromJson("{\"#type\": null, \"tableName\": \"orders\"}", TableMutation.class));
    }
}
