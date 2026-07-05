package grok_connect.utils;

import com.google.gson.*;
import grok_connect.connectors_info.DataQuery;
import grok_connect.table_mutation.DeleteRows;
import grok_connect.table_mutation.InsertRows;
import grok_connect.table_mutation.MutationBatch;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_mutation.UpdateRows;
import grok_connect.table_mutation.UpsertRows;
import grok_connect.table_query.TableQuery;
import java.lang.reflect.Type;

public class DataQueryDeserializer implements JsonDeserializer<DataQuery> {
    private static final Gson gson = new Gson();

    static Class<? extends TableMutation> getMutationClass(String type) {
        switch (type) {
            case "InsertRows": return InsertRows.class;
            case "UpdateRows": return UpdateRows.class;
            case "DeleteRows": return DeleteRows.class;
            case "UpsertRows": return UpsertRows.class;
            case "MutationBatch": return MutationBatch.class;
            default: return null;
        }
    }

    @Override
    public DataQuery deserialize(JsonElement jsonElement, Type type, JsonDeserializationContext jsonDeserializationContext) throws JsonParseException {
        String myType = jsonElement.getAsJsonObject().get("#type").getAsString();
        Class<? extends TableMutation> mutationClass = getMutationClass(myType);
        if (mutationClass != null)
            return jsonDeserializationContext.deserialize(jsonElement, mutationClass);
        switch (myType) {
            case "TableQuery": return jsonDeserializationContext.deserialize(jsonElement, TableQuery.class);
            case "DataQuery": return gson.fromJson(jsonElement, DataQuery.class);
            default: return null;
        }
    }
}
