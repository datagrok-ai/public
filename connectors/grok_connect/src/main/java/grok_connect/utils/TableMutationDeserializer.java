package grok_connect.utils;

import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonParseException;
import grok_connect.table_mutation.TableMutation;
import java.lang.reflect.Type;

/**
 * Polymorphic adapter for fields declared as TableMutation (MutationBatch.operations);
 * Gson does not reuse the DataQuery adapter for subclasses of DataQuery.
 */
public class TableMutationDeserializer implements JsonDeserializer<TableMutation> {
    @Override
    public TableMutation deserialize(JsonElement jsonElement, Type type, JsonDeserializationContext jsonDeserializationContext) throws JsonParseException {
        JsonElement typeElement = jsonElement.getAsJsonObject().get("#type");
        if (typeElement == null || typeElement.isJsonNull())
            throw new JsonParseException("Missing #type in TableMutation JSON");
        String myType = typeElement.getAsString();
        Class<? extends TableMutation> mutationClass = DataQueryDeserializer.getMutationClass(myType);
        if (mutationClass == null)
            throw new JsonParseException("Unknown TableMutation type: " + myType);
        return jsonDeserializationContext.deserialize(jsonElement, mutationClass);
    }
}
