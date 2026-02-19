package serialization;

import java.util.Map;

public interface Taggable {
    Map<String, String> getTags();

    default void setTag(String key, String value) {
        getTags().put(key, value);
    }

    default String getTag(String key) {
        return getTags().get(key);
    }

    default boolean hasTag(String key) {
        return getTags().containsKey(key);
    }

    default void removeTag(String key) {
        getTags().remove(key);
    }
}
