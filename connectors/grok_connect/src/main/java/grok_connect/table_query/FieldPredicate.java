package grok_connect.table_query;

import grok_connect.utils.GrokConnectUtil;
import grok_connect.utils.PatternMatcher;

import java.util.Arrays;

public class FieldPredicate {
    public String field;
    public String dataType;
    public String pattern;
    public PatternMatcher matcher;

    public FieldPredicate(String field, String pattern, String dataType) {
        this.field = field;
        this.pattern = pattern;
        this.dataType = dataType;
    }

    public String getParamName() {
        String[] split = field
                .replaceAll("\\.", "_")
                .replaceAll(" ", "_")
                .toLowerCase()
                .split("_");
        if (split.length > 1)
            for (int i = 1; i < split.length; i++)
                split[i] = GrokConnectUtil.capitalize(split[i]);

        return String.join("", split);
    }

    @Override
    public String toString() {
        return field + " " + pattern;
    }
}
