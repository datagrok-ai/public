package grok_connect.utils;

import java.util.*;
import grok_connect.connectors_info.*;

public class PatternMatcherResult {
    public List<FuncParam> params;
    public String query;

    public PatternMatcherResult() {
        params = new ArrayList<>();
    }

    public List<FuncParam> getParams() {
        return params;
    }

    public void addParam(FuncParam param) {
        this.params.add(param);
    }

    public String getQuery() {
        return query;
    }

    public void setQuery(String query) {
        this.query = query;
    }
}
