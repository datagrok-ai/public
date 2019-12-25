package grok_connect.utils;

import java.util.*;
import grok_connect.connectors_info.*;


public class PatternMatcherResult {
    public List<FuncParam> params;
    public String query;

    public PatternMatcherResult() {
        params = new ArrayList<>();
    }
}
