package grok_connect.connectors_info;

import java.util.*;
import org.apache.commons.lang3.text.*;


public class DataQuery
{
    public String id;
    public String name;
    public String query;
    public String connectionId;
    public DataConnection connection;
    private List<FuncParam> params;

    public Map<String, Object> parameters = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);

    public String get(String parameter) {
        return (String)parameters.get(WordUtils.capitalize(parameter));
    }

    public DataQuery() {
    }

    public DataQuery(String sql) {
        query = sql;
    }

    public List<FuncParam> getInputParams() {
        List<FuncParam> inputParams = new ArrayList<>();

        for (FuncParam param : params)
            if (param.isInput)
                inputParams.add(param);

        return inputParams;
    }

    public int numInputParams() {
        int cnt = 0;

        if (params != null)
            for (FuncParam param : params)
                if (param.isInput)
                    cnt++;

        return cnt;
    }

    public FuncParam getParam(String name) {
        FuncParam param = null;

        for (FuncParam p : params)
            if (p.name.equals(name))
                param = p;

        return param;
    }

    public void removeParam(FuncParam param) {
        params.remove(param);
    }

    public void addParams(List<FuncParam> params) {
        this.params.addAll(params);
    }
}
