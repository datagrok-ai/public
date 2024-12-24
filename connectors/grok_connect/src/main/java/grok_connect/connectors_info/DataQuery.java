package grok_connect.connectors_info;

import java.util.*;

import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.text.*;


public class DataQuery {
    @SerializedName("#type")
    public String type;
    public String id;
    public String name;
    public String query;
    public String connectionId;
    public DataConnection connection;
    public List<FuncParam> params;
    public Map<String, Object> options;
    public Map<String, Object> aux;

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

    public int inputParamsCount() {
        int cnt = 0;

        if (params != null)
            for (FuncParam param : params)
                if (param.isInput)
                    cnt++;

        return cnt;
    }

    public FuncParam getParam(String name) {
         for (FuncParam p : params)
            if (p.name.equals(name))
                return p;
         return null;
    }

    public boolean existsParam(String name) {
        for (FuncParam p : params)
            if (p.name.equals(name))
                return true;
        return false;
    }

    public void removeParam(FuncParam param) {
        params.remove(param);
    }

    public void removeParams(List<FuncParam> params) {
        this.params.removeAll(params);
    }

    public void addParams(List<FuncParam> params) {
        this.params.addAll(params);
    }
}
