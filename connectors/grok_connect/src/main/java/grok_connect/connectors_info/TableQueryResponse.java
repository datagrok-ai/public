package grok_connect.connectors_info;

import java.util.List;

public class TableQueryResponse {
    public String query;
    public List<FuncParam> params;

    public TableQueryResponse(String query, List<FuncParam> params) {
        this.query = query;
        this.params = params;
    }
}
