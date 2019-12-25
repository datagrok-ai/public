package grok_connect.connectors_info;

import com.google.api.client.util.*;


public class DataRequest
{
    public String userId;
    public String sessionId;
    public DateTime requestedOn;
    public DateTime completedOn;

    public DataSource source;
    public DataConnection connection;
    public DataQuery query;
}
