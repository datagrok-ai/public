package grok_connect.connectors_info;

import com.google.gson.annotations.SerializedName;


public class DataQueryRunResult {
    @SerializedName("#type")
    public String type;

    public String timeStamp;
    public double execTime;
    public int columns;
    public int rows;
    public int blobLength;
    public String errorMessage;
    public String errorStackTrace;
    public String log;

    public transient byte[] blob;

    public DataQueryRunResult() {
        type = "DataQueryRunResult";
    }
}
