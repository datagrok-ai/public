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
    // §6.2 post-hoc raw-write detection (connector-writes WO-B13): true when a no-result-set
    // statement produced an update count while allowRawWrites auditing was requested.
    public boolean rawWriteDetected;

    public transient byte[] blob;

    public DataQueryRunResult() {
        type = "DataQueryRunResult";
    }
}
