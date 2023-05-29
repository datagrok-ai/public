package grok_connect.log;

import org.slf4j.Marker;
import org.slf4j.MarkerFactory;

public enum EventType {
    MISC("MISC"),
    CONNECTION_RECEIVING("CONNECTION_RECEIVING"),
    QUERY_PARSING("QUERY_PARSING"),
    QUERY_INTERPOLATING("QUERY_INTERPOLATING"),
    STATEMENT_PARAMETERS_REPLACING("STATEMENT_PARAMETERS_REPLACING"),
    STATEMENT_EXECUTING("STATEMENT_EXECUTING"),
    CHUNK_SENDING("CHUNK_SENDING"),
    DATAFRAME_FILLING("DATAFRAME_FILLING"),
    DATAFRAME_TO_BYTEARRAY_CONVERTING("DATAFRAME_TO_BYTEARRAY_CONVERTING"),
    CHECKSUM_SENDING("CHECKSUM_SENDING"),
    RESULT_SET_INIT("RESULT_SET_INIT"),
    RESULT_SET_PROCESSING("RESULT_SET_PROCESSING"),
    SCHEME_INFO_INIT("SCHEME_INFO_INIT"),
    ERROR("ERROR"),
    FETCH_SIZE_CHANGING("FETCH_SIZE_CHANGING"),
    SOCKET_MESSAGE_PROCESSING("GROK_CONNECT_SOCKET_MESSAGE_PROCESSING"),
    LOG_PROCESSING("LOG_PROCESSING");

    private final String name;
    private final Marker marker;

    EventType(String name) {
        this.name = name;
        this.marker = MarkerFactory.getMarker(name);
    }

    public Marker getMarker() {
        return marker;
    }

    public Marker getMarkerNumbered(Integer stageNumber) {
        return MarkerFactory.getMarker(String.format("%s #%s", name, stageNumber));
    }
}
