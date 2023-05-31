package grok_connect.log;

import org.slf4j.Marker;
import org.slf4j.MarkerFactory;

public enum EventType {
    MISC("MISC", null),

    SCHEME_INFO_INIT("SCHEME_INFO_INIT", null),

    ERROR("ERROR", null),

    LOG_PROCESSING("LOG_PROCESSING", null),

    DATAFRAME_PROCESSING("DATAFRAME_PROCESSING", null),
    CHECKSUM_SENDING("CHECKSUM_SENDING", DATAFRAME_PROCESSING),
    DATAFRAME_TO_BYTEARRAY_CONVERTING("DATAFRAME_TO_BYTEARRAY_CONVERTING", DATAFRAME_PROCESSING),
    RESULT_SET_PROCESSING("RESULT_SET_PROCESSING", DATAFRAME_PROCESSING),
    COLUMN_FILLING("COLUMN_FILLING", RESULT_SET_PROCESSING),
    DATAFRAME_FILLING("DATAFRAME_FILLING", RESULT_SET_PROCESSING),
    CHUNK_SENDING("CHUNK_SENDING", DATAFRAME_PROCESSING),

    RESULT_SET_INIT("RESULT_SET_INIT", null),
    CONNECTION_RECEIVING("CONNECTION_RECEIVING", RESULT_SET_INIT),
    QUERY_PARSING("QUERY_PARSING", RESULT_SET_INIT),
    QUERY_INTERPOLATING("QUERY_INTERPOLATING", RESULT_SET_INIT),
    STATEMENT_PARAMETERS_REPLACING("STATEMENT_PARAMETERS_REPLACING", RESULT_SET_INIT),
    STATEMENT_EXECUTING("STATEMENT_EXECUTING", RESULT_SET_INIT);

    private final String name;
    private final Marker marker;
    private final EventType parent;

    EventType(String name, EventType parent) {
        this.name = name;
        this.parent = parent;
        this.marker = MarkerFactory.getMarker(String.format("%s| | ", this));
    }

    public Marker getMarker() {
        return marker;
    }

    public Marker getMarker(Stage stage) {
        String formattedName = String.format("%s| |%s", this, stage);
        return MarkerFactory.getMarker(formattedName);
    }

    public Marker getMarker(Integer dfNumber) {
        String formattedName = String.format("%s|%s| ", this, dfNumber);
        return MarkerFactory.getMarker(formattedName);
    }

    public Marker getMarker(Integer dfNumber, Stage stage) {
        String formattedName = String.format("%s|%s|%s", this, dfNumber, stage);
        return MarkerFactory.getMarker(formattedName);
    }

    public String getName() {
        return name;
    }

    @Override
    public String toString() {
        String parentStr = parent == null ? "" : String.format("%s.", parent);
        return String.format("%s%s", parentStr, name);
    }

    public enum Stage {
        START("START"),
        INTERMEDIATE("INTERMEDIATE"),
        END("END");

        private final String name;

        Stage(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }
    }
}
