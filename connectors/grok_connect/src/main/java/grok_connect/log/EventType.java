package grok_connect.log;

import org.slf4j.Marker;
import org.slf4j.MarkerFactory;

public enum EventType {
    MISC("MISC"),
    SCHEME_INFO_INIT("SCHEME_INFO_INIT"),

    ERROR("ERROR"),
    LOG_SENDING("LOG_SENDING"),
    CHECKSUM_SENDING("CHECKSUM_SENDING"),
    DATAFRAME_TO_BYTEARRAY_CONVERTING("DATAFRAME_TO_BYTEARRAY_CONVERTING"),
    COLUMN_FILLING("COLUMN_FILLING"),
    DATA_SENDING("DATA_SENDING"),
    CONNECTION_RECEIVING("CONNECTION_RECEIVING"),
    QUERY_PARSING("QUERY_PARSING"),
    QUERY_INTERPOLATING("QUERY_INTERPOLATING"),
    STATEMENT_PARAMETERS_REPLACING("STATEMENT_PARAMETERS_REPLACING"),
    STATEMENT_EXECUTING("STATEMENT_EXECUTING");

    private final String name;
    private final Marker marker;

    EventType(String name) {
        this.name = name;
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
        return name;
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
