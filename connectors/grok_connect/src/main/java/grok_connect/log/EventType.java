package grok_connect.log;

import org.slf4j.Marker;
import org.slf4j.MarkerFactory;

public enum EventType {
    MISC("MISC"),
    ERROR("ERROR"),
    CHECKSUM_SEND("CHECKSUM SEND"),
    DATAFRAME_TO_BYTEARRAY_CONVERSION("DATAFRAME TO BYTEARRAY CONVERSION"),
    RESULT_SET_PROCESSING_WITH_DATAFRAME_FILL("RESULT SET PROCESSING WITH DATAFRAME FILL"),
    RESULT_SET_PROCESSING_WITHOUT_DATAFRAME_FILL("DRY RUN"),
    DRY_RUN("DRY RUN TOTAL DURATION"),
    SOCKET_BINARY_DATA_EXCHANGE("SOCKET BINARY DATA EXCHANGE"),
    CONNECTION_RECEIVE("CONNECTION RECEIVE"),
    QUERY_PARSE("SQL QUERY PATTERNS PARSE"),
    QUERY_INTERPOLATION("SQL QUERY INTERPOLATION"),
    STATEMENT_PARAMETERS_REPLACEMENT("STATEMENT PARAMETERS REPLACEMENT"),
    STATEMENT_EXECUTION("STATEMENT EXECUTION");

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
