package grok_connect.log;

import org.slf4j.Marker;
import org.slf4j.MarkerFactory;

public enum EventType {
    MISC("MISC"),
    SCHEME_INFO_INIT("SCHEME_INFO_INIT"),

    ERROR("ERROR"),
    LOG_SEND("LOG_SEND"),
    CHECKSUM_SEND("CHECKSUM_SEND"),
    DATAFRAME_TO_BYTEARRAY_CONVERSION("DATAFRAME_TO_BYTEARRAY_CONVERSION"),
    RESULT_SET_PROCESSING_WITH_DATAFRAME_FILL("RESULT_SET_PROCESSING_WITH_DATAFRAME_FILL"),
    RESULT_SET_PROCESSING_WITHOUT_DATAFRAME_FILL("DRY_RUN"),
    DRY_RUN("DRY_RUN_TOTAL_DURATION"),
    SOCKET_BINARY_DATA_EXCHANGE("SOCKET_BINARY_DATA_EXCHANGE"),
    CONNECTION_RECEIVE("CONNECTION_RECEIVE"),
    QUERY_PARSE("QUERY_PARSE"),
    QUERY_INTERPOLATION("QUERY_INTERPOLATION"),
    STATEMENT_PARAMETERS_REPLACEMENT("STATEMENT_PARAMETERS_REPLACEMENT"),
    STATEMENT_EXECUTION("STATEMENT_EXECUTION");

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
