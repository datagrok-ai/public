package grok_connect.log;

import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.AppenderBase;
import org.slf4j.Marker;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.util.ArrayList;
import java.util.List;

public class QueryLoggerAppender extends AppenderBase<ILoggingEvent> {
    private static final int EVENT_TIME_INDEX = 0;
    private static final int COMPONENT_INDEX = 1;
    private static final int EVENT_LEVEL_INDEX = 2;
    private static final int EVENT_TYPE_INDEX = 3;
    private static final int EVENT_STAGE_INDEX = 4;
    private static final int EVENT_DF_NUMBER_INDEX = 5;
    private static final int EVENT_SUB_DF_NUMBER_INDEX = 6;
    private static final int EVENT_DESTINATION_INDEX = 7;
    private static final int EVENT_MESSAGE_INDEX = 8;
    private static final int EVENT_DURATION_INDEX = 9;
    private static final String COMPONENT_NAME = "GrokConnect";
    private static final String DESTINATION = "GrokConnect -> DatagrokServer";
    private boolean writeLog = true;
    private final List<ILoggingEvent> logs;

    public QueryLoggerAppender() {
        logs = new ArrayList<>();
    }

    @Override
    protected void append(ILoggingEvent iLoggingEvent) {
        if (writeLog)
            logs.add(iLoggingEvent);
    }

    @SuppressWarnings("unchecked")
    public DataFrame getLog() {
        DataFrame df = getDataFrame();
        for (ILoggingEvent iLoggingEvent: logs) {
            df.columns.get(EVENT_TIME_INDEX).add(iLoggingEvent.getTimeStamp() * 1000.0);
            df.columns.get(COMPONENT_INDEX).add(COMPONENT_NAME);
            df.columns.get(EVENT_LEVEL_INDEX).add(iLoggingEvent.getLevel().toString());
            addMarkerColumns(df, iLoggingEvent.getMarker());
            df.columns.get(EVENT_MESSAGE_INDEX).add(iLoggingEvent.getFormattedMessage());
            df.columns.get(EVENT_DURATION_INDEX).add(0);
            df.rowCount++;
        }
        return df;
    }

    public void setWriteLog(boolean writeLog) {
        this.writeLog = writeLog;
    }

    @SuppressWarnings("unchecked")
    private void addMarkerColumns(DataFrame logs, Marker marker) {
        String[] split = marker.getName().split("\\|");
        String type = split[0];
        String dfNumber = split[1];
        logs.columns.get(EVENT_TYPE_INDEX).add(type);
        logs.columns.get(EVENT_DF_NUMBER_INDEX).add(dfNumber.equals(" ") ? null : Integer.parseInt(dfNumber));
        logs.columns.get(EVENT_SUB_DF_NUMBER_INDEX).add(null);
        logs.columns.get(EVENT_STAGE_INDEX).add(split[2]);
        String destination = type.equals(EventType.CHECKSUM_SEND.toString())
                || type.equals(EventType.DATA_SEND.toString()) || type.equals(EventType.LOG_SEND.toString()) ? DESTINATION : "";
        logs.columns.get(EVENT_DESTINATION_INDEX).add(destination);
    }

    private DataFrame getDataFrame() {
        DateTimeColumn timeStampColumn = new DateTimeColumn();
        timeStampColumn.name = "Date";
        StringColumn prefixColumn = new StringColumn();
        prefixColumn.name = "Component";
        StringColumn levelColumn = new StringColumn();
        levelColumn.name = "Level";
        StringColumn typeColumn = new StringColumn();
        typeColumn.name = "Type";
        StringColumn stageColumn = new StringColumn();
        stageColumn.name = "Stage";
        IntColumn dfNumberColumn = new IntColumn();
        dfNumberColumn.name = "DF #";
        IntColumn subDfNumberColumn = new IntColumn();
        subDfNumberColumn.name = "SubDF #";
        StringColumn destColumn = new StringColumn();
        destColumn.name = "Destination";
        StringColumn messageColumn = new StringColumn();
        messageColumn.name = "Message";
        IntColumn durationColumn = new IntColumn();
        durationColumn.name = "Duration, ms";
        DataFrame logs = new DataFrame();
        logs.addColumn(timeStampColumn);
        logs.addColumn(prefixColumn);
        logs.addColumn(levelColumn);
        logs.addColumn(typeColumn);
        logs.addColumn(stageColumn);
        logs.addColumn(dfNumberColumn);
        logs.addColumn(subDfNumberColumn);
        logs.addColumn(destColumn);
        logs.addColumn(messageColumn);
        logs.addColumn(durationColumn);
        return logs;
    }
}
