package grok_connect.log;

import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.AppenderBase;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.StringColumn;
import java.time.Instant;
import java.time.LocalDateTime;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;

public class DataFrameAppender extends AppenderBase<ILoggingEvent> {
    private static final int EVENT_TIME_INDEX = 0;
    private static final int COMPONENT_INDEX = 1;
    private static final int EVENT_LEVEL_INDEX = 2;
    private static final int EVENT_MESSAGE_INDEX = 3;
    private static final String COMPONENT_NAME = "GrokConnect";
    private static final String DATE_FORMAT = "MM/dd/yyyy HH:mm:ss.SSS";
    private final DataFrame logs;

    public DataFrameAppender() {
        DateTimeColumn timeStampColumn = new DateTimeColumn();
        timeStampColumn.name = "Date";
        StringColumn prefixColumn = new StringColumn();
        prefixColumn.name = "Component";
        StringColumn levelColumn = new StringColumn();
        levelColumn.name = "Level";
        StringColumn messageColumn = new StringColumn();
        messageColumn.name = "Message";
        logs = new DataFrame();
        logs.addColumn(timeStampColumn);
        logs.addColumn(prefixColumn);
        logs.addColumn(levelColumn);
        logs.addColumn(messageColumn);
    }

    @Override
    @SuppressWarnings("unchecked")
    protected void append(ILoggingEvent iLoggingEvent) {
        logs.columns.get(EVENT_TIME_INDEX).add(iLoggingEvent.getTimeStamp() * 1000.0);
        logs.columns.get(COMPONENT_INDEX).add(COMPONENT_NAME);
        logs.columns.get(EVENT_LEVEL_INDEX).add(iLoggingEvent.getLevel().toString());
        logs.columns.get(EVENT_MESSAGE_INDEX).add(iLoggingEvent.getFormattedMessage());
        logs.rowCount++;
    }

    public DataFrame getLog() {
        return logs;
    }

    private String getFormattedDate(long timeStamp) {
        return LocalDateTime.ofInstant(Instant.ofEpochMilli(timeStamp), ZoneId.of("UTC"))
                .format(DateTimeFormatter.ofPattern(DATE_FORMAT));
    }
}
