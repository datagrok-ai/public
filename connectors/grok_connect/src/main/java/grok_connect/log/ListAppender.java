package grok_connect.log;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.AppenderBase;
import java.util.ArrayList;
import java.util.List;

public class ListAppender extends AppenderBase<ILoggingEvent> {
    private static final String PATTERN_LAYOUT = "%s GrokConnect %s %s";
    private final List<String> log;

    public ListAppender() {
        this.log = new ArrayList<>();
    }

    @Override
    protected void append(ILoggingEvent iLoggingEvent) {
        log.add(formatMessage(iLoggingEvent));
    }

    public List<String> getLog() {
        return log;
    }

    private String formatMessage(ILoggingEvent iLoggingEvent) {
        Level level = iLoggingEvent.getLevel();
        String message = iLoggingEvent.getFormattedMessage();
        return String.format(PATTERN_LAYOUT, iLoggingEvent.getTimeStamp(), level, message);
    }
}
