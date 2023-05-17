package grok_connect.log;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.AppenderBase;
import org.apache.commons.collections4.MultiValuedMap;
import org.apache.commons.collections4.multimap.ArrayListValuedHashMap;
import java.time.Instant;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.Collection;
import java.util.Map;
import java.util.TimeZone;

public class MapAppender extends AppenderBase<ILoggingEvent> {
    private static final String PATTERN_LAYOUT = "%s (GrokConnect) [%s] - %s";
    private static final String TIME_FORMAT = "HH:mm:ss.SSS";
    private final MultiValuedMap<Long, String> log;

    public MapAppender() {
        this.log = new ArrayListValuedHashMap<>();
    }

    @Override
    protected void append(ILoggingEvent iLoggingEvent) {
        log.put(iLoggingEvent.getTimeStamp(), formatMessage(iLoggingEvent));
    }

    public Map<Long, Collection<String>> getLog() {
        return log.asMap();
    }

    private String formatMessage(ILoggingEvent iLoggingEvent) {
        String timeStamp = LocalDateTime
                .ofInstant(Instant.ofEpochMilli(iLoggingEvent.getTimeStamp()), TimeZone.getDefault().toZoneId())
                .format(DateTimeFormatter.ofPattern(TIME_FORMAT));
        Level level = iLoggingEvent.getLevel();
        String message = iLoggingEvent.getFormattedMessage();
        return String.format(PATTERN_LAYOUT, timeStamp, level, message);
    }
}
