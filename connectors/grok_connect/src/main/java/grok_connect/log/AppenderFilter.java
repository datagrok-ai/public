package grok_connect.log;

import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.filter.Filter;
import ch.qos.logback.core.spi.FilterReply;
import java.util.List;

public class AppenderFilter extends Filter<ILoggingEvent> {
    private final List<String> levels;

    public AppenderFilter(List<String> levels) {
        this.levels = levels;
    }

    @Override
    public FilterReply decide(ILoggingEvent iLoggingEvent) {
        boolean contains = levels.contains(iLoggingEvent.getLevel().levelStr.toLowerCase());
        return contains ? FilterReply.ACCEPT : FilterReply.DENY;
    }
}
