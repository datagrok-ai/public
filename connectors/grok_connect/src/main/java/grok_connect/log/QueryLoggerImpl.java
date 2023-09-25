package grok_connect.log;

import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import org.slf4j.LoggerFactory;
import serialization.DataFrame;
import java.util.List;
import java.util.UUID;

public class QueryLoggerImpl implements QueryLogger<DataFrame> {
    private final QueryLoggerAppender appender;
    private final Logger logger;

    public QueryLoggerImpl(List<String> logLevels) {
        LoggerContext context = (LoggerContext) LoggerFactory.getILoggerFactory();
        appender = new QueryLoggerAppender();
        UUID uuid = UUID.randomUUID();
        appender.setName(uuid.toString());
        appender.setContext(context);
        appender.addFilter(new AppenderFilter(logLevels));
        appender.start();
        logger = context.getLogger(uuid.toString());
        logger.addAppender(appender);
        logger.setAdditive(true);
    }

    @Override
    public Logger getLogger() {
        return logger;
    }

    @Override
    public void closeLogger() {
        logger.detachAndStopAllAppenders();
    }

    @Override
    public DataFrame dumpLogMessages() {
        return appender.getLog();
    }

    @Override
    public void writeLog(boolean write) {
        appender.setWriteLog(write);
    }
}
