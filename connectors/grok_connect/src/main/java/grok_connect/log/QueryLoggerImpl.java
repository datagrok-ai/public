package grok_connect.log;

import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import com.google.gson.Gson;
import org.slf4j.LoggerFactory;
import java.util.List;
import java.util.UUID;

public class QueryLoggerImpl implements QueryLogger {
    private final DataFrameAppender appender;
    private final Logger logger;

    public QueryLoggerImpl(List<String> logLevels) {
        LoggerContext context = (LoggerContext) LoggerFactory.getILoggerFactory();
        appender = new DataFrameAppender();
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
    public String dumpLogMessages() {
        return appender.getLog().toCsv();
    }
}
