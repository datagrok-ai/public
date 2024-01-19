package grok_connect.log;

import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import org.eclipse.jetty.websocket.api.Session;
import org.slf4j.LoggerFactory;
import java.util.UUID;

public class QueryLoggerImpl implements QueryLogger {
    private final QueryStreamAppender appender;
    private final Logger logger;

    public QueryLoggerImpl(Session session) {
        LoggerContext context = (LoggerContext) LoggerFactory.getILoggerFactory();
        appender = new QueryStreamAppender(session);
        String name = UUID.randomUUID().toString();
        appender.setName(name);
        appender.setContext(context);
        appender.addFilter(new AppenderFilter());
        appender.start();
        logger = context.getLogger(name);
        logger.addAppender(appender);
        logger.setAdditive(true); // set true if root logger should log too
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
    public void writeLog(boolean write) {
        if (write)
            logger.addAppender(appender);
        else
            logger.detachAppender(appender);
    }
}
