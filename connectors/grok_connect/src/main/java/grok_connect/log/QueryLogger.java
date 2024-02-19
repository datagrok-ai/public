package grok_connect.log;

import org.slf4j.Logger;

public interface QueryLogger {
    Logger getLogger();

    void closeLogger();

    void writeLog(boolean write);
}
