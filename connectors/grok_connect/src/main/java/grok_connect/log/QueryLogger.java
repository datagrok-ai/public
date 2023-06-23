package grok_connect.log;

import org.slf4j.Logger;

public interface QueryLogger<T> {
    Logger getLogger();

    void closeLogger();

    T dumpLogMessages();

    void writeLog(boolean write);
}
