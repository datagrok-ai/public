package grok_connect.log;

import org.slf4j.Logger;

import java.util.List;

public interface QueryLogger {
    Logger getLogger();

    void closeLogger();

    void writeLog(boolean write);
}
