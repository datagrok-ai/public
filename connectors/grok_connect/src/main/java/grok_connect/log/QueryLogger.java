package grok_connect.log;

import org.slf4j.Logger;
import serialization.DataFrame;

public interface QueryLogger {
    Logger getLogger();

    void closeLogger();

    DataFrame dumpLogMessages();
}
