package grok_connect.handlers;

import org.eclipse.jetty.websocket.api.Session;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import java.nio.ByteBuffer;
import java.util.Map;

import org.apache.log4j.Logger;

import grok_connect.GrokConnect;
import grok_connect.utils.NoSettingsException;
import grok_connect.utils.QueryChunkNotSent;
import grok_connect.utils.QueryManager;
import grok_connect.utils.QueryType;
import grok_connect.utils.SettingsManager;
import serialization.DataFrame;

public class SessionHandler {
    static int rowsPerChunk = 10000;
    
    Session session;
    QueryType queryType;
    DataFrame dataFrame;

    Boolean firstTry = true;

    Boolean oneDfSent = false;

    byte[] bytes;

    QueryManager qm;

    int queryMessages = 0;
    static String okResponse = "DATAFRAME PART OK";
    static String endMessage = "the end";
    static String sizeRecievedMessage = "DATAFRAME PART SIZE RECEIVED";

    SessionHandler(Session session, QueryType queryType) {
        this.session = session;
        this.queryType = queryType;
    }

    public void onMessage(String message) throws Throwable {
        if (SettingsManager.getInstance().settings != null) {
            try {
                if (message.startsWith("QUERY")) {
                    System.out.print("query");
                    if (queryMessages++ > 0 ) {
                        System.out.println("ALERT SECOND QUERY TO SAME SOCKET WTF??" + queryMessages);
                        System.out.println(queryMessages);
                    }
                    message = message.substring(6);

                    Logger logger = Logger.getLogger(GrokConnect.class.getName());
                    qm = new QueryManager(message, queryType);
                    qm.getResultSet();
                    qm.initScheme();

                    dataFrame = qm.getSubDF(rowsPerChunk);
                } else if (message.startsWith(sizeRecievedMessage)) {
                    System.out.print("sending bytes");
                    session.getRemote().sendBytes(ByteBuffer.wrap(bytes));
                    return;
                } else {
                    if (!message.equals(okResponse)) {
                        System.out.print("not ok response?" + message);
                        if (!firstTry)
                            throw new QueryChunkNotSent();
                        else {
                            firstTry = false;
                        }
                    }
                    else {
                        System.out.print("ok response?" + message);
                        firstTry = true;
                        oneDfSent = true;
                        dataFrame = qm.getSubDF(rowsPerChunk);
                    }
                }
                System.out.println("sending df info");
                if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {
                    
                    System.out.println("rows: ");
                    System.out.println(dataFrame.rowCount);
                    bytes = dataFrame.toByteArray();

                    session.getRemote().sendString(checksumMessage(bytes.length));
                    return;
                } else {
                    System.out.println("df empty, end");
                    session.getRemote().sendString(endMessage);
                    session.close();
                }

            } catch (Throwable ex) {
                if (ex instanceof OutOfMemoryError)
                    GrokConnect.needToReboot = true;
                session.getRemote().sendString(socketErrorMessage(ex));
                session.close();
            }

        }
        else {
            // result.errorMessage = NoSettingsException.class.getName();
            System.out.println(NoSettingsException.class.getName());
            session.getRemote().sendString(NoSettingsException.class.getName());
            session.close();
            // buffer = new BufferAccessor();
        }

    }

    static String checksumMessage(int i) {
        return String.format("DATAFRAME PART SIZE: %d", i);
    }

    static String socketErrorMessage(Throwable th) {
        Gson gson = new GsonBuilder().create();
        return "Error: ".concat(gson.toJson(GrokConnect.printError(th), new TypeToken<Map<String, String>>() { }.getType()));
    }
}
