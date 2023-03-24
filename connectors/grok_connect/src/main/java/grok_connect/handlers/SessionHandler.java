package grok_connect.handlers;

import org.eclipse.jetty.websocket.api.Session;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import java.nio.ByteBuffer;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Map;
import java.util.concurrent.*;

import grok_connect.GrokConnect;
import grok_connect.utils.*;
import serialization.DataFrame;

public class SessionHandler {
    static int rowsPerChunk = 10000;
    Session session;
    ExecutorService threadpool = Executors.newCachedThreadPool();
    Future<DataFrame> fdf;
    DataFrame dataFrame;

    Boolean firstTry = true;

    Boolean oneDfSent = false;

    byte[] bytes;

    QueryManager qm;

    int queryMessages = 0;
    static String okResponse = "DATAFRAME PART OK";
    static String endMessage = "EOF";
    static String sizeRecievedMessage = "DATAFRAME PART SIZE RECEIVED";
    static String logRecievedMessage = "LOG RECEIVED";

    SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");
    SessionHandler(Session session) {
        this.session = session;
    }

    public void onError(Throwable err) throws Throwable {
        session.getRemote().sendString(socketErrorMessage(err));
        session.close();
        if (qm != null)
            qm.closeConnection();
    }

    public void onMessage(String message) throws Throwable {
        if (SettingsManager.getInstance().settings != null) {
            try {
                if (message.startsWith("QUERY")) {
                    message = message.substring(6);
                    qm = new QueryManager(message);
                    if (qm.query.debugQuery) {
                        qm.query.log += "Getting resultSet, " + sdf.format(new Date());
                        System.out.println("Getting resultSet, " + sdf.format(new Date()));
                    }
                    qm.getResultSet();
                    if (qm.query.debugQuery) {
                        qm.query.log += "Got resultSet, " + sdf.format(new Date());
                        System.out.println("Got resultSet, " + sdf.format(new Date()));
                    }
                    if (qm.resultSet != null) {
                        if (qm.query.debugQuery) {
                            qm.query.log += "initting scheme, " + sdf.format(new Date());
                            System.out.println("initting scheme, " + sdf.format(new Date()));
                        }
                        qm.initScheme();
                        if (qm.query.debugQuery) {
                            qm.query.log += "inited scheme, " + sdf.format(new Date());
                            System.out.println("inited scheme, " + sdf.format(new Date()));
                        }
                        dataFrame = qm.getSubDF(100);
                    } else {
                        dataFrame = new DataFrame(); 
                    }

                } else if (message.startsWith(logRecievedMessage)) {
                    session.getRemote().sendString(endMessage);
                    session.close();
                    return;
                } else if (message.startsWith(sizeRecievedMessage)) {
                    fdf = threadpool.submit(() -> qm.getSubDF(rowsPerChunk));
                    session.getRemote().sendBytes(ByteBuffer.wrap(bytes));
                    return;
                } else {
                    if (!message.equals(okResponse)) {
                        if (!firstTry)
                            throw new QueryChunkNotSent();
                        else {
                            firstTry = false;
                        }
                    }
                    else {
                        firstTry = true;
                        oneDfSent = true;
                        if (qm.resultSet != null)
                            dataFrame = fdf.get();
                    }
                }
                if (dataFrame != null && (dataFrame.rowCount != 0 || !oneDfSent)) {   
                    bytes = dataFrame.toByteArray();
                    session.getRemote().sendString(checksumMessage(bytes.length));
                    return;
                } else {
                    qm.closeConnection();

                    if (qm.query.debugQuery) {
                        session.getRemote().sendString(socketLogMessage(qm.query.log));
                        return;
                    }
                    
                    session.getRemote().sendString(endMessage);
                    session.close();
                }

            } catch (Throwable ex) {
                if (ex instanceof OutOfMemoryError)
                    GrokConnect.needToReboot = true;
                if (ex.getCause() instanceof QueryCancelledByUser)
                    ex = ex.getCause();
                session.getRemote().sendString(socketErrorMessage(ex));
                session.close();
                qm.closeConnection();
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
        return "ERROR: ".concat(gson.toJson(GrokConnect.printError(th), new TypeToken<Map<String, String>>() { }.getType()));
    }

    static String socketLogMessage(String s) {
        return "LOG: " + s;
    }
}
