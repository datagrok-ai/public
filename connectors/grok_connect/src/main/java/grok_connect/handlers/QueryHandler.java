package grok_connect.handlers;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.sql.ResultSet;

import org.eclipse.jetty.websocket.api.*;
import org.eclipse.jetty.websocket.api.annotations.*;
import org.joda.time.*;
import com.google.gson.*;
import com.mysql.cj.xdevapi.Result;

import org.apache.log4j.*;

import serialization.*;
import grok_connect.utils.*;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.*;
import grok_connect.providers.JdbcDataProvider;

// TODO add logger

@WebSocket
public class QueryHandler {
    
    @OnWebSocketConnect
    public void connected(Session session) {
        System.out.print("connented socket");
    }

    @OnWebSocketClose
    public void closed(Session session, int statusCode, String reason) {
        
    }

    @OnWebSocketMessage
    public void message(Session session, String message) throws Throwable {
        System.out.print(message);

        Logger logger = Logger.getLogger(GrokConnect.class.getName());
        logger.setLevel(Level.INFO);

        ProviderManager providerManager = new ProviderManager(logger);

        BufferAccessor buffer;

        FuncCall call = null;
        if (SettingsManager.getInstance().settings != null) {
            try {
                Gson gson = new GsonBuilder()
                    .registerTypeAdapter(Property.class, new PropertyAdapter())
                    .create();
                call = gson.fromJson(message, FuncCall.class);
                call.log = "";
                call.setParamValues();
                call.afterDeserialization();
                System.out.println(call.func.query);
 
                System.out.println("sending first df");
                session.getRemote().sendString("sending first df");

                DateTime startTime = DateTime.now();
                JdbcDataProvider provider = providerManager.getByName(call.func.connection.dataSource);

                QueryManager qm = new QueryManager(call, provider);

                qm.getResultSet();

                qm.initScheme();
                DataFrame dataFrame = qm.getSubDF(10);
                while (dataFrame != null && dataFrame.rowCount != 0) {
                    
                    // DataFrame dataFrame = provider.execute(call);
                    double execTime = (DateTime.now().getMillis() - startTime.getMillis()) / 1000.0;
                    
                    System.out.println("rows: ");
                    System.out.println(dataFrame.rowCount);

                    session.getRemote().sendBytes(ByteBuffer.wrap(dataFrame.toByteArray()));

                    dataFrame = qm.getSubDF(10);
                }

                session.close();
                // buffer = new BufferAccessor(result.blob);
                // buffer.bufPos = result.blob.length;

            } catch (Throwable ex) {
                if (ex instanceof OutOfMemoryError)
                    GrokConnect.needToReboot = true;
                else {
                    System.out.println(ex);
                    throw ex;
                }
            }
            finally {
                session.close();
                // if (call != null)
                    // result.log += call.log;
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

}
