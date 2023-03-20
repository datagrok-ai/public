package grok_connect;

import serialization.BufferAccessor;
import serialization.DataFrame;
import spark.*;
import java.io.*;
import java.util.*;
import org.joda.time.*;
import javax.servlet.*;
import com.google.gson.*;
import org.apache.log4j.*;

import static spark.Spark.*;
import org.restlet.data.Status;
import javax.ws.rs.core.MediaType;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;
import grok_connect.handlers.QueryHandler;


public class GrokConnect {

    private static ProviderManager providerManager;

    public static ProviderManager getProviderManager() {
        return providerManager;
    }

    private static Logger logger;
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    public static boolean needToReboot = false;

    public static void main(String[] args) {
        int port = 1234;
        String uri = "http://localhost:" + port;

        try {
            BasicConfigurator.configure();
            logger = Logger.getLogger(GrokConnect.class.getName());
            logger.setLevel(Level.INFO);

            logMemory();

            providerManager = new ProviderManager(logger);
            port(port);
            connectorsModule();

            logger.info("grok_connect with Hikari pool");
            logger.info("grok_connect: Running on %s\n" + uri);
            logger.info("grok_connect: Connectors: " + String.join(", ",
            providerManager.getAllProvidersTypes()));
        } catch (Throwable ex) {
            System.out.println("ERROR: " + ex.toString());
            System.out.print("STACK TRACE: ");
            ex.printStackTrace(System.out);
        }
    }

    private static void connectorsModule() {
        webSocket("/query_socket", new QueryHandler());
        // webSocket("/query_table", new QueryHandler(QueryType.tableQuery));

        post("/query", (request, response) -> {
            logMemory();

            BufferAccessor buffer;
            DataQueryRunResult result = new DataQueryRunResult();
            result.log = "";// use builder instead

            FuncCall call = null;
            if (SettingsManager.getInstance().settings != null) {
                try {
                    call = gson.fromJson(request.body(), FuncCall.class);
                    call.log = "";
                    call.setParamValues();
                    call.afterDeserialization();
                    System.out.println(call.func.query);
                    DateTime startTime = DateTime.now();
                    DataProvider provider = providerManager.getByName(call.func.connection.dataSource);
                    DataFrame dataFrame = provider.execute(call);
                    double execTime = (DateTime.now().getMillis() - startTime.getMillis()) / 1000.0;

                    result.blob = dataFrame.toByteArray();
                    result.blobLength = result.blob.length;
                    result.timeStamp = startTime.toString("yyyy-MM-dd hh:mm:ss");
                    result.execTime = execTime;
                    result.columns = dataFrame.columns.size();
                    result.rows = dataFrame.rowCount;
                    // TODO Write to result log there

                    String logString = String.format("%s: Execution time: %f s, Columns/Rows: %d/%d, Blob size: %d bytes\n",
                            result.timeStamp,
                            result.execTime,
                            result.columns,
                            result.rows,
                            result.blobLength);

                    if (call.debugQuery) {
                        result.log += logMemory();
                        result.log += logString;
                    }
                    logger.info(logString);

                    buffer = new BufferAccessor(result.blob);
                    buffer.bufPos = result.blob.length;

                } catch (Throwable ex) {
                    buffer = packException(result,ex);
                    if (ex instanceof OutOfMemoryError)
                        needToReboot = true;
                }
                finally {
                    if (call != null)
                        result.log += call.log;
                }
            }
            else {
                result.errorMessage = NoSettingsException.class.getName();
                buffer = new BufferAccessor();
            }

            try {
                buffer.insertStringHeader(gson.toJson(result));
                buildResponse(response, buffer.toUint8List());
            } catch (Throwable ex) {
                buildExceptionResponse(response, printError(ex));
            }

            return response;
        });

        post("/test", (request, response) -> {
            if (SettingsManager.getInstance().settings == null)
                return NoSettingsException.class.getName();

            DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
            DataProvider provider = providerManager.getByName(connection.dataSource);
            response.type(MediaType.TEXT_PLAIN);
            return provider.testConnection(connection);
        });

        post("/query_table_sql", (request, response) -> {
            logMemory();
            if (SettingsManager.getInstance().settings == null)
                return NoSettingsException.class.getName();

            String result = "";
            try {
                DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                TableQuery query = gson.fromJson(connection.get("queryTable"), TableQuery.class);
                DataProvider provider = providerManager.getByName(connection.dataSource);
                result = provider.queryTableSql(connection, query);
            } catch (Throwable ex) {
                buildExceptionResponse(response, printError(ex));
            }

            logMemory();
            return result;
        });

        post("/schemas", (request, response) -> {
            logMemory();
            BufferAccessor buffer;
            DataQueryRunResult result = new DataQueryRunResult();
            if (SettingsManager.getInstance().settings != null) {
                try {
                    DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                    DataProvider provider = providerManager.getByName(connection.dataSource);
                    DataFrame dataFrame = provider.getSchemas(connection);
                    buffer = packDataFrame(result, dataFrame);
                } catch (Throwable ex) {
                    buffer = packException(result, ex);
                }
            }
            else {
                result.errorMessage = NoSettingsException.class.getName();
                buffer = new BufferAccessor();
            }
            prepareResponse(result, response, buffer);

            logMemory();
            return response;
        });

        post("/schema", (request, response) -> {
            logMemory();

            BufferAccessor buffer;
            DataQueryRunResult result = new DataQueryRunResult();
            if (SettingsManager.getInstance().settings != null) {
                try {
                    DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                    DataProvider provider = providerManager.getByName(connection.dataSource);
                    DataFrame dataFrame = provider.getSchema(connection, connection.get("schema"), connection.get("table"));
                    buffer = packDataFrame(result, dataFrame);
                } catch (Throwable ex) {
                    buffer = packException(result, ex);
                }
            }
            else {
                result.errorMessage = NoSettingsException.class.getName();
                buffer = new BufferAccessor();
            }
            prepareResponse(result, response, buffer);

            logMemory();
            return response;
        });

        get("/conn", (request, response) -> {
            List<DataSource> dataSources = new ArrayList<>();
            for (DataProvider provider : providerManager.providers)
                dataSources.add(provider.descriptor);
            response.type(MediaType.APPLICATION_JSON);
            return gson.toJson(dataSources);
        });

        before((request, response) -> {
            System.out.printf("%s: %s - %s\n", DateTime.now().toString("yyyy-MM-dd hh:mm:ss"),
                    request.requestMethod(), request.pathInfo());
        });

        get("/info", (request, response) -> {
            GrokConnect s = new GrokConnect();
            System.out.println(s.getClass().getPackage().getSpecificationVersion());
            System.out.println(s.getClass().getPackage().getImplementationVersion());
            response.type(MediaType.APPLICATION_JSON);
            return "{\n" +
                    "    \"name\": \"GrokConnect server\",\n" +
                    "    \"version\": \"1.0.3\"\n" +
                    "}";
        });

        get("/log_memory", (request, response) -> logMemory());

        post("/cancel", (request, response) -> {
            FuncCall call = gson.fromJson(request.body(), FuncCall.class);
            providerManager.getQueryMonitor().cancelStatement(call.id);
            providerManager.getQueryMonitor().cancelResultSet(call.id);
            return null;
        });
// how it works, who sends this request?
        post("/set_settings", (request, response) -> {
            System.out.println(request.body());
            Settings settings = gson.fromJson(request.body(), Settings.class);
            SettingsManager.getInstance().setSettings(settings);
            ConnectionPool.getInstance().setTimer();
            return null;
        });
// maybe here would be better to return something in json format, like ReponseEntity(status, body)
        get("/health", (request, response) -> {
            int status;
            String body;
            if (needToReboot) {
                status = 500;
                body = "Grok connect needs a reboot";
            }
            else {
                status = 200;
                body = "OK";
            }

            response.status(status);
            return body;
        });
    }

    private static BufferAccessor packDataFrame(DataQueryRunResult result, DataFrame dataFrame) {
        result.blob = dataFrame.toByteArray();
        result.blobLength = result.blob.length;
        result.columns = dataFrame.columns.size();
        result.rows = dataFrame.rowCount;

        BufferAccessor buffer = new BufferAccessor(result.blob);
        buffer.bufPos = result.blob.length;
        return buffer;
    }

    public static BufferAccessor packException(DataQueryRunResult result, Throwable ex) {
        Map<String, String> exception = printError(ex);
        result.errorMessage = exception.get("errorMessage");
        result.errorStackTrace = exception.get("errorStackTrace");
        return new BufferAccessor();
    }

    private static void prepareResponse(DataQueryRunResult result, Response response, BufferAccessor buffer) {
        try {
            buffer.insertStringHeader(gson.toJson(result));
            buildResponse(response, buffer.toUint8List());
        } catch (Throwable ex) {
            buildExceptionResponse(response, printError(ex));
        }
    }

    private static void buildResponse(Response response, byte[] bytes) throws Throwable {
        response.type(MediaType.APPLICATION_OCTET_STREAM);
        response.raw().setContentLength(bytes.length);
        response.status(Status.SUCCESS_OK.getCode());

        final ServletOutputStream os = response.raw().getOutputStream();
        os.write(bytes);
        os.close();
    }

    private static void buildExceptionResponse(Response response, Map<String, String> exception) {
        response.type(MediaType.TEXT_PLAIN);
        response.body(exception.get("errorMessage") + "\n" + exception.get("errorStackTrace"));
        response.status(Status.SERVER_ERROR_INTERNAL.getCode());
    }

    public static Map<String, String> printError(Throwable ex) {
        String errorMessage = ex.toString();
        StringWriter stackTrace = new StringWriter();
        ex.printStackTrace(new PrintWriter(stackTrace));
        String errorStackTrace = stackTrace.toString();

        System.out.println("ERROR: \n" + errorMessage);
        System.out.print("STACK TRACE: \n" + errorStackTrace);

        return new HashMap<String, String>() {{
            put("errorMessage", errorMessage);
            put("errorStackTrace", errorStackTrace);
        }};
    }

    public static String logMemory() {
        long used = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        long free = Runtime.getRuntime().maxMemory() - used;
        long total = Runtime.getRuntime().maxMemory();

        String str = String.format("Memory: free: %s(%.2f%%), used: %s", free, 100.0 * free/total, used);

        logger.info(str);
        return str;
    }
}
