package grok_connect;

import static spark.Spark.*;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.DataQueryRunResult;
import grok_connect.connectors_info.FuncCall;
import grok_connect.table_query.TableQuery;
import grok_connect.utils.*;
import org.slf4j.LoggerFactory;
import serialization.BufferAccessor;
import serialization.DataFrame;
import java.io.*;
import java.net.HttpURLConnection;
import java.time.Instant;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Properties;
import javax.servlet.ServletOutputStream;
import javax.ws.rs.core.MediaType;
import grok_connect.handlers.QueryHandler;
import spark.Response;

public class GrokConnect {
    private static final int DEFAULT_PORT = 1234;
    private static final String VERSION = "version";
    private static final String NAME = "artifactId";
    private static final String DEFAULT_URI = String.format("http://localhost:%s", DEFAULT_PORT);
    private static final String DEFAULT_LOG_EXCEPTION_MESSAGE = "An exception was thrown";
    private static final Logger PARENT_LOGGER = (Logger) LoggerFactory.getLogger(GrokConnect.class);
    private static final Gson gson = new GsonBuilder()
            .registerTypeAdapter(Property.class, new PropertyAdapter())
            .create();
    private static final String LOG_LEVEL_PREFIX = "LogLevel";
    public static ProviderManager providerManager;
    public static Properties properties;

    public static void main(String[] args) {
        try {
            properties = getInfo();
            setGlobalLogLevel();
            PARENT_LOGGER.info("{} - version: {}", properties.get(NAME), properties.get(VERSION));
            PARENT_LOGGER.info("Grok Connect initializing");
            PARENT_LOGGER.info(getStringLogMemory());
            providerManager = new ProviderManager();
            threadPool(300, 25, -1);
            port(DEFAULT_PORT);
            System.setProperty("com.zaxxer.hikari.housekeeping.periodMs", "60000");
            connectorsModule();
            PARENT_LOGGER.info("grok_connect with Hikari pool");
            PARENT_LOGGER.info("grok_connect: Running on {}", DEFAULT_URI);
            PARENT_LOGGER.info("grok_connect: Connectors: {}", providerManager.getAllProvidersTypes());
        } catch (Throwable ex) {
            PARENT_LOGGER.error(DEFAULT_LOG_EXCEPTION_MESSAGE, ex);
        }
    }

    private static void connectorsModule() {
        webSocket("/query_socket", new QueryHandler());

        before((request, response) -> {
            PARENT_LOGGER.debug("Endpoint {} was called", request.pathInfo());
            PARENT_LOGGER.debug(getStringLogMemory());
        });

        exception(QueryCancelledByUser.class, (exception, request, response) -> {
            FuncCall call = gson.fromJson(request.body(), FuncCall.class);
            QueryMonitor.getInstance().removeResultSet(call.id);
        });

        exception(Exception.class, (exception, request, response) -> {
            if (request.raw().getRequestURI().equals("/test")) {
                StringWriter errors = new StringWriter();
                errors.write("ERROR:\n" + exception.getLocalizedMessage() + "\n\nSTACK TRACE:\n");
                exception.printStackTrace(new PrintWriter(errors));
                response.body(errors.toString());
            } else if (request.raw().getRequestURI().equals("/query_table_sql")) {
                buildExceptionResponse(response, printError(exception));
            } else if (request.raw().getRequestURI().equals("/schema")
                    || request.raw().getRequestURI().equals("/schemas")
                    || request.raw().getRequestURI().equals("/query")) {
                DataQueryRunResult result = new DataQueryRunResult();
                BufferAccessor bufferAccessor = packException(result, exception);
                prepareResponse(result, response, bufferAccessor);
            }
        });

        post("/query", (request, response) -> {
            BufferAccessor buffer;
            DataQueryRunResult result = new DataQueryRunResult();

            FuncCall call = gson.fromJson(request.body(), FuncCall.class);
            call.setParamValues();
            call.afterDeserialization();
            PARENT_LOGGER.debug("Query: {}", call.func.query);
            long startTime = System.currentTimeMillis();
            DataProvider provider = providerManager.getByName(call.func.connection.dataSource);
            DataFrame dataFrame = provider.executeCall(call);
            double execTime = (System.currentTimeMillis() - startTime) / 1000.0;
            result.blob = dataFrame.toByteArray();
            result.blobLength = result.blob.length;
            result.timeStamp = Instant.ofEpochMilli(startTime)
                    .atZone(ZoneId.systemDefault())
                    .toLocalDate().format(DateTimeFormatter.ofPattern("yyyy-MM-dd hh:mm:ss"));
            result.execTime = execTime;
            result.columns = dataFrame.columns.size();
            result.rows = dataFrame.rowCount;
            String logString = String.format("%s: Execution time: %f s, Columns/Rows: %d/%d, Blob size: %d bytes\n",
                    result.timeStamp,
                    result.execTime,
                    result.columns,
                    result.rows,
                    result.blobLength);
            PARENT_LOGGER.debug(logString);
            buffer = new BufferAccessor(result.blob);
            buffer.bufPos = result.blob.length;
            prepareResponse(result, response, buffer);
            return response;
        });

        post("/test", (request, response) -> {
            DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
            DataProvider provider = providerManager.getByName(connection.dataSource);
            response.type(MediaType.TEXT_PLAIN);
            provider.testConnection(connection);
            return DataProvider.CONN_AVAILABLE;
        });

        post("/query_table_sql", (request, response) -> {
            DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
            TableQuery tableQuery = gson.fromJson(connection.get("queryTable"), TableQuery.class);
            DataProvider provider = providerManager.getByName(connection.dataSource);
            return provider.queryTableSql(connection, tableQuery);
        });

        post("/schemas", (request, response) -> {
            DataQueryRunResult result = new DataQueryRunResult();
            DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
            DataProvider provider = providerManager.getByName(connection.dataSource);
            DataFrame dataFrame = provider.getSchemas(connection);
            BufferAccessor buffer = packDataFrame(result, dataFrame);
            prepareResponse(result, response, buffer);
            return response;
        });

        post("/schema", (request, response) -> {
            DataQueryRunResult result = new DataQueryRunResult();
            DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
            DataProvider provider = providerManager.getByName(connection.dataSource);
            DataFrame dataFrame = provider.getSchema(connection, connection.get("schema"), connection.get("table"));
            BufferAccessor buffer = packDataFrame(result, dataFrame);
            prepareResponse(result, response, buffer);
            return response;
        });

        get("/conn", (request, response) -> {
            response.type(MediaType.APPLICATION_JSON);
            return providerManager.getAllDescriptors();
        }, gson::toJson);

        get("/info", (request, response) -> {
            response.type(MediaType.APPLICATION_JSON);
            return properties;
        }, gson::toJson);

        get("/log_memory", (request, response) -> getStringLogMemory());

        post("/cancel", (request, response) -> {
            FuncCall call = gson.fromJson(request.body(), FuncCall.class);
            QueryMonitor.getInstance().cancelStatement(call.id);
            QueryMonitor.getInstance().addCancelledResultSet(call.id);
            response.status(HttpURLConnection.HTTP_OK);
            return response;
        });

        post("/set_settings", (request, response) -> {
            Settings settings = gson.fromJson(request.body(), Settings.class);
            SettingsManager.getInstance().setSettings(settings);
            ConnectionPool.getInstance().setTimer();
            response.status(HttpURLConnection.HTTP_OK);
            return response;
        });

        get("/health", (request, response) -> {
            int status;
            String body;
            if (Runtime.getRuntime().freeMemory() < Runtime.getRuntime().maxMemory() * 0.1) {
                status = HttpURLConnection.HTTP_INTERNAL_ERROR;
                body = "Grok connect needs a reboot";
            } else {
                status = HttpURLConnection.HTTP_OK;
                body = "OK";
            }
            response.status(status);
            return body;
        });

        after(((request, response) -> {
            PARENT_LOGGER.debug("Endpoint {} call was proceeded", request.pathInfo());
            PARENT_LOGGER.debug(getStringLogMemory());
        }));
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


    public static Map<String, String> printError(Throwable ex) {
        String errorMessage = ex.toString();
        StringWriter stackTrace = new StringWriter();
        ex.printStackTrace(new PrintWriter(stackTrace));
        String errorStackTrace = stackTrace.toString();
        return new HashMap<String, String>() {{
            put("errorMessage", errorMessage);
            put("errorStackTrace", errorStackTrace);
        }};
    }

    public static String getStringLogMemory() {
        long used = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
        long free = Runtime.getRuntime().maxMemory() - used;
        long total = Runtime.getRuntime().maxMemory();
        return String.format("Memory: free: %s(%.2f%%), used: %s; active threads: %s", free, 100.0 * free/total, used, activeThreadCount());
    }

    private static void prepareResponse(DataQueryRunResult result, Response response, BufferAccessor buffer) {
        try {
            buffer.insertStringHeader(gson.toJson(result));
            buildResponse(response, buffer.toUint8List());
        } catch (Exception ex) {
            buildExceptionResponse(response, printError(ex));
        }
    }

    private static void buildResponse(Response response, byte[] bytes) throws IOException {
        response.type(MediaType.APPLICATION_OCTET_STREAM);
        response.raw().setContentLength(bytes.length);
        response.status(HttpURLConnection.HTTP_OK);
        ServletOutputStream os = response.raw().getOutputStream();
        os.write(bytes);
        os.close();
    }

    private static void buildExceptionResponse(Response response, Map<String, String> exception) {
        response.type(MediaType.TEXT_PLAIN);
        response.body(exception.get("errorMessage") + "\n" + exception.get("errorStackTrace"));
        response.status(HttpURLConnection.HTTP_INTERNAL_ERROR);
    }

    public static Properties getInfo() {
        try {
            InputStream resourceAsStream = Thread.currentThread()
                    .getContextClassLoader()
                    .getResourceAsStream("app.properties");
            Properties properties = new Properties();
            if (resourceAsStream != null) {
                properties.load(new InputStreamReader(resourceAsStream));
            }
            return properties;
        } catch (IOException e) {
            throw new RuntimeException("Something went wrong when getting info", e);
        }
    }

    private static void setGlobalLogLevel() {
        Optional<String> level = System.getProperties().stringPropertyNames().stream()
                .filter(name -> name.startsWith(LOG_LEVEL_PREFIX))
                .findFirst();
        if (level.isPresent()) {
            String property = System.getProperty(level.get());
            PARENT_LOGGER.setLevel(Level.toLevel(property));
        }
    }
}
