package grok_connect;

import spark.*;
import java.io.*;
import java.util.*;
import org.joda.time.*;
import javax.servlet.*;
import com.google.gson.*;
import org.apache.log4j.*;
import static spark.Spark.*;
import serialization.*;
import org.restlet.data.Status;
import javax.ws.rs.core.MediaType;
import grok_connect.utils.*;
import grok_connect.table_query.*;
import grok_connect.connectors_info.*;


public class GrokConnect {
    public static void main(String[] args) {
        int port = 1234;
        String uri = "http://localhost:" + port;

        try {
            BasicConfigurator.configure();
            Logger logger = Logger.getLogger("org.eclipse.jetty");
            logger.setLevel(Level.ERROR);

            port(port);
            connectorsModule();

            System.out.printf("grok_connect: Running on %s\n", uri);
            System.out.println("grok_connect: Connectors: " + String.join(", ",
                    DataProvider.getAllProvidersTypes()));
        } catch (Throwable ex) {
            System.out.println("ERROR: " + ex.toString());
            System.out.print("STACK TRACE: ");
            ex.printStackTrace(System.out);
        }
    }

    private static void connectorsModule() {
        Gson gson = new GsonBuilder()
                .registerTypeAdapter(Property.class, new PropertyAdapter())
                .create();

        post("/query", (request, response) -> {

            BufferAccessor buffer;
            DataQueryRunResult result = new DataQueryRunResult();

            try {
                FuncCall call = gson.fromJson(request.body(), FuncCall.class);
                call.setParamValues();
                System.out.println(call.func.query);
                DateTime startTime = DateTime.now();
                DataProvider provider = DataProvider.getByName(call.func.connection.dataSource);
                DataFrame dataFrame = provider.execute(call);
                double execTime = (DateTime.now().getMillis() - startTime.getMillis()) / 1000.0;

                result.blob = dataFrame.toByteArray();
                result.blobLength = result.blob.length;
                result.timeStamp = startTime.toString("yyyy-MM-dd hh:mm:ss");
                result.execTime = execTime;
                result.columns = dataFrame.columns.size();
                result.rows = dataFrame.rowCount;
                // TODO Write to result log there

                System.out.printf("%s: Execution time: %f s, Columns/Rows: %d/%d, Blob size: %d bytes\n",
                        result.timeStamp,
                        result.execTime,
                        result.columns,
                        result.rows,
                        result.blobLength);

                buffer = new BufferAccessor(result.blob);
                buffer.bufPos = result.blob.length;

            } catch (Throwable ex) {
                Map<String, String> exception = printError(ex);
                result.errorMessage = exception.get("errorMessage");
                result.errorStackTrace = exception.get("errorStackTrace");

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
            DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
            DataProvider provider = DataProvider.getByName(connection.dataSource);
            response.type(MediaType.TEXT_PLAIN);
            return provider.testConnection(connection);
        });

        post("/query_table", (request, response) -> {
            try {
                DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                TableQuery query = gson.fromJson((String)connection.parameters.get("queryTable"), TableQuery.class);
                DataProvider provider = DataProvider.getByName(connection.dataSource);
                DataFrame result = provider.queryTable(connection, query);
                buildResponse(response, result.toByteArray());
            } catch (Throwable ex) {
                buildExceptionResponse(response, printError(ex));
            }

            return response;
        });

        post("/query_table_sql", (request, response) -> {
            String result = "";
            try {
                DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                TableQuery query = gson.fromJson(connection.get("queryTable"), TableQuery.class);
                DataProvider provider = DataProvider.getByName(connection.dataSource);
                result = provider.queryTableSql(connection, query);
            } catch (Throwable ex) {
                buildExceptionResponse(response, printError(ex));
            }

            return result;
        });

        post("/schemas", (request, response) -> {
            try {
                DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                DataProvider provider = DataProvider.getByName(connection.dataSource);
                DataFrame result = provider.getSchemas(connection);
                buildResponse(response, result.toByteArray());
            } catch (Throwable ex) {
                buildExceptionResponse(response, printError(ex));
            }

            return response;
        });

        post("/schema", (request, response) -> {
            try {
                DataConnection connection = gson.fromJson(request.body(), DataConnection.class);
                DataProvider provider = DataProvider.getByName(connection.dataSource);
                DataFrame result = provider.getSchema(connection, connection.get("schema"), connection.get("table"));
                buildResponse(response, result.toByteArray());
            } catch (Throwable ex) {
                buildExceptionResponse(response, printError(ex));
            }

            return response;
        });

        get("/conn", (request, response) -> {
            List<DataSource> dataSources = new ArrayList<>();
            for (DataProvider provider : DataProvider.Providers)
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

    private static Map<String, String> printError(Throwable ex) {
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
}
