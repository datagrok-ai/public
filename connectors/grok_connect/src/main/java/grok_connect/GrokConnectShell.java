package grok_connect;

import java.io.*;
import java.nio.file.*;
import java.nio.charset.*;
import org.joda.time.*;
import com.google.gson.*;
import java.sql.SQLException;
import org.apache.commons.cli.*;

import serialization.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class GrokConnectShell {
    public static void main(String[] args) throws ParseException, IOException, SQLException,
            java.text.ParseException, ClassNotFoundException {

        Options options = new Options();

        Option input = new Option("q", "query", true, "Query JSON file path");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "Output CSV file path");
        output.setRequired(false);
        options.addOption(output);

        Option help = new Option("h", "help", false, "Help");
        help.setRequired(false);
        options.addOption(output);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = parser.parse(options, args);

        if (cmd.getOptionValue("help") != null) {
            formatter.printHelp("utility-name", options);
            System.exit(0);
        }

        Gson gson = new GsonBuilder()
                .registerTypeAdapter(Property.class, new PropertyAdapter())
                .create();

        FuncCall call = gson.fromJson(new String(Files.readAllBytes(Paths.get(cmd.getOptionValue("query"))), StandardCharsets.UTF_8), FuncCall.class);
        call.setParamValues();
        DateTime startTime = DateTime.now();
        DataProvider provider = DataProvider.getByName(call.func.connection.dataSource);
        DataFrame table = provider.execute(call);
        double execTime = (DateTime.now().getMillis() - startTime.getMillis()) / 1000.0;

        System.out.printf("\n%s: Execution time: %f s, Columns/Rows: %d/%d\n\n",
                startTime.toString("yyyy-MM-dd hh:mm:ss"),
                execTime,
                table.columns.size(),
                table.rowCount);

        String csv = table.toCsv();
        String csvPath = cmd.getOptionValue("output");

        if (csvPath != null)
            Files.write(Paths.get(csvPath), csv.getBytes(StandardCharsets.UTF_8));
        else
            System.out.println(csv);
    }
}
