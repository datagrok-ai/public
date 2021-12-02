package grok_connect;

import java.io.*;
import java.nio.file.*;
import java.nio.charset.*;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.joda.time.*;
import com.google.gson.*;
import java.sql.SQLException;
import org.apache.commons.cli.*;

import serialization.*;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class GrokConnectShell {
    public static void main(String[] args) throws ParseException, IOException, SQLException,
            java.text.ParseException, ClassNotFoundException, QueryCancelledByUser, GrokConnectException {

        Options options = new Options();

        Option input = new Option("q", "query", true, "Query JSON file path");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "Output CSV file path");
        output.setRequired(false);
        options.addOption(output);

        Option binaryOutput = new Option("bo", "binary_output", true, "Output binary file path");
        binaryOutput.setRequired(false);
        options.addOption(binaryOutput);

        Option originalCsv = new Option("oc", "original_csv", true, "Output original CSV file path");
        originalCsv.setRequired(false);
        options.addOption(originalCsv);

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
        SettingsManager.getInstance().initSettingsWithDefaults();
        DateTime startTime = DateTime.now();
        BasicConfigurator.configure();
        Logger logger = Logger.getLogger(GrokConnect.class.getName());
        logger.setLevel(Level.INFO);
        DataProvider provider = new ProviderManager(logger).getByName(call.func.connection.dataSource);
        provider.outputCsv = cmd.getOptionValue("original_csv");
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

        String binaryPath = cmd.getOptionValue("binary_output");

        if (binaryPath != null) {
            try (FileOutputStream outputStream = new FileOutputStream(new File(binaryPath))) {
                outputStream.write(table.toByteArray());
            }
        }

    }
}
