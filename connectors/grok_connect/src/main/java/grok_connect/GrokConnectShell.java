package grok_connect;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.FuncCall;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.Property;
import grok_connect.utils.PropertyAdapter;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryCancelledByUser;
import grok_connect.utils.SettingsManager;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.time.Instant;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.DataFrame;

public class GrokConnectShell {
    private static final Logger LOGGER = LoggerFactory.getLogger(GrokConnectShell.class);

    public static void main(String[] args) throws ParseException, IOException, SQLException,
            java.text.ParseException, ClassNotFoundException, QueryCancelledByUser, GrokConnectException {
        LOGGER.debug("Received args: {}", Arrays.toString(args));
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
        long startTime = System.currentTimeMillis();
        DataProvider provider = new ProviderManager().getByName(call.func.connection.dataSource);
        provider.outputCsv = cmd.getOptionValue("original_csv");
        DataFrame table = provider.execute(call);
        double execTime = (System.currentTimeMillis() - startTime) / 1000.0;

        LOGGER.debug("{}: Execution time: {} s, Columns/Rows: {}/{}\n",
                Instant.ofEpochMilli(startTime)
                        .atZone(ZoneId.systemDefault())
                        .toLocalDate().format(DateTimeFormatter.ofPattern("yyyy-MM-dd hh:mm:ss")),
                execTime,
                table.getColumnCount(),
                table.rowCount);

        String csv = table.toCsv();
        String csvPath = cmd.getOptionValue("output");

        if (csvPath != null)
            Files.write(Paths.get(csvPath), csv.getBytes(StandardCharsets.UTF_8));
        else
            LOGGER.debug(csv);

        String binaryPath = cmd.getOptionValue("binary_output");

        if (binaryPath != null) {
            try (FileOutputStream outputStream = new FileOutputStream(binaryPath)) {
                outputStream.write(table.toByteArray());
            }
        }

    }
}
