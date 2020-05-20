package utils;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

import com.google.gson.*;
import org.xml.sax.SAXException;
import org.apache.tika.parser.*;
import org.apache.tika.metadata.Metadata;
import org.apache.tika.exception.TikaException;
import org.apache.tika.sax.BodyContentHandler;


public class TikaExtractor {
    public static void main(String[] args)  throws IOException, TikaException, SAXException {
        String input = args[0];
        String output = args[1];

        Gson gson = new GsonBuilder().create();
        File inputFile = new File(input);
        Parser parser = new AutoDetectParser();
        BodyContentHandler handler = new BodyContentHandler(-1);
        Metadata metadata = new Metadata();
        FileInputStream stream = new FileInputStream(inputFile);
        ParseContext context = new ParseContext();
        parser.parse(stream, handler, metadata, context);

        Map<String, Object> map = new HashMap<>();
        for (String name: metadata.names())
            map.put(name, metadata.get(name));
        FileWriter outputFile = new FileWriter(output);
        outputFile.write(gson.toJson(map));
        outputFile.close();
    }
}
