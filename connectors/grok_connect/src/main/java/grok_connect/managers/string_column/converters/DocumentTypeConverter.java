package grok_connect.managers.string_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Node;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.StringWriter;

public class DocumentTypeConverter implements Converter<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ClobTypeConverter.class);

    @Override
    public String convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        try {
            StringWriter writer = new StringWriter();
            TransformerFactory tf = TransformerFactory.newInstance();
            Transformer transformer = tf.newTransformer();
            transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "no");
            transformer.setOutputProperty(OutputKeys.METHOD, "xml");
            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
            transformer.transform(new DOMSource((Node) value), new StreamResult(writer));
            String converted = writer.toString();
            converted  = converted
                    .replaceAll("&lt;", "<")
                    .replaceAll("&gt;", ">");
            return converted;
        } catch (TransformerException exception) {
            throw new RuntimeException("Something went wrong when converting xml to string", exception);
        }
    }
}
