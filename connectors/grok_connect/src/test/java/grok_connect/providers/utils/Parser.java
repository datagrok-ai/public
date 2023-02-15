package grok_connect.providers.utils;

public interface Parser {
    Double[] parseDatesToDoubles(String pattern, String... dateToParse);

    Double parseDateToDouble(String pattern, String parseDate);
}
