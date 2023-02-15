package grok_connect.providers.utils;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

/**
 * Class for parsing dates from mock data into Doubles[] (required for building DateTime Column for DataFrame)
 */
public class DateParser implements Parser {
    public Double[] parseDatesToDoubles(String pattern, String... dateToParse) {
        return Arrays.stream(dateToParse)
                .map(date -> parseDateToDouble(pattern, date))
                .toArray(Double[]::new);
    }

    public Double parseDateToDouble(String pattern, String parseDate) {
        Date date;
        SimpleDateFormat formatter = new SimpleDateFormat(pattern);
        try {
            date = formatter.parse(parseDate);
        } catch (ParseException e) {
            throw new RuntimeException("Something went wrong when parsing " + parseDate, e);
        }
        return date == null ? null : date.getTime() * 1000.0; // taken from execute method of JdbcProvider
    }
}
