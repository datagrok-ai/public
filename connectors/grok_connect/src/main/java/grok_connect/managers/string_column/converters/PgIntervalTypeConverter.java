package grok_connect.managers.string_column.converters;

import grok_connect.managers.Converter;
import org.postgresql.util.PGInterval;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;

// pgjdbc 42.7.13 (upstream PR 3866) rewrote PGInterval.getValue() to skip zero-valued units, so
// relying on toString() silently changed every interval on the wire — "1 years 5 mons 5 days
// 0 hours 0 mins 0.0 secs" became "1 years 5 mons 5 days" on a driver bump alone. Render the six
// units here so the value does not depend on the driver version.
public class PgIntervalTypeConverter implements Converter<String> {
    @Override
    public String convert(Object value) {
        PGInterval interval = (PGInterval) value;
        DecimalFormat df = (DecimalFormat) NumberFormat.getInstance(Locale.US);
        df.applyPattern("0.0#####");
        return String.format(Locale.ROOT, "%d years %d mons %d days %d hours %d mins %s secs",
                interval.getYears(), interval.getMonths(), interval.getDays(),
                interval.getHours(), interval.getMinutes(), df.format(interval.getSeconds()));
    }
}
