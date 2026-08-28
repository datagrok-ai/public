package grok_connect.managers.string_column;

import grok_connect.managers.ColumnManager;
import grok_connect.resultset.ColumnMeta;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.postgresql.util.PGInterval;

// Guards the six-unit interval rendering against a driver bump: pgjdbc 42.7.13 dropped the
// zero-valued units from PGInterval.toString(), which changed every interval on the wire.
class PgIntervalTypeConverterTest {
    private static final ColumnMeta META = new ColumnMeta(java.sql.Types.OTHER, "interval", 0, 0, "interval", 1);

    private String convert(PGInterval interval) {
        ColumnManager<String> manager = new DefaultStringColumnManager();
        return manager.convert(interval, META);
    }

    @Test
    void keepsZeroValuedUnits() {
        Assertions.assertEquals("1 years 5 mons 5 days 0 hours 0 mins 0.0 secs",
                convert(new PGInterval(1, 5, 5, 0, 0, 0.0)));
    }

    @Test
    void rendersEveryUnit() {
        Assertions.assertEquals("2 years 3 mons 4 days 5 hours 6 mins 7.5 secs",
                convert(new PGInterval(2, 3, 4, 5, 6, 7.5)));
    }

    @Test
    void rendersZeroInterval() {
        Assertions.assertEquals("0 years 0 mons 0 days 0 hours 0 mins 0.0 secs",
                convert(new PGInterval(0, 0, 0, 0, 0, 0.0)));
    }

    @Test
    void rendersNegativeInterval() {
        Assertions.assertEquals("-1 years -2 mons -3 days -4 hours -5 mins -6.5 secs",
                convert(new PGInterval(-1, -2, -3, -4, -5, -6.5)));
    }
}
