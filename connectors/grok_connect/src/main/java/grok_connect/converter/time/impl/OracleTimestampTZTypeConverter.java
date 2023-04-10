package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import oracle.sql.TIMESTAMPTZ;
import oracle.sql.ZONEIDMAP;
import java.time.LocalDateTime;
import java.time.OffsetDateTime;
import java.time.ZoneId;
import java.time.ZoneOffset;
import java.util.Date;

public class OracleTimestampTZTypeConverter implements Converter<Date> {
    private static final byte REGION_ID_BIT = (byte) 0b1000_0000;

    @Override
    public Date convert(Object value) {
        OffsetDateTime offsetDateTime = timestamptzToOffsetDateTime((TIMESTAMPTZ) value);
        return Date.from(offsetDateTime.toInstant());
    }

    private OffsetDateTime timestamptzToOffsetDateTime(TIMESTAMPTZ dbData) {
        if (dbData == null) {
            return null;
        }
        byte[] bytes = dbData.toBytes();
        OffsetDateTime utc = extractUtc(bytes);
        if (isFixedOffset(bytes)) {
            ZoneOffset offset = extractOffset(bytes);
            return utc.withOffsetSameInstant(offset);
        }
        ZoneId zoneId = extractZoneId(bytes);
        return utc.atZoneSameInstant(zoneId).toOffsetDateTime();
    }

    private static OffsetDateTime extractUtc(byte[] bytes) {
        return OffsetDateTime.of(extractLocalDateTime(bytes), ZoneOffset.UTC);
    }

    private static boolean isFixedOffset(byte[] bytes) {
        return (bytes[11] & REGION_ID_BIT) == 0;
    }

    private static ZoneOffset extractOffset(byte[] bytes) {
        int hours = bytes[11] - 20;
        int minutes = bytes[12] - 60;
        if ((hours == 0) && (minutes == 0)) {
            return ZoneOffset.UTC;
        }
        return ZoneOffset.ofHoursMinutes(hours, minutes);
    }

    private static ZoneId extractZoneId(byte[] bytes) {
        // high order bits
        int regionCode = (bytes[11] & 0b1111111) << 6;
        // low order bits
        regionCode += (bytes[12] & 0b11111100) >> 2;
        String regionName = ZONEIDMAP.getRegion(regionCode);
        return ZoneId.of(regionName);
    }

    private static LocalDateTime extractLocalDateTime(byte[] bytes) {
        int year = ((Byte.toUnsignedInt(bytes[0]) - 100) * 100) + (Byte.toUnsignedInt(bytes[1]) - 100);
        int month = bytes[2];
        int dayOfMonth = bytes[3];
        int hour = bytes[4] - 1;
        int minute = bytes[5] - 1;
        int second = bytes[6] - 1;
        int nanoOfSecond = (Byte.toUnsignedInt(bytes[7]) << 24)
                | (Byte.toUnsignedInt(bytes[8]) << 16)
                | (Byte.toUnsignedInt(bytes[9]) << 8)
                | Byte.toUnsignedInt(bytes[10]);
        return LocalDateTime.of(year, month, dayOfMonth, hour, minute, second, nanoOfSecond);
    }
}
