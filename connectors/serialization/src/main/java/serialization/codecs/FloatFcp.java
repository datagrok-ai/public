package serialization.codecs;

import serialization.BufferAccessor;

// Decode-only port of ddt/lib/src/serialization/float_fcp.dart:5-323 (FCM/DFCM
// predictor + decode path only). Single-precision (float32) only.
//
// All predictor state is kept in Java `int` (32-bit), which truncates exactly
// like Dart's Int32List: the encode/decode only ever observe the low 32 bits of
// each value (both `_intToSingle` and the predictor tables reinterpret / store
// as int32), and Java `>>` is arithmetic like Dart `>>`, so the ports match the
// Dart arbitrary-precision-then-truncate arithmetic bit-for-bit.
public class FloatFcp {
    static final int logOfTableSize = 17;

    // Ports ddt float_fcp.dart:5-25.
    private static final class FcmPredictor {
        final int[] table;
        int hash = 0;
        final int logSize;

        FcmPredictor(int logSize) {
            this.logSize = logSize;
            table = new int[1 << logSize];
        }

        int predict() {
            return table[hash];
        }

        void update(int trueValue) {
            table[hash] = trueValue;
            hash = ((hash << 6) ^ (trueValue >> (32 - logSize))) & (table.length - 1);
        }
    }

    // Ports ddt float_fcp.dart:29-54.
    private static final class DfcmPredictor {
        final int[] table;
        int hash = 0;
        int lastValue = 0;
        final int logSize;

        DfcmPredictor(int logSize) {
            this.logSize = logSize;
            table = new int[1 << logSize];
        }

        int predict() {
            return table[hash] + lastValue;
        }

        void update(int trueValue) {
            table[hash] = trueValue - lastValue;
            hash = ((hash << 2) ^ ((trueValue - lastValue) >> (32 - logSize - 3))) & (table.length - 1);
            lastValue = trueValue;
        }
    }

    private FcmPredictor fcm;
    private DfcmPredictor dfcm;

    // Ports ddt float_fcp.dart:124-134.
    public float[] decode(BufferAccessor buf) {
        float[] values = new float[buf.readInt32()];
        fcm = new FcmPredictor(logOfTableSize);
        dfcm = new DfcmPredictor(logOfTableSize);
        for (int i = 0; i < values.length; i += 2)
            decodePair(buf, values, i);
        return values;
    }

    // Ports ddt float_fcp.dart:231-274.
    private void decodePair(BufferAccessor buf, float[] values, int i) {
        int prediction;
        int header = buf.readInt8() & 0xFF;

        // First value.
        if ((header & 0x80) != 0)
            prediction = dfcm.predict();
        else
            prediction = fcm.predict();

        int numZeroBytes = (header & 0x70) >> 4;
        int diff = readInt(buf, 4 - numZeroBytes);
        int actual = prediction ^ diff;

        fcm.update(actual);
        dfcm.update(actual);

        values[i] = Float.intBitsToFloat(actual);

        // Second value.
        if ((header & 0x08) != 0)
            prediction = dfcm.predict();
        else
            prediction = fcm.predict();

        numZeroBytes = header & 0x07;
        diff = readInt(buf, 4 - numZeroBytes);

        if (numZeroBytes == 3 && diff == 0)
            return;

        actual = prediction ^ diff;

        fcm.update(actual);
        dfcm.update(actual);

        values[i + 1] = Float.intBitsToFloat(actual);
    }

    // Ports ddt float_fcp.dart:289-322 (_readBytes + _toInt): little-endian bytes.
    private static int readInt(BufferAccessor buf, int length) {
        byte[] bytes = new byte[length];
        for (int n = 0; n < length; n++)
            bytes[n] = (byte) buf.readInt8();
        int result = 0;
        for (int i = length; i > 0; i--) {
            result = result << 8;
            result |= bytes[i - 1] & 0xff;
        }
        return result;
    }
}
