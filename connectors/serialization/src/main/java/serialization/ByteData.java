package serialization;


// Class for convenient for casting content to different types.
public class ByteData
{
    byte[] buf;

    public ByteData(byte[] buf)
    {
        this.buf = buf;
    }

    public void setInt8(int idx, byte value)
    {
        buf[idx] = value;
    }

    public void setInt16(int idx, short value)
    {
        buf[idx] = (byte)value;
        buf[idx + 1] = (byte)(value >> 8);
    }

    public void setInt32(int idx, int value)
    {
        buf[idx] = (byte)value;
        buf[idx + 1] = (byte)(value >> 8);
        buf[idx + 2] = (byte)(value >> 16);
        buf[idx + 3] = (byte)(value >> 24);
    }

    public void setFloat32(int idx, float value)
    {
        int bits = Float.floatToIntBits(value);

        buf[idx] = (byte)bits;
        buf[idx + 1] = (byte)(bits >> 8);
        buf[idx + 2] = (byte)(bits >> 16);
        buf[idx + 3] = (byte)(bits >> 24);
    }

    public void setFloat64(int idx, double value)
    {
        long bits = Double.doubleToLongBits(value);

        buf[idx] = (byte)bits;
        buf[idx + 1] = (byte)(bits >> 8);
        buf[idx + 2] = (byte)(bits >> 16);
        buf[idx + 3] = (byte)(bits >> 24);
        buf[idx + 4] = (byte)(bits >> 32);
        buf[idx + 5] = (byte)(bits >> 40);
        buf[idx + 6] = (byte)(bits >> 48);
        buf[idx + 7] = (byte)(bits >> 56);
    }

    // Little-endian read accessors (mirror the set* methods above).

    public byte getInt8(int idx)
    {
        return buf[idx];
    }

    public int getUint8(int idx)
    {
        return buf[idx] & 0xFF;
    }

    public short getInt16(int idx)
    {
        return (short)((buf[idx] & 0xFF) | ((buf[idx + 1] & 0xFF) << 8));
    }

    public int getUint16(int idx)
    {
        return (buf[idx] & 0xFF) | ((buf[idx + 1] & 0xFF) << 8);
    }

    public int getInt32(int idx)
    {
        return (buf[idx] & 0xFF)
                | ((buf[idx + 1] & 0xFF) << 8)
                | ((buf[idx + 2] & 0xFF) << 16)
                | ((buf[idx + 3] & 0xFF) << 24);
    }

    public float getFloat32(int idx)
    {
        return Float.intBitsToFloat(getInt32(idx));
    }

    public double getFloat64(int idx)
    {
        long bits = 0;
        for (int k = 0; k < 8; k++)
            bits |= (long)(buf[idx + k] & 0xFF) << (8 * k);
        return Double.longBitsToDouble(bits);
    }

    // Reinterprets a raw byte array (little-endian) as a typed array.
    // Used by the archive-flag (zlib) decode branches.

    public static int[] toInt32List(byte[] b)
    {
        int[] r = new int[b.length / 4];
        for (int i = 0; i < r.length; i++)
            r[i] = (b[i * 4] & 0xFF)
                    | ((b[i * 4 + 1] & 0xFF) << 8)
                    | ((b[i * 4 + 2] & 0xFF) << 16)
                    | ((b[i * 4 + 3] & 0xFF) << 24);
        return r;
    }

    public static int[] toUint32List(byte[] b)
    {
        return toInt32List(b);
    }

    public static short[] toInt16List(byte[] b)
    {
        short[] r = new short[b.length / 2];
        for (int i = 0; i < r.length; i++)
            r[i] = (short)((b[i * 2] & 0xFF) | ((b[i * 2 + 1] & 0xFF) << 8));
        return r;
    }

    public static float[] toFloat32List(byte[] b)
    {
        int[] ints = toInt32List(b);
        float[] r = new float[ints.length];
        for (int i = 0; i < ints.length; i++)
            r[i] = Float.intBitsToFloat(ints[i]);
        return r;
    }

    public static double[] toFloat64List(byte[] b)
    {
        double[] r = new double[b.length / 8];
        for (int i = 0; i < r.length; i++) {
            long bits = 0;
            for (int k = 0; k < 8; k++)
                bits |= (long)(b[i * 8 + k] & 0xFF) << (8 * k);
            r[i] = Double.longBitsToDouble(bits);
        }
        return r;
    }
}
