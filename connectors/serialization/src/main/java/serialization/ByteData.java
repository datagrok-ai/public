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
}
