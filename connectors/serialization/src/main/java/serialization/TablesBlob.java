package serialization;

import java.util.*;
import com.google.gson.*;
import com.google.gson.annotations.SerializedName;


/// Table binary file structure.
public class TablesBlob {
    @SerializedName("#type")
    public String type;

    /// Binary file version.
    public String version;

    /// Tables information.
    public TableInfo[] tables;

    /// Buffer.
    transient BufferAccessor buffer;

    /// Tables offsets in buffer.
    transient int[] tablesOffsets;

    /// Columns offsets in buffer, for each table.
    transient int[][] columnsOffsets;

    /// Total length of blob, in bytes.
    transient int length;

    /// No-arg constructor used by [fromByteArray].
    public TablesBlob() {
    }

    /// Class constructors.
    public TablesBlob(DataFrame[] tables) {
        type = "TablesBlob";
        version = BufferAccessor.VERSION;

        buffer = new BufferAccessor();
        this.tables = new TableInfo[tables.length];
        tablesOffsets = new int[tables.length];
        columnsOffsets = new int[tables.length][];

        for (int n = 0; n < tables.length; n++) {
            tablesOffsets[n] = buffer.bufPos;
            columnsOffsets[n] = buffer.writeDataFrame(tables[n]);
            this.tables[n] = new TableInfo(tables[n]);
        }

        length = buffer.bufPos;
    }

    /// Serializes the [TablesBlob] into array of bytes [Uint8List].
    public byte[] toByteArray() {
        GsonBuilder builder = new GsonBuilder();
        Gson gson = builder.create();

        int infoOffset = buffer.bufPos;

        // JSON with TablesBlob class
        buffer.writeString(gson.toJson(this));

        // Offsets (tables, columns)
        buffer.writeInt32List(tablesOffsets);
        for (int n = 0; n < columnsOffsets.length; n++)
            buffer.writeInt32List(columnsOffsets[n]);

        // Version string, 24 bytes
        String versionStr = (version.length() >= 14)
                ? version.substring(0, 14)
                : (version + (String.join("", Collections.nCopies(14 - version.length(), "_"))));
        buffer.writeString(versionStr);

        // Length of blob without tail, in bytes, 4 bytes
        buffer.writeInt32(buffer.bufPos - infoOffset);
        return buffer.toUint8List();
    }

    // Accepted version prefixes: 0.6.0 (current Dart) and 0.1.0 (Java writer).
    private static final String[] ACCEPTED_VERSIONS = {"0.6.0", "0.1.0"};

    /// Deserializes a [TablesBlob] from an array of bytes (tail-first).
    /// A truncated / malformed / out-of-range buffer fails with a structured
    /// RuntimeException rather than a raw ArrayIndexOutOfBoundsException - the
    /// WS mutation/query protocol sends one TablesBlob per binary frame, so a
    /// corrupt frame must surface cleanly.
    public static TablesBlob fromByteArray(byte[] bytes) {
        try {
            return fromByteArrayImpl(bytes);
        } catch (IndexOutOfBoundsException | NegativeArraySizeException e) {
            throw new RuntimeException("Malformed or truncated d42 buffer", e);
        }
    }

    private static TablesBlob fromByteArrayImpl(byte[] bytes) {
        TablesBlob blob = new TablesBlob();
        blob.buffer = new BufferAccessor(bytes);

        // Tail length lives in the last 4 bytes; it is the size of the metadata
        // region, so the JSON starts at length = total - (tail + 4).
        blob.buffer.bufPos = bytes.length - 4;
        int tail = blob.buffer.readInt32();
        blob.length = bytes.length - (tail + 4);
        blob.buffer.bufPos = blob.length;

        // JSON metadata with the TablesBlob class.
        TablesBlob parsed = new Gson().fromJson(blob.buffer.readString(), TablesBlob.class);
        blob.type = parsed.type;
        blob.version = parsed.version;
        blob.tables = parsed.tables;

        boolean accepted = false;
        if (blob.version != null)
            for (String v : ACCEPTED_VERSIONS)
                if (blob.version.startsWith(v)) {
                    accepted = true;
                    break;
                }
        if (!accepted)
            throw new RuntimeException("Unsupported d42 version " + blob.version
                    + "; this GrokConnect reads 0.6.0/0.1.0");

        // Offsets (tables, then columns per table).
        blob.tablesOffsets = blob.buffer.readInt32List();
        blob.columnsOffsets = new int[blob.tablesOffsets.length][];
        for (int n = 0; n < blob.tablesOffsets.length; n++)
            blob.columnsOffsets[n] = blob.buffer.readInt32List();

        return blob;
    }

    /// Returns the i-th [DataFrame] in the blob. Wraps out-of-range offsets /
    /// truncated column data into a clean RuntimeException (see [fromByteArray]).
    public DataFrame getTable(int idx) {
        try {
            buffer.bufPos = tablesOffsets[idx];
            return buffer.readDataFrame(columnsOffsets[idx]);
        } catch (IndexOutOfBoundsException | NegativeArraySizeException e) {
            throw new RuntimeException("Malformed or truncated d42 buffer", e);
        }
    }
}
