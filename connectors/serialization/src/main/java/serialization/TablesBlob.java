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
}
