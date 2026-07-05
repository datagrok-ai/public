package serialization;

/// Column encoder additional archive types.
/// Deprecated on the write side (GROK-15721) - current writers always emit
/// ARCHIVE_TYPE_NONE. The reader still honors ARCHIVE_TYPE_ZLIB for old blobs.
public class ColumnEncoderArchiveType {
    public static final int ARCHIVE_TYPE_NONE = 0;
    public static final int ARCHIVE_TYPE_ZLIB = 1;
}
