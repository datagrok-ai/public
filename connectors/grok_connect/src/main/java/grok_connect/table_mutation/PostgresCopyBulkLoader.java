package grok_connect.table_mutation;

import java.sql.Connection;
import java.sql.SQLException;

import org.postgresql.PGConnection;
import org.postgresql.copy.CopyIn;
import org.postgresql.copy.CopyManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.DataFrame;

/**
 * Postgres bulk-insert fast path via {@code COPY <table> (<cols>) FROM STDIN (FORMAT csv)}: each decoded
 * d42 chunk is formatted into COPY-csv text ({@link CopyCsvFormatter}) and streamed straight into pgjdbc's
 * {@link CopyManager}. Insert-only — the provider selects it only for {@code mode == "insert"} &&
 * {@code allOrNothing} ({@link PostgresDataProvider}).
 */
public class PostgresCopyBulkLoader implements BulkLoader {
    private static final Logger LOGGER = LoggerFactory.getLogger(PostgresCopyBulkLoader.class);

    private final CopyIn copyIn;

    public PostgresCopyBulkLoader(Connection connection, String copySql) throws SQLException {
        CopyManager copyManager = connection.unwrap(PGConnection.class).getCopyAPI();
        this.copyIn = copyManager.copyIn(copySql);
    }

    @Override
    public void feed(DataFrame chunk) throws SQLException {
        byte[] csv = CopyCsvFormatter.format(chunk);
        copyIn.writeToCopy(csv, 0, csv.length);
    }

    @Override
    public MutationResult finish() throws SQLException {
        long rows = copyIn.endCopy();
        MutationResult result = new MutationResult();
        result.affectedRows = (int) rows;
        return result;
    }

    @Override
    public void abort() {
        if (!copyIn.isActive())
            return;
        try {
            copyIn.cancelCopy();
        } catch (SQLException e) {
            LOGGER.warn("Failed to cancel COPY operation", e);
        }
    }
}
