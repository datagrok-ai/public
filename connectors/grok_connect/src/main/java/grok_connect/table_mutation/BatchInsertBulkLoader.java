package grok_connect.table_mutation;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.List;

import grok_connect.providers.JdbcDataProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Default bulk loader for any prepared-statement JDBC provider: parses the CSV stream with
 * {@link CsvChunkParser}, binds each row via {@link MutationRunner#bindValue} using the declared dg
 * types, and {@code executeBatch}es every {@value #BATCH_SIZE} rows and at {@link #finish()}. Holds at
 * most one pending batch plus the parser's partial record, so heap stays flat regardless of payload size.
 */
public class BatchInsertBulkLoader implements BulkLoader {
    private static final Logger LOGGER = LoggerFactory.getLogger(BatchInsertBulkLoader.class);
    private static final int BATCH_SIZE = 10_000;

    private final JdbcDataProvider provider;
    private final InsertRows mutation;
    private final int columnCount;
    private final PreparedStatement statement;
    private final CsvChunkParser parser = new CsvChunkParser();
    private int pending;
    private int affected;

    public BatchInsertBulkLoader(JdbcDataProvider provider, Connection connection, InsertRows m) throws SQLException {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("Bulk insert requires a non-empty columns list");
        if (m.columnTypes == null || m.columnTypes.size() != m.columns.size())
            throw new MutationValidationException("Bulk insert requires columnTypes parallel to columns");
        this.provider = provider;
        this.mutation = m;
        this.columnCount = m.columns.size();
        this.statement = connection.prepareStatement(provider.insertSql(m));
    }

    @Override
    public void feed(byte[] csvChunk) throws SQLException {
        for (List<String> row : parser.feed(csvChunk))
            addRow(row);
    }

    @Override
    public MutationResult finish() throws SQLException {
        for (List<String> row : parser.finish())
            addRow(row);
        flush();
        closeQuietly();
        MutationResult result = new MutationResult();
        result.affectedRows = affected;
        return result;
    }

    @Override
    public void abort() {
        closeQuietly();
    }

    private void addRow(List<String> row) throws SQLException {
        if (row.size() != columnCount)
            throw new MutationValidationException("CSV row has " + row.size() + " field(s), expected " + columnCount);
        for (int c = 0; c < columnCount; c++)
            MutationRunner.bindValue(provider, statement, c + 1, row.get(c), mutation.columnTypes.get(c));
        statement.addBatch();
        if (++pending == BATCH_SIZE)
            flush();
    }

    private void flush() throws SQLException {
        if (pending == 0)
            return;
        affected += MutationRunner.sumBatchCounts(statement.executeBatch());
        pending = 0;
    }

    private void closeQuietly() {
        try {
            statement.close();
        } catch (SQLException e) {
            LOGGER.warn("Failed to close bulk insert statement", e);
        }
    }
}
