package grok_connect.table_mutation;

import grok_connect.providers.JdbcDataProvider;
import grok_connect.providers.MsSqlDataProvider;
import grok_connect.providers.MySqlDataProvider;
import grok_connect.providers.OracleDataProvider;
import grok_connect.providers.PostgresDataProvider;
import grok_connect.providers.SnowflakeDataProvider;
import grok_connect.providers.ClickHouseProvider;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Upsert SQL-emission tests (connector-writes WO-4) — no Docker. Pins the exact dialect SQL per
 * provider, the chunk-size seam ({@link JdbcDataProvider#upsertBatchRows}), the matchKeys safety
 * rail, and the capability flags surfaced through descriptors.
 */
public class UpsertSqlTest {

    private UpsertRows sales() {
        UpsertRows m = new UpsertRows();
        m.tableName = "sales";
        m.columns = Arrays.asList("id", "region", "amount");
        m.columnTypes = Arrays.asList("int", "string", "double");
        m.matchKeys = Arrays.asList("id", "region");
        return m;
    }

    @DisplayName("Postgres: INSERT ... ON CONFLICT DO UPDATE, single row (addBatch)")
    @Test
    public void postgres() {
        UpsertRows m = sales();
        m.schema = "public";
        Assertions.assertEquals(
                "INSERT INTO \"public\".\"sales\" (\"id\", \"region\", \"amount\") VALUES (?, ?, ?) "
                        + "ON CONFLICT (\"id\", \"region\") DO UPDATE SET \"amount\" = EXCLUDED.\"amount\"",
                new PostgresDataProvider().upsertSql(m, 1));
    }

    @DisplayName("Postgres: all columns are keys -> DO NOTHING")
    @Test
    public void postgres_allKeys_doNothing() {
        UpsertRows m = sales();
        m.matchKeys = Arrays.asList("id", "region", "amount");
        Assertions.assertEquals(
                "INSERT INTO \"sales\" (\"id\", \"region\", \"amount\") VALUES (?, ?, ?) "
                        + "ON CONFLICT (\"id\", \"region\", \"amount\") DO NOTHING",
                new PostgresDataProvider().upsertSql(m, 1));
    }

    @DisplayName("MySQL: INSERT ... ON DUPLICATE KEY UPDATE, VALUES(col)")
    @Test
    public void mysql() {
        Assertions.assertEquals(
                "INSERT INTO `sales` (`id`, `region`, `amount`) VALUES (?, ?, ?) "
                        + "ON DUPLICATE KEY UPDATE `amount` = VALUES(`amount`)",
                new MySqlDataProvider().upsertSql(sales(), 1));
    }

    @DisplayName("MS SQL: MERGE over multi-row VALUES with a trailing semicolon")
    @Test
    public void mssql() {
        UpsertRows m = sales();
        m.schema = "dbo";
        Assertions.assertEquals(
                "MERGE INTO [dbo].[sales] AS t USING (VALUES (?, ?, ?), (?, ?, ?)) AS src ([id], [region], [amount]) "
                        + "ON (t.[id] = src.[id] AND t.[region] = src.[region]) "
                        + "WHEN MATCHED THEN UPDATE SET t.[amount] = src.[amount] "
                        + "WHEN NOT MATCHED THEN INSERT ([id], [region], [amount]) "
                        + "VALUES (src.[id], src.[region], src.[amount]);",
                new MsSqlDataProvider().upsertSql(m, 2));
    }

    @DisplayName("Oracle: MERGE ... USING (SELECT ? AS col FROM dual), one row")
    @Test
    public void oracle() {
        Assertions.assertEquals(
                "MERGE INTO \"sales\" t USING (SELECT ? AS \"id\", ? AS \"region\", ? AS \"amount\" FROM dual) src "
                        + "ON (t.\"id\" = src.\"id\" AND t.\"region\" = src.\"region\") "
                        + "WHEN MATCHED THEN UPDATE SET t.\"amount\" = src.\"amount\" "
                        + "WHEN NOT MATCHED THEN INSERT (\"id\", \"region\", \"amount\") "
                        + "VALUES (src.\"id\", src.\"region\", src.\"amount\")",
                new OracleDataProvider().upsertSql(sales(), 1));
    }

    @DisplayName("Snowflake: MERGE over multi-row VALUES (MS SQL shape, no semicolon)")
    @Test
    public void snowflake() {
        Assertions.assertEquals(
                "MERGE INTO \"sales\" AS t USING (VALUES (?, ?, ?), (?, ?, ?)) AS src (\"id\", \"region\", \"amount\") "
                        + "ON (t.\"id\" = src.\"id\" AND t.\"region\" = src.\"region\") "
                        + "WHEN MATCHED THEN UPDATE SET t.\"amount\" = src.\"amount\" "
                        + "WHEN NOT MATCHED THEN INSERT (\"id\", \"region\", \"amount\") "
                        + "VALUES (src.\"id\", src.\"region\", src.\"amount\")",
                new SnowflakeDataProvider().upsertSql(sales(), 2));
    }

    @DisplayName("Chunk seam: addBatch dialects emit 1 row/statement; MERGE dialects chunk under the param limit")
    @Test
    public void upsertBatchRows() {
        Assertions.assertEquals(1, new PostgresDataProvider().upsertBatchRows(3));
        Assertions.assertEquals(1, new MySqlDataProvider().upsertBatchRows(3));
        Assertions.assertEquals(1, new OracleDataProvider().upsertBatchRows(3));
        Assertions.assertEquals(500, new MsSqlDataProvider().upsertBatchRows(3)); // min(500, 2000/3)
        Assertions.assertEquals(4, new MsSqlDataProvider().upsertBatchRows(500)); // 2000/500
        Assertions.assertEquals(1, new MsSqlDataProvider().upsertBatchRows(5000)); // clamped to >= 1
        Assertions.assertEquals(500, new SnowflakeDataProvider().upsertBatchRows(3));
    }

    @DisplayName("Safety rail: empty/absent matchKeys and a key not in columns are refused")
    @Test
    public void matchKeys_validation() {
        UpsertRows noKeys = sales();
        noKeys.matchKeys = null;
        Assertions.assertThrows(MutationValidationException.class, () -> new PostgresDataProvider().upsertSql(noKeys, 1));
        UpsertRows strayKey = sales();
        strayKey.matchKeys = Arrays.asList("id", "missing");
        Assertions.assertThrows(MutationValidationException.class, () -> new MsSqlDataProvider().upsertSql(strayKey, 2));
        UpsertRows hostile = sales();
        hostile.matchKeys = Arrays.asList("id\") OR 1=1 --", "region");
        hostile.columns = Arrays.asList("id\") OR 1=1 --", "region", "amount");
        Assertions.assertThrows(MutationValidationException.class, () -> new PostgresDataProvider().upsertSql(hostile, 1));
    }

    @DisplayName("Capability: a provider with no upsert override throws (mapped to a capability 400)")
    @Test
    public void noOverride_capabilityError() {
        Assertions.assertThrows(UnsupportedOperationException.class,
                () -> new ClickHouseProvider().upsertSql(sales(), 1));
    }

    @DisplayName("Descriptor flags: supportsWrite defaulted centrally, denylist stays read-only, upsert per dialect")
    @Test
    public void descriptorFlags() {
        ProviderManager pm = new ProviderManager();
        // supportsWrite: auto-interpolation JDBC providers, denylist excluded
        Assertions.assertTrue(pm.getByName("Postgres").descriptor.supportsWrite);
        Assertions.assertTrue(pm.getByName("MySQL").descriptor.supportsWrite);
        Assertions.assertTrue(pm.getByName("MariaDB").descriptor.supportsWrite);
        Assertions.assertFalse(pm.getByName("Athena").descriptor.supportsWrite);
        Assertions.assertFalse(pm.getByName("BigQuery").descriptor.supportsWrite);
        Assertions.assertFalse(pm.getByName("Impala").descriptor.supportsWrite);
        Assertions.assertFalse(pm.getByName("PI").descriptor.supportsWrite);
        // supportsUpsert: set per dialect (MariaDB inherits from MySQL's shared descriptor)
        for (String type : new String[] {"Postgres", "MySQL", "MariaDB", "MS SQL", "Oracle", "Snowflake"})
            Assertions.assertTrue(pm.getByName(type).descriptor.supportsUpsert, type + " should support upsert");
        Assertions.assertTrue(pm.getByName("Postgres").descriptor.supportsGeneratedKeys);
        Assertions.assertFalse(pm.getByName("Oracle").descriptor.supportsGeneratedKeys); // flag not set for Oracle
    }
}
