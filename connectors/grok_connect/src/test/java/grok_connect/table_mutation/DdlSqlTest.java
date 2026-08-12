package grok_connect.table_mutation;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataSource;
import grok_connect.providers.ClickHouseProvider;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.providers.MariaDbDataProvider;
import grok_connect.providers.MsSqlDataProvider;
import grok_connect.providers.MySqlDataProvider;
import grok_connect.providers.OracleDataProvider;
import grok_connect.providers.PostgresDataProvider;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * DDL SQL-emission tests (connector-writes WO-B6) — no Docker. Pins the exact per-dialect SQL for the
 * seven DDL ops, the dgToNativeType reverse maps and capability flags (ARCHITECTURE §3.3), the
 * DEFAULT-literal discipline, and the identifier/literal injection rails. Execution is WO-B7.
 */
public class DdlSqlTest {
    private final PostgresDataProvider postgres = new PostgresDataProvider();
    private final MySqlDataProvider mysql = new MySqlDataProvider();
    private final MariaDbDataProvider mariadb = new MariaDbDataProvider();
    private final MsSqlDataProvider mssql = new MsSqlDataProvider(); // default "[]" nameBrackets
    private final OracleDataProvider oracle = new OracleDataProvider();

    private ColumnSpec column(String name, String type, boolean nullable, String defaultValue) {
        ColumnSpec c = new ColumnSpec();
        c.name = name;
        c.type = type;
        c.nullable = nullable;
        c.defaultValue = defaultValue;
        return c;
    }

    private IndexSpec index(String name, boolean unique, String... columns) {
        IndexSpec i = new IndexSpec();
        i.name = name;
        i.unique = unique;
        i.columns = Arrays.asList(columns);
        return i;
    }

    private AlterTable alter(String schema, String action) {
        AlterTable m = new AlterTable();
        m.tableName = "orders";
        m.schema = schema;
        m.action = action;
        m.columnName = "qty";
        return m;
    }

    // ---- reverse type maps + capability flags ----

    private void assertNativeTypes(JdbcDataProvider p, String string, String intType, String bigint,
                                   String floatType, String bool, String datetime) {
        Assertions.assertEquals(string, p.nativeType("string"));
        Assertions.assertEquals(intType, p.nativeType("int"));
        Assertions.assertEquals(bigint, p.nativeType("bigint"));
        Assertions.assertEquals(floatType, p.nativeType("float"));
        Assertions.assertEquals(bool, p.nativeType("bool"));
        Assertions.assertEquals(datetime, p.nativeType("datetime"));
    }

    @DisplayName("dgToNativeType: the §3.3 table, pinned per dialect")
    @Test
    public void nativeTypes_pinnedPerDialect() {
        assertNativeTypes(postgres, "text", "int", "int8", "float8", "bool", "timestamp without time zone");
        assertNativeTypes(mysql, "text", "int", "bigint", "double", "boolean", "datetime(6)");
        assertNativeTypes(mssql, "nvarchar(max)", "int", "bigint", "float", "bit", "datetime2");
        assertNativeTypes(oracle, "varchar2(4000)", "number(10)", "number(19)", "binary_double", "number(1)", "timestamp");
    }

    @DisplayName("MariaDB inherits the MySQL map, flags, and MODIFY emission")
    @Test
    public void mariadb_inheritsMySql() {
        Assertions.assertEquals(mysql.descriptor.dgToNativeType, mariadb.descriptor.dgToNativeType);
        Assertions.assertTrue(mariadb.descriptor.supportsDdl);
        Assertions.assertFalse(mariadb.descriptor.supportsTransactionalDdl);
        AlterTable m = alter(null, "changeType");
        m.newType = "bigint";
        m.nullable = Boolean.TRUE;
        Assertions.assertEquals("ALTER TABLE `orders` MODIFY `qty` bigint NULL", mariadb.alterTableSql(m));
        Assertions.assertEquals("'\\\\'' inert'", mariadb.ddlLiteral("string", "\\' inert"));
    }

    @DisplayName("nativeType: missing map entry and DDL-less provider are validation errors")
    @Test
    public void nativeType_missing_refused() {
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.nativeType("qnum"));
        Assertions.assertThrows(MutationValidationException.class, () -> new ClickHouseProvider().nativeType("string"));
    }

    @DisplayName("ProviderManager: supportsDdl ⇒ dgToNativeType, v1 provider set + transactional flags")
    @Test
    public void providerManager_ddlInvariant() {
        ProviderManager manager = new ProviderManager();
        for (DataSource d : manager.getAllDescriptors())
            if (d.supportsDdl)
                Assertions.assertNotNull(d.dgToNativeType, "supportsDdl without dgToNativeType: " + d.type);
        Assertions.assertTrue(manager.getByName("Postgres").descriptor.supportsDdl);
        Assertions.assertTrue(manager.getByName("Postgres").descriptor.supportsTransactionalDdl);
        Assertions.assertTrue(manager.getByName("MySQL").descriptor.supportsDdl);
        Assertions.assertFalse(manager.getByName("MySQL").descriptor.supportsTransactionalDdl);
        Assertions.assertTrue(manager.getByName("MariaDB").descriptor.supportsDdl);
        Assertions.assertTrue(manager.getByName("MS SQL").descriptor.supportsDdl);
        Assertions.assertTrue(manager.getByName("MS SQL").descriptor.supportsTransactionalDdl);
        Assertions.assertTrue(manager.getByName("Oracle").descriptor.supportsDdl);
        Assertions.assertFalse(manager.getByName("Oracle").descriptor.supportsTransactionalDdl);
        Assertions.assertFalse(manager.getByName("MongoDB").descriptor.supportsDdl);
        Assertions.assertNull(manager.getByName("MongoDB").descriptor.dgToNativeType);
    }

    @DisplayName("/conn wire shape: PG carries dgToNativeType/supportsDdl/supportsTransactionalDdl; a DDL-less provider has no map")
    @Test
    public void descriptorWire_ddlFields() {
        JsonObject pg = JsonParser.parseString(GrokConnect.gson.toJson(postgres.descriptor)).getAsJsonObject();
        Assertions.assertTrue(pg.get("supportsDdl").getAsBoolean());
        Assertions.assertTrue(pg.get("supportsTransactionalDdl").getAsBoolean());
        Assertions.assertEquals(6, pg.getAsJsonObject("dgToNativeType").size());
        Assertions.assertEquals("int8", pg.getAsJsonObject("dgToNativeType").get("bigint").getAsString());
        JsonObject clickhouse = JsonParser.parseString(
                GrokConnect.gson.toJson(new ClickHouseProvider().descriptor)).getAsJsonObject();
        Assertions.assertFalse(clickhouse.get("supportsDdl").getAsBoolean());
        Assertions.assertFalse(clickhouse.has("dgToNativeType"));
    }

    // ---- ddlLiteral ----

    @DisplayName("ddlLiteral: parse-validated numerics, bools, quoted strings/datetimes, NULL")
    @Test
    public void ddlLiteral_typedValues() {
        Assertions.assertEquals("NULL", postgres.ddlLiteral("string", null));
        Assertions.assertEquals("0", postgres.ddlLiteral("int", "0"));
        Assertions.assertEquals("9007199254740995", postgres.ddlLiteral("bigint", "9007199254740995"));
        Assertions.assertEquals("0.5", postgres.ddlLiteral("float", "0.5"));
        Assertions.assertEquals("true", postgres.ddlLiteral("bool", "TRUE"));
        Assertions.assertEquals("false", postgres.ddlLiteral("bool", "yes"));
        Assertions.assertEquals("'it''s ok'", postgres.ddlLiteral("string", "it's ok"));
        Assertions.assertEquals("'2026-01-01T00:00:00.000Z'", postgres.ddlLiteral("datetime", "2026-01-01T00:00:00.000Z"));
    }

    @DisplayName("ddlLiteral: MS SQL / Oracle emit 1/0 for bool")
    @Test
    public void ddlLiteral_boolDialects() {
        Assertions.assertEquals("1", mssql.ddlLiteral("bool", "true"));
        Assertions.assertEquals("0", mssql.ddlLiteral("bool", "false"));
        Assertions.assertEquals("1", oracle.ddlLiteral("bool", "true"));
        Assertions.assertEquals("0", oracle.ddlLiteral("bool", "false"));
    }

    @DisplayName("ddlLiteral: hostile numeric refused; hostile string stays an inert quoted literal")
    @Test
    public void ddlLiteral_injection() {
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.ddlLiteral("int", "0; DROP TABLE orders"));
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.ddlLiteral("float", "Infinity"));
        Assertions.assertEquals("'''); DROP TABLE orders; --'",
                postgres.ddlLiteral("string", "'); DROP TABLE orders; --"));
    }

    @DisplayName("ddlLiteral: MySQL doubles backslashes too — \\' cannot escape out of the literal; PG stays quote-doubling only")
    @Test
    public void ddlLiteral_mysqlBackslashEscaping() {
        String hostile = "\\' ; DROP TABLE x -- ";
        // backslash doubled AND quote doubled: MySQL reads \\ as a literal backslash, '' as a literal quote
        Assertions.assertEquals("'\\\\'' ; DROP TABLE x -- '", mysql.ddlLiteral("string", hostile));
        // Postgres (standard_conforming_strings): backslash is not an escape — quote doubling alone suffices
        Assertions.assertEquals("'\\'' ; DROP TABLE x -- '", postgres.ddlLiteral("string", hostile));
    }

    // ---- CreateTable ----

    @DisplayName("createTableSql: Postgres — columns with defaults, PRIMARY KEY, per-IndexSpec statements, auto-named index")
    @Test
    public void createTable_postgres() {
        CreateTable m = new CreateTable();
        m.tableName = "orders";
        m.schema = "public";
        m.columns = new ArrayList<>();
        m.columns.add(column("id", "int", false, null));
        m.columns.add(column("big_id", "bigint", true, null));
        m.columns.add(column("amount", "float", true, "0.5"));
        m.columns.add(column("active", "bool", false, "true"));
        m.columns.add(column("note", "string", true, "it's ok"));
        m.columns.add(column("created", "datetime", true, "2026-01-01T00:00:00.000Z"));
        m.primaryKey = Collections.singletonList("id");
        m.indexes = Arrays.asList(index(null, false, "created"), index("ux_orders_note_created", true, "note", "created"));
        Assertions.assertEquals(Arrays.asList(
                "CREATE TABLE \"public\".\"orders\" (\"id\" int NOT NULL, \"big_id\" int8, "
                        + "\"amount\" float8 DEFAULT 0.5, \"active\" bool DEFAULT true NOT NULL, "
                        + "\"note\" text DEFAULT 'it''s ok', "
                        + "\"created\" timestamp without time zone DEFAULT '2026-01-01T00:00:00.000Z', "
                        + "PRIMARY KEY (\"id\"))",
                "CREATE INDEX \"ix_orders_created\" ON \"public\".\"orders\" (\"created\")",
                "CREATE UNIQUE INDEX \"ux_orders_note_created\" ON \"public\".\"orders\" (\"note\", \"created\")"),
                postgres.createTableSql(m));
    }

    @DisplayName("createTableSql: ifNotExists on Postgres/MySQL; MS SQL and Oracle emit a plain CREATE")
    @Test
    public void createTable_ifNotExists() {
        CreateTable m = new CreateTable();
        m.tableName = "orders";
        m.ifNotExists = true;
        m.columns = Collections.singletonList(column("id", "int", true, null));
        Assertions.assertEquals(Collections.singletonList("CREATE TABLE IF NOT EXISTS \"orders\" (\"id\" int)"),
                postgres.createTableSql(m));
        Assertions.assertEquals(Collections.singletonList("CREATE TABLE IF NOT EXISTS `orders` (`id` int)"),
                mysql.createTableSql(m));
        Assertions.assertEquals(Collections.singletonList("CREATE TABLE [orders] ([id] int)"),
                mssql.createTableSql(m));
        Assertions.assertEquals(Collections.singletonList("CREATE TABLE \"orders\" (\"id\" number(10))"),
                oracle.createTableSql(m));
    }

    @DisplayName("createTableSql: empty columns refused; hostile column name refused before SQL is built")
    @Test
    public void createTable_validation() {
        CreateTable m = new CreateTable();
        m.tableName = "orders";
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.createTableSql(m));
        m.columns = Collections.singletonList(column("x\"; drop table t; --", "int", true, null));
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.createTableSql(m));
    }

    // ---- AlterTable ----

    @DisplayName("alterTableSql addColumn: PG ADD COLUMN with DEFAULT + NOT NULL; MS SQL / Oracle ADD without COLUMN")
    @Test
    public void alterTable_addColumn() {
        AlterTable m = alter("public", "addColumn");
        m.column = column("qty", "int", false, "0");
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" ADD COLUMN \"qty\" int DEFAULT 0 NOT NULL",
                postgres.alterTableSql(m));
        AlterTable ms = alter("dbo", "addColumn");
        ms.column = column("note", "string", true, null);
        Assertions.assertEquals("ALTER TABLE [dbo].[orders] ADD [note] nvarchar(max)", mssql.alterTableSql(ms));
        AlterTable ora = alter(null, "addColumn");
        ora.column = column("note", "string", true, null);
        Assertions.assertEquals("ALTER TABLE \"orders\" ADD \"note\" varchar2(4000)", oracle.alterTableSql(ora));
    }

    @DisplayName("alterTableSql dropColumn: Postgres")
    @Test
    public void alterTable_dropColumn() {
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" DROP COLUMN \"qty\"",
                postgres.alterTableSql(alter("public", "dropColumn")));
    }

    @DisplayName("alterTableSql renameColumn: PG/MySQL RENAME COLUMN (MySQL 8+/MariaDB 10.5+); MS SQL sp_rename")
    @Test
    public void alterTable_renameColumn() {
        AlterTable m = alter("public", "renameColumn");
        m.newName = "quantity";
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" RENAME COLUMN \"qty\" TO \"quantity\"",
                postgres.alterTableSql(m));
        AlterTable my = alter(null, "renameColumn");
        my.newName = "quantity";
        Assertions.assertEquals("ALTER TABLE `orders` RENAME COLUMN `qty` TO `quantity`", mysql.alterTableSql(my));
        AlterTable ms = alter("dbo", "renameColumn");
        ms.newName = "quantity";
        Assertions.assertEquals("EXEC sp_rename 'dbo.orders.qty', 'quantity', 'COLUMN'", mssql.alterTableSql(ms));
    }

    @DisplayName("alterTableSql changeType: ALTER..TYPE (PG), MODIFY (MySQL/Oracle), ALTER COLUMN (MS SQL); PG/Oracle preserve nullability natively — no nullable needed")
    @Test
    public void alterTable_changeType() {
        AlterTable m = alter("public", "changeType");
        m.newType = "bigint";
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" ALTER COLUMN \"qty\" TYPE int8",
                postgres.alterTableSql(m));
        AlterTable ora = alter(null, "changeType");
        ora.newType = "bigint";
        Assertions.assertEquals("ALTER TABLE \"orders\" MODIFY \"qty\" number(19)", oracle.alterTableSql(ora));
    }

    @DisplayName("alterTableSql changeType on MySQL / MS SQL restates the column definition: requires explicit nullable, refused without it")
    @Test
    public void alterTable_changeType_nullabilityRestatingDialects() {
        AlterTable my = alter(null, "changeType");
        my.newType = "bigint";
        Assertions.assertThrows(MutationValidationException.class, () -> mysql.alterTableSql(my));
        my.nullable = Boolean.FALSE;
        Assertions.assertEquals("ALTER TABLE `orders` MODIFY `qty` bigint NOT NULL", mysql.alterTableSql(my));
        my.nullable = Boolean.TRUE;
        Assertions.assertEquals("ALTER TABLE `orders` MODIFY `qty` bigint NULL", mysql.alterTableSql(my));
        AlterTable ms = alter("dbo", "changeType");
        ms.newType = "bigint";
        Assertions.assertThrows(MutationValidationException.class, () -> mssql.alterTableSql(ms));
        ms.nullable = Boolean.FALSE;
        Assertions.assertEquals("ALTER TABLE [dbo].[orders] ALTER COLUMN [qty] bigint NOT NULL", mssql.alterTableSql(ms));
        ms.nullable = Boolean.TRUE;
        Assertions.assertEquals("ALTER TABLE [dbo].[orders] ALTER COLUMN [qty] bigint NULL", mssql.alterTableSql(ms));
    }

    @DisplayName("alterTableSql setNullable: PG SET/DROP NOT NULL; Oracle MODIFY without type")
    @Test
    public void alterTable_setNullable() {
        AlterTable m = alter("public", "setNullable");
        m.nullable = Boolean.FALSE;
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" ALTER COLUMN \"qty\" SET NOT NULL",
                postgres.alterTableSql(m));
        m.nullable = Boolean.TRUE;
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" ALTER COLUMN \"qty\" DROP NOT NULL",
                postgres.alterTableSql(m));
        AlterTable ora = alter(null, "setNullable");
        ora.nullable = Boolean.FALSE;
        Assertions.assertEquals("ALTER TABLE \"orders\" MODIFY \"qty\" NOT NULL", oracle.alterTableSql(ora));
    }

    @DisplayName("alterTableSql setNullable on MySQL / MS SQL restates the type: requires newType, refused without it")
    @Test
    public void alterTable_setNullable_typeRestatingDialects() {
        AlterTable my = alter(null, "setNullable");
        my.nullable = Boolean.FALSE;
        Assertions.assertThrows(MutationValidationException.class, () -> mysql.alterTableSql(my));
        my.newType = "int";
        Assertions.assertEquals("ALTER TABLE `orders` MODIFY `qty` int NOT NULL", mysql.alterTableSql(my));
        AlterTable ms = alter("dbo", "setNullable");
        ms.nullable = Boolean.TRUE;
        Assertions.assertThrows(MutationValidationException.class, () -> mssql.alterTableSql(ms));
        ms.newType = "int";
        Assertions.assertEquals("ALTER TABLE [dbo].[orders] ALTER COLUMN [qty] int NULL", mssql.alterTableSql(ms));
    }

    @DisplayName("alterTableSql: unknown action and missing per-action fields are validation errors")
    @Test
    public void alterTable_validation() {
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.alterTableSql(alter("public", "explode")));
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.alterTableSql(alter("public", null)));
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.alterTableSql(alter("public", "addColumn"))); // column == null
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.alterTableSql(alter("public", "renameColumn"))); // newName == null
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.alterTableSql(alter("public", "changeType"))); // newType == null
        Assertions.assertThrows(MutationValidationException.class,
                () -> postgres.alterTableSql(alter("public", "setNullable"))); // nullable == null — absent != false
    }

    // ---- CreateIndex / DropIndex ----

    @DisplayName("createIndexSql: auto-named and unique named indexes")
    @Test
    public void createIndex() {
        CreateIndex m = new CreateIndex();
        m.tableName = "orders";
        m.schema = "public";
        m.columns = Arrays.asList("status", "created");
        Assertions.assertEquals("CREATE INDEX \"ix_orders_status_created\" ON \"public\".\"orders\" (\"status\", \"created\")",
                postgres.createIndexSql(m));
        m.indexName = "ux_orders_status";
        m.unique = true;
        Assertions.assertEquals("CREATE UNIQUE INDEX \"ux_orders_status\" ON \"public\".\"orders\" (\"status\", \"created\")",
                postgres.createIndexSql(m));
    }

    @DisplayName("createIndexSql: hostile index name and empty columns refused")
    @Test
    public void createIndex_validation() {
        CreateIndex m = new CreateIndex();
        m.tableName = "orders";
        m.columns = Collections.singletonList("status");
        m.indexName = "ix\"; drop table t; --";
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.createIndexSql(m));
        m.indexName = null;
        m.columns = new ArrayList<>();
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.createIndexSql(m));
    }

    @DisplayName("dropIndexSql: schema-qualified (PG), DROP INDEX .. ON table (MySQL / MS SQL)")
    @Test
    public void dropIndex() {
        DropIndex m = new DropIndex();
        m.tableName = "orders";
        m.schema = "public";
        m.indexName = "ix_orders_status";
        Assertions.assertEquals("DROP INDEX \"public\".\"ix_orders_status\"", postgres.dropIndexSql(m));
        DropIndex my = new DropIndex();
        my.tableName = "orders";
        my.indexName = "ix_orders_status";
        Assertions.assertEquals("DROP INDEX `ix_orders_status` ON `orders`", mysql.dropIndexSql(my));
        DropIndex ms = new DropIndex();
        ms.tableName = "orders";
        ms.schema = "dbo";
        ms.indexName = "ix_orders_status";
        Assertions.assertEquals("DROP INDEX [ix_orders_status] ON [dbo].[orders]", mssql.dropIndexSql(ms));
    }

    // ---- AddKey ----

    @DisplayName("addKeySql: primary auto-named pk_<table>; foreign auto-named fk_<table>_<col> with REFERENCES")
    @Test
    public void addKey() {
        AddKey pk = new AddKey();
        pk.tableName = "orders";
        pk.schema = "public";
        pk.keyType = "primary";
        pk.columns = Collections.singletonList("id");
        Assertions.assertEquals("ALTER TABLE \"public\".\"orders\" ADD CONSTRAINT \"pk_orders\" PRIMARY KEY (\"id\")",
                postgres.addKeySql(pk));
        AddKey fk = new AddKey();
        fk.tableName = "order_items";
        fk.schema = "public";
        fk.keyType = "foreign";
        fk.columns = Collections.singletonList("order_id");
        fk.refTable = "public.orders";
        fk.refColumns = Collections.singletonList("id");
        Assertions.assertEquals("ALTER TABLE \"public\".\"order_items\" ADD CONSTRAINT \"fk_order_items_order_id\" "
                        + "FOREIGN KEY (\"order_id\") REFERENCES \"public\".\"orders\" (\"id\")",
                postgres.addKeySql(fk));
    }

    @DisplayName("addKeySql: bad keyType, missing columns, foreign without refTable/refColumns, hostile refTable refused")
    @Test
    public void addKey_validation() {
        AddKey m = new AddKey();
        m.tableName = "orders";
        m.keyType = "unique";
        m.columns = Collections.singletonList("id");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.addKeySql(m));
        m.keyType = "primary";
        m.columns = new ArrayList<>();
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.addKeySql(m));
        m.keyType = "foreign";
        m.columns = Collections.singletonList("order_id");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.addKeySql(m)); // no refTable
        m.refTable = "orders\"; drop table t; --";
        m.refColumns = Collections.singletonList("id");
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.addKeySql(m));
    }

    // ---- DropTable / TruncateTable ----

    @DisplayName("dropTableSql / truncateTableSql: fully-qualified addressing, hostile table refused")
    @Test
    public void dropAndTruncate() {
        DropTable drop = new DropTable();
        drop.tableName = "orders";
        drop.schema = "public";
        Assertions.assertEquals("DROP TABLE \"public\".\"orders\"", postgres.dropTableSql(drop));
        TruncateTable truncate = new TruncateTable();
        truncate.tableName = "orders";
        truncate.schema = "public";
        Assertions.assertEquals("TRUNCATE TABLE \"public\".\"orders\"", postgres.truncateTableSql(truncate));
        drop.tableName = "orders\"; DROP TABLE orders; --";
        Assertions.assertThrows(MutationValidationException.class, () -> postgres.dropTableSql(drop));
    }
}
