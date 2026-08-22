# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**GrokConnect** is a Java-based REST API server that bridges the Datagrok data analytics platform with 30+ databases and data sources. It handles database connections, schema introspection, query execution, and efficient binary data serialization.

**Location:** `public/connectors`

## Common Commands

### Building

```bash
# Build without tests (fastest)
mvn -Dmaven.test.skip=true package

# Build with tests
mvn package

# Clean build
mvn clean package

# Using build scripts
./grok_connect.sh        # Unix/Linux
grok_connect.cmd         # Windows
```

### Running

```bash
# Start REST server (port 1234, 4GB heap). Keep the jar FIRST on the classpath: with lib/* first,
# an SLF4J 1.x bundled in a driver jar shadows the shaded 2.x and GrokConnect.<clinit> crashes
# (grok_connect.cmd and the .sh shell-mode line still put lib/* first — watch out on a lib/ with
# extra driver jars). Run on JDK 8-17, not 18+ (see Notes).
java -Xmx4g -classpath "grok_connect/target/grok_connect.jar:grok_connect/lib/*" grok_connect.GrokConnect

# Shell mode (CLI for testing queries)
java -classpath "grok_connect/target/grok_connect.jar:grok_connect/lib/*" grok_connect.GrokConnectShell
```

### Testing

Use `/test-connectors` to run tests. Targeted suites that need no Docker:

- `mvn -pl serialization test` — d42 writer/reader, including the Java → Dart goldens
  (`D42JavaFixtureWriterTest`, compare mode: `src/test/resources/d42-java/*.d42` +
  `.expected.json` must match byte-for-byte; `-Dd42.regenerate=true` rewrites them after an
  intentional writer change — then run the Dart side,
  `cd core/shared/ddt && pub run test test/serialization/java_d42_fixture_test.dart`) and the
  Dart → Java goldens (`D42DartFixtureTest`).
- `mvn -pl grok_connect test -Dtest='SessionHandler*,GzipFrameTest,FastReadTest,QueryManagerTest'` —
  `SessionHandlerGzipChunkTest` drives the real `SessionHandler` state machine against a
  file-backed SQLite database (its driver ships with grok_connect; `SQLiteDataProvider` is put
  into `ProviderManager` and the static pool is set by reflection) with a mocked Jetty `Session`
  capturing announcements and frames. Reuse that pattern for any handler-level test that must run
  without a database container.

### Docker

```bash
# Build the main image (all providers on lean/patched drivers)
docker build -t grok_connect .

# Build the extended image (CVE-quarantined drivers: Neptune, Impala)
docker build --build-arg FLAVOR=extended --build-arg GROK_CONNECT_PROVIDERS="Neptune,Impala" -t grok_connect_extended .

# Run container
docker run -p 1234:1234 grok_connect
```

## Architecture Overview

### Core Pattern: Provider System

All database connectors extend `JdbcDataProvider` (or `DataProvider` for non-JDBC sources):

```
DataProvider (abstract)
├── descriptor          # Provider metadata (name, type, connection template)
├── getConnection()     # Establish database connection
├── execute()           # Run queries
├── getSchemas()        # Schema introspection
└── testConnection()    # Connectivity verification

JdbcDataProvider extends DataProvider
├── JDBC-specific implementations
├── Column type managers
└── Result set processing
```

### Two-Module Structure

```
connectors/
├── grok_connect/       # Main REST server application
│   ├── src/main/java/
│   │   ├── grok_connect/
│   │   │   ├── GrokConnect.java         # Main entry point, Spark routes
│   │   │   ├── GrokConnectShell.java    # CLI testing mode
│   │   │   ├── connectors_info/         # Data models (DataQuery, DataConnection, etc.)
│   │   │   ├── providers/               # 32 database provider implementations
│   │   │   ├── handlers/                # Request handlers (QueryHandler, etc.)
│   │   │   ├── managers/                # Column type managers
│   │   │   ├── resultset/               # ResultSet processing
│   │   │   ├── table_query/             # Structured query execution
│   │   │   └── utils/                   # Utilities
│   │   └── src/main/kotlin/             # Kotlin providers (SAP HANA)
│   ├── lib/                             # Pre-built JDBC drivers
│   └── pom.xml
│
├── serialization/      # Binary (d42) DataFrame serialization module: read + write
│   ├── src/main/java/serialization/
│   │   ├── DataFrame.java               # Columnar container; toByteArray() writes, fromByteArray() reads
│   │   ├── Column.java                  # Typed column: encode() writes, decode() reads
│   │   ├── BufferAccessor.java          # d42 primitives — write* and read* halves
│   │   ├── Types.java                   # Type constants
│   │   ├── codecs/                      # BitIntList/IntRle/IntSequencePattern read+write; FloatFcp/FloatRle/StringSquash decode-only
│   │   └── BigIntColumn.java, etc.      # Type-specific columns
│   └── pom.xml
│
└── pom.xml             # Parent POM
```

### Provider Implementation Pattern

Each database provider follows this structure:

```java
public class PostgresDataProvider extends JdbcDataProvider {
    // 1. Descriptor with connection parameters
    public DataSource descriptor = new DataSource(
        "Postgres",
        "jdbc:postgresql://{server}:{port}/{db}",
        PostgresDataProvider.class
    );

    // 2. Connection string building
    @Override
    public String getConnectionString(DataConnection conn) {
        return "jdbc:postgresql://" + conn.getServer() + ":" + conn.getPort() + "/" + conn.getDb();
    }

    // 3. Schema introspection (optional override)
    @Override
    public DataFrame getSchemas(DataConnection conn) { ... }

    // 4. Provider-specific SQL handling (optional)
    @Override
    protected String limitToSql(String query, int limit) { ... }
}
```

### Column Type Managers

Type-specific handling for JDBC result sets:

```
managers/
├── ColumnManager.java           # Base interface
├── bigint_column/               # BigInt handling
├── bool_column/                 # Boolean handling
├── datetime_column/             # DateTime handling
├── float_column/                # Float/Double handling
├── int_column/                  # Integer handling
├── string_column/               # String handling
└── complex_column/              # Complex types (JSON, arrays)
```

**Hierarchy:**

```
ColumnManager (interface)
├── DefaultIntColumnManager     — INTEGER, TINYINT, SMALLINT
├── DefaultBigIntColumnManager  — BIGINT, plus provider-specific variants
├── DefaultFloatColumnManager   — FLOAT, DOUBLE, DECIMAL, NUMERIC
├── DefaultDateTimeColumnManager — DATE, TIMESTAMP, TIME
├── DefaultBoolColumnManager    — BOOLEAN
├── DefaultStringColumnManager  — VARCHAR, CHAR, TEXT, UUID
└── DefaultComplexColumnManager — JSON, ARRAY, nested objects
```

**Fallback strategy:** If a type converter throws, the column manager silently downgrades
to `StringColumn`, preserving values at the cost of type information.

**Fast reads:** `ColumnManager.canReadFast(meta)` / `readFast(rs, i, column, meta)` let a manager
read a cell with the typed JDBC getter (`getInt`, `getDouble` + narrow / `getFloat`, `getBoolean`,
`getString` + `wasNull`) straight into the column's primitive `add` overload, skipping `getObject`
and the converter chain. `DefaultResultSetManager.init` resolves it once per column; the default
and the datetime / bigint / complex managers return `false`. The row loop in
`JdbcDataProvider.getResultSetSubDf` checks cancellation once per row. A converter failure turns
the column into a `StringColumn` and switches its fast path off.

**Provider-specific overrides:**
- `OracleBigIntColumnManager` — handles Oracle `NUMBER` types
- `SnowflakeBigIntColumnManager` — Snowflake-specific integer handling
- `MySqlMssqlBoolColumnManager` — MySQL/MSSQL `BIT` type (stored as 0/1)
- `SQLiteDateTimeColumnManager` — SQLite text-based dates
- Various DateTime type converters for vendor-specific temporal types
  (Oracle `TIMESTAMP WITH TIME ZONE`, Microsoft `DateTimeOffset`, etc.)

### SQL Type → Datagrok Type Mapping

Each provider declares a `typesMap` in its `DataSource` descriptor. Example (Postgres):

| SQL Type                                 | Datagrok Type       |
|------------------------------------------|---------------------|
| `int`, `int4`, `serial4`                 | `int`               |
| `int2`, `smallint`, `serial2`            | `int`               |
| `int8`, `bigint`, `bigserial`            | `bigint`            |
| `numeric`, `decimal`, `float4`, `float8` | `double`            |
| `varchar`, `text`, `char`, `uuid`        | `string`            |
| `boolean`                                | `bool`              |
| `date`, `timestamp`, `timestamptz`       | `datetime`          |
| `json`, `jsonb`, `ARRAY`                 | Complex → flattened |

### Query Parameter Handling

**Two interpolation modes:**

**A) Auto-Interpolation** (default, when `autoInterpolation()` returns `true`):

SQL uses `@paramName` placeholders, converted to `?` for `PreparedStatement`:

| Param Type     | JDBC Method       | Notes                                             |
|----------------|-------------------|---------------------------------------------------|
| `int`          | `setInt()`        |                                                   |
| `double`       | `setDouble()`     |                                                   |
| `bool`         | `setBoolean()`    |                                                   |
| `string`       | `setString()`     | UUID format → `setObject(UUID)` on some providers |
| `datetime`     | `setTimestamp()`  | With UTC Calendar                                 |
| `bigint`       | `setLong()`       |                                                   |
| `list(string)` | `createArrayOf()` | For `IN` clauses                                  |
| null           | `setNull()`       | With appropriate SQL type code                    |

**B) Manual Interpolation** (when `autoInterpolation()` returns `false`):
Values are escaped and interpolated directly into the SQL string. Used by some
NoSQL providers and specialized connectors.

**Pattern Matching Parameters:**

`@paramName(columnName)` syntax enables server-side filter construction:

| Pattern Type   | Converter                   | Generated SQL             |
|----------------|-----------------------------|---------------------------|
| Numeric range  | `numericPatternConverter`   | `col >= ? AND col <= ?`   |
| String match   | `stringPatternConverter`    | `col LIKE ?` or `col ~ ?` |
| DateTime range | `dateTimePatternConverter`  | `col >= ? AND col <= ?`   |
| Boolean        | `boolPatternConverter`      | `col = ?`                 |
| IN list        | (detected from list value)  | `col IN (?, ?, ?)`        |
| IS NULL        | (detected from null marker) | `col IS NULL`             |

### REST API Endpoints

Defined in `GrokConnect.java` using Spark Java:

| Endpoint           | Method    | Description                            |
|--------------------|-----------|----------------------------------------|
| `/query`           | POST      | Execute query, return binary DataFrame |
| `/query_socket`    | WebSocket | Streaming query results                |
| `/connectors`      | GET       | List available providers               |
| `/schema`          | POST      | Get database schema                    |
| `/test_connection` | POST      | Test connectivity                      |
| `/cancel`          | POST      | Cancel running query                   |

### Adding a New Auth Method to a Provider

1. **Define constant** in the provider (e.g., `private static final String MY_METHOD = "My Method"`)
2. **Add credential fields** to `descriptor.credentialsTemplate` with the method as category (4th arg)
3. **Handle in `getConnection()`** — branch on `#chosen-auth-method`
4. **Server side** (`credentials_service.dart`) — add branch in `readCredentials()` if server
   needs to transform credentials (e.g., mint tokens, exchange tokens)
5. **Server side** (`grok_server.dart`) — add method to decrypt condition in
   `getCredentialsForEntity()` if the flow needs credential resolution
6. **Connection pool** — bypass `ConnectionPool` for short-lived per-user tokens
   (use `DriverManager.getConnection()` directly, requires `Class.forName(driverClassName)`)

Credential field properties:
- `new Prop("password")` — masked input
- `new Prop("rsa")` with accept `.pem,.der` — file upload for keys
- `new Prop("textarea")` — multiline text
- Category parameter (4th arg) controls which auth method tab shows the field
- Same field name in different categories will collide — use unique names

### Request/Response Flow

```
HTTP Request → Spark Route → Handler → Provider → JDBC → ResultSet
                                                            ↓
HTTP Response ← Binary Serialization ← DataFrame ← Column Managers
```

### Streaming pipeline (`/query_socket`, `SessionHandler`)

`QUERY` → `QueryManager` → chunk 1 is read and serialized on the WebSocket thread and, with the
pipeline on (`PIPELINE`, see Runtime flags), `fetch(2)` = `getSubDF(2, bytesPerRow)` is submitted
to the pool right away (`GrokConnect.submitPoolTask` via `async`, MDC context propagated). Two
stages then run while chunk k is in flight: on `DATAFRAME PART SIZE RECEIVED` for chunk k the WS
thread chains `serialize(k+1)` onto the pending fetch (`fetchFuture.thenCompose`) and, unless that
fetch came back empty, `fetch(k+2)` after it, then sends chunk k's bytes; `serialize` =
`DataFrame.toBlob().toByteArray()` → `reportSerialized` → optional gzip → `Frame{bytes, gzipped,
rowCount}`. On `PART OK` the WS thread joins the serialize future and announces
`DATAFRAME PART SIZE: <n> [gzip=true]`; a non-`PART OK` reply re-announces the same frame once
(`QueryChunkNotSent` after that) and schedules nothing. `fetch(k+2)` is chained only by the next
SIZE RECEIVED, which bounds live data to two DataFrames and two Frames. With the pipeline off, one
pool task per chunk does fetch + serialize (the previous sequencing; the bytes are identical).

- Fetch sizing: the WS thread snapshots `QueryManager.getWireBytesPerRow()` (serialized bytes/row
  of the last chunk `serialize` reported) when it schedules a fetch and passes it as
  `getSubDF(n, wireBytesPerRow)`, which sizes the chunk after n from it (`memoryInBytes()/rowCount`
  only while it is 0) — so the size does not depend on which serialize happened to finish first.
- Separate data per stage: `getSubDF` for every chunk after the first calls
  `ResultSetManager.detach(nextColSize)` before reading — `DefaultResultSetManager` replaces each
  column with a fresh instance of the same class (`Column.getColumnForType`, class asserted;
  `BigIntColumn` is rebuilt as `bigint` and carries `setDowncastAllowed` + its sticky flag), so
  `serialize(k+1)` encodes columns the running `fetch(k+2)` never touches.
- Cancellation: `cancelTasks()` (from `onClose` and from the error cleanup) marks the query
  cancelled in `QueryMonitor` so the row loop throws `QueryCancelledByUser`, waits up to 5 s
  (`CANCEL_JOIN_TIMEOUT_MS`) for the running fetch, then `cancel(true)`s both futures; only after
  that does `QueryManager.close` close the result set and connection under them (`cancel(true)`
  alone never interrupts a running pool task).
- With `debugQuery`, `serialize` logs `COLUMN SIZES [...]` (`MISC` marker, per column
  `{name, type, enc, bytes, gz}`) right after the `DATAFRAME TO BYTEARRAY CONVERSION` END marker;
  the `columnTimings` option adds `ms` per column (re-encode of that column alone).

Datlas-side contract: `core/docs/DB_QUERY_FLOW.md` Steps 6-7. `SessionHandlerGzipChunkTest` (see
Testing) runs every chunk scenario with the pipeline off and on and asserts identical announce
tokens and bytes; `cancelledFetchPropagatesThroughThePipeline` and
`closeMidStreamCancelsBothStagesAndReleasesTheConnection` (a `GatedProvider` whose column reads
park on a latch) cover the cancel path: `onClose` waits for the parked row, the fetch ends in
`QueryCancelledByUser`, and the connection is closed.

Options Datlas puts into the `QUERY` JSON (`QueryManager`; all absent on old Datlas builds and
on the non-streaming `/query` path):

| Option                 | Effect                                                                                                                                                                                    |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `int8AsInt32`          | `BigIntColumn` reports and encodes as `int` while every value fits int32; sticky, lossless upgrade to `bigint` on the first overflow, plus the `ALLOW_COL_TYPE_CHANGE` tag on every chunk |
| `compressChunks`       | gzip each chunk inside the pool task and announce `gzip=true`; chunk 1 only when `initConnectFetchSize` is given and is not the number 100 (mirrors Datlas)                               |
| `gzipLevel`            | level for that gzip (default 1, clamped 1-9)                                                                                                                                              |
| `connectFetchSize`     | rows per chunk after the first, or a `"N MB"` byte target                                                                                                                                 |
| `initConnectFetchSize` | rows of chunk 1 (default 100)                                                                                                                                                             |
| `columnTimings`        | with `debugQuery`, adds a per-column `ms` (re-encode of that column alone) to the `COLUMN SIZES` debug log line                                                                           |

### Runtime flags (JVM system properties)

| Property                                | Default | Effect when `false`                                                                                                              |
|-----------------------------------------|---------|----------------------------------------------------------------------------------------------------------------------------------|
| `grok.connect.advancedEncoders`         | `true`  | `IntColumn` writes `int:raw` (id 1) only — also string category indices, the bigint int mode and, in effect, datetime components |
| `grok.connect.dateTimeComponentEncoder` | `true`  | `DateTimeColumn` writes `dateTime:microseconds` (id 3) only                                                                      |
| `grok.connect.pipelineFetch`            | `true`  | one pool task per chunk does fetch + serialize instead of the two-stage pipeline (same bytes on the wire)                        |

The encoder flags exist to roll a stand back to the pre-optimization wire bytes without a rebuild;
no reader needs them (every id was already decoded everywhere, `core/docs/D42_BINARY_FORMAT.md`).
`pipelineFetch` changes only task scheduling, never the bytes.

## Directory Structure

```
grok_connect/src/main/java/grok_connect/
├── GrokConnect.java              # Main class, REST routes
├── GrokConnectShell.java         # CLI mode
├── connectors_info/              # Data models
│   ├── DataConnection.java       # Connection parameters
│   ├── DataQuery.java            # Query with parameters
│   ├── DataSource.java           # Provider descriptor
│   ├── DbCredentials.java        # Authentication
│   └── FuncCall.java             # Function call wrapper
├── providers/                    # Database providers (32 implementations)
│   ├── JdbcDataProvider.java     # Base JDBC provider
│   ├── PostgresDataProvider.java
│   ├── MsSqlDataProvider.java
│   ├── OracleDataProvider.java
│   ├── MySqlDataProvider.java
│   ├── SnowflakeDataProvider.java
│   ├── BigQueryDataProvider.java
│   ├── DatabricksDataProvider.java
│   └── ...
├── handlers/                     # Request processing
│   ├── QueryHandler.java         # Query execution
│   └── SessionHandler.java       # Session management
├── managers/                     # Column type handlers
├── resultset/                    # ResultSet utilities
├── table_query/                  # Structured queries
├── utils/                        # Helpers
└── log/                          # Logging (QueryStreamAppender)
```

## Supported Databases

**Relational:** PostgreSQL, MySQL, MariaDB, Oracle, MS SQL Server, Teradata, Firebird, Vertica, ClickHouse, Denodo, SQLite, MS Access, HSQLDB

**Cloud Data Warehouses:** Snowflake, Google BigQuery, Amazon Redshift, Databricks, Amazon Athena

**NoSQL/Graph:** MongoDB, Cassandra, Neo4j, Amazon Neptune, Virtuoso, OrientDB, HBase, DynamoDB

**Big Data:** Hive, Hive2, Impala

**Enterprise:** SAP HANA (Kotlin), PI Data

## Key Technologies

| Category | Technology      | Version         |
|----------|-----------------|-----------------|
| Language | Java            | 8               |
| Language | Kotlin          | 1.6.21          |
| Build    | Maven           | 3+              |
| REST     | Spark Java      | 2.9.4           |
| HTTP     | Jetty           | 9.4.x           |
| JSON     | Gson            | 2.11.0          |
| Testing  | JUnit 5         | 5.9.2           |
| Testing  | TestContainers  | 1.17.6          |
| Logging  | SLF4J + Logback | 2.0.17 / 1.3.16 |
| AWS      | AWS SDK         | 2.52.0 (bom)    |

## Adding a New Database Provider

1. **Create provider class** in `providers/`:

    ```java
    public class MyDbDataProvider extends JdbcDataProvider {
        public DataSource descriptor = new DataSource(
            "MyDB",                          // Display name
            "jdbc:mydb://{server}:{port}/{db}", // Connection template
            MyDbDataProvider.class
        );

        public MyDbDataProvider() {
            descriptor.type = "MyDB";
            descriptor.defaultSchema = "public";

            // Define connection parameters
            descriptor.connectionProperties = Arrays.asList(
                new Property("server", Property.STRING_TYPE),
                new Property("port", Property.INT_TYPE, "3306"),
                new Property("db", Property.STRING_TYPE)
            );
        }

        @Override
        public String getConnectionString(DataConnection conn) {
            return "jdbc:mydb://" + conn.getServer() + ":" + conn.getPort() + "/" + conn.getDb();
        }

        // Override methods as needed for provider-specific behavior
    }
    ```

2. **Add JDBC driver** to `lib/` directory (or as a Maven dependency)

3. **Register provider** in `ProviderManager.java` — add its fully-qualified class name to the
   `PROVIDER_CLASSES` array. Registration is reflective: the provider is advertised on `/conn`
   only when its driver class is present and it passes the `providers.conf` allowlist baked
   into the image at build time (from the `GROK_CONNECT_PROVIDERS` build arg).

4. **Add tests** in `src/test/java/`:

    ```java
    public class MyDbDataProviderTest extends DataProviderTest {
        @Container
        public static GenericContainer<?> myDb = new GenericContainer<>("mydb:latest")
            .withExposedPorts(3306);

        // Test methods
    }
    ```

5. **Update CHANGELOG.md**

## Testing with TestContainers

Integration tests use TestContainers for real database instances:

```java
@Testcontainers
public class PostgresDataProviderTest {
    @Container
    public static PostgreSQLContainer<?> postgres = new PostgreSQLContainer<>("postgres:14");

    @Test
    void testQuery() {
        DataConnection conn = new DataConnection();
        conn.setServer(postgres.getHost());
        conn.setPort(postgres.getMappedPort(5432));
        // ...
    }
}
```

## Configuration Files

- `pom.xml` - Parent Maven configuration
- `grok_connect/pom.xml` - Main module dependencies
- `serialization/pom.xml` - Serialization module
- `Dockerfile` - Container build definition
- `grok_connect.sh` / `grok_connect.cmd` - Build/run scripts

## Key Design Principles

1. **Provider abstraction** - All databases implement common interface
2. **Efficient serialization** - Custom binary format for DataFrame transfer
3. **Type-specific handling** - Column managers for proper type conversion
4. **Connection pooling** - Reuse connections where appropriate
5. **Parameterized queries** - Prevent SQL injection
6. **Schema introspection** - Auto-discover tables, columns, types
7. **Streaming support** - WebSocket for large result sets

## Debugging Query Execution

### Debug Mode

Set `queryCall.options['debugFunc'] = true` to enable:
- GrokConnect runs 2 dry runs before the actual query (measures timing)
- Detailed timing logs for connection, fetch, serialization
- Log messages sent to Datlas via `LOG <json>` WebSocket messages

### Key Log Markers

| `EventType` marker                             | Log name                                    | What it measures                                                                                                           |
|------------------------------------------------|---------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|
| `CONNECTION_RECEIVE`                           | `CONNECTION RECEIVE`                        | Time to obtain JDBC connection                                                                                             |
| `RESULT_SET_PROCESSING_WITH_DATAFRAME_FILL`    | `RESULT SET PROCESSING WITH DATAFRAME FILL` | Row loop per chunk (driver fetch + column fill), with the df number                                                        |
| `DATAFRAME_TO_BYTEARRAY_CONVERSION`            | `DATAFRAME TO BYTEARRAY CONVERSION`         | Serialization time per chunk                                                                                               |
| `MISC` `COLUMN SIZES [...]`                    | `MISC`                                      | Per-column `{name, type, enc, bytes, gz}` of the serialized chunk (debug only; the `gz` re-compression is debug-only cost) |
| `COMPRESSION`                                  | `GZIP COMPRESSION`                          | Gzip per chunk when `compressChunks` is on; END message carries `<raw> -> <gz> bytes`                                      |
| `CHECKSUM_SEND`                                | `CHECKSUM SEND`                             | Time to send size announcement (`Data size:` is the gzipped length when gzipped)                                           |
| `SOCKET_BINARY_DATA_EXCHANGE`                  | `SOCKET BINARY DATA EXCHANGE`               | Time to send binary data                                                                                                   |
| `RESULT_SET_PROCESSING_WITHOUT_DATAFRAME_FILL` | `DRY RUN`                                   | Row loop of a dry run (`getObject` per cell, no fill)                                                                      |
| `DRY_RUN`                                      | `DRY RUN TOTAL DURATION`                    | Both dry runs (debug mode only)                                                                                            |

### Fetch Size Tuning

Override via query options:
- `connectFetchSize`: Override all chunks (rows or `"N MB"` for byte-based)
- `initConnectFetchSize`: Override first chunk only (rows)

Example: `queryCall.options['connectFetchSize'] = '5 MB'` targets 5 MB chunks.

Default adaptive behavior: initial 100 rows, then auto-calculated per chunk targeting
10 MB (`MAX_CHUNK_SIZE_BYTES`) from the serialized bytes/row of the previous chunk
(`lastWireBytesPerRow`; the in-memory estimate is used only before any chunk was serialized),
capped at 500,000 rows (`MAX_FETCH_SIZE`).

## Runtime Configuration

Default settings (can be overridden):

- **Port:** 1234
- **Heap:** 4GB (`-Xmx4g`)
- **JDBC Drivers:** `lib/` directory
- **Logging:** Logback (configurable via `logback.xml`)

## Notes

- Source/target level is Java 8 (`pom.xml`; some JDBC drivers don't support newer versions) and
  the Docker image runs on `datagrok/openjdk:8`. For local runs use JDK 8-17. **JDK 18+ breaks
  every hostname lookup**: the shade config excludes `META-INF/versions/**`
  (`grok_connect/pom.xml:164`) but keeps dnsjava's
  `META-INF/services/java.net.spi.InetAddressResolverProvider`, so the JVM tries to load a
  multi-release class that is no longer in the jar, `InetAddress` resolution throws and every
  JDBC connect times out (~30 s, pool `total=0`). Ticket candidate: exclude that service file (or
  keep `META-INF/versions`) in the shade config.
- Two image flavors from one codebase: `datagrok/grok_connect` (main) and `datagrok/grok_connect_extended`
  (opt-in; ships only the CVE-quarantined drivers — Amazon Neptune 3.0.3, Cloudera Impala). The `FLAVOR`
  build arg prunes `lib/`, and the `GROK_CONNECT_PROVIDERS` build arg (comma-separated `descriptor.type`
  values, empty = all) is baked into the image as `providers.conf` next to the jar — the allowlist is
  fixed at build time and cannot be changed with a runtime env var. ProviderManager also probes each
  provider's driver class and skips providers whose driver jar is absent, so `/conn` never advertises
  a provider that cannot connect.
- JDBC drivers in `lib/` are not managed by Maven (pre-built)
- Kotlin is used only for SAP HANA provider and utilities
- TestContainers tests require Docker to be running
