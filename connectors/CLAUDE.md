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
# Start REST server (port 1234, 4GB heap)
java -Xmx4g -classpath "grok_connect/target/grok_connect.jar:grok_connect/lib/*" grok_connect.GrokConnect

# Shell mode (CLI for testing queries)
java -classpath "grok_connect/target/grok_connect.jar:grok_connect/lib/*" grok_connect.GrokConnectShell
```

### Testing

```bash
# Run all tests
mvn test

# Run specific test class
mvn -Dtest=PostgresDataProviderTest test

# Run with TestContainers (requires Docker)
mvn verify
```

### Docker

```bash
# Build image
docker build -t grok_connect .

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
├── serialization/      # Binary data serialization module
│   ├── src/main/java/serialization/
│   │   ├── DataFrame.java               # Columnar data container
│   │   ├── Column.java                  # Typed column with data
│   │   ├── Types.java                   # Type constants
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

### REST API Endpoints

Defined in `GrokConnect.java` using Spark Java:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/query` | POST | Execute query, return binary DataFrame |
| `/query_socket` | WebSocket | Streaming query results |
| `/connectors` | GET | List available providers |
| `/schema` | POST | Get database schema |
| `/test_connection` | POST | Test connectivity |
| `/cancel` | POST | Cancel running query |

### Request/Response Flow

```
HTTP Request → Spark Route → Handler → Provider → JDBC → ResultSet
                                                            ↓
HTTP Response ← Binary Serialization ← DataFrame ← Column Managers
```

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

| Category | Technology | Version |
|----------|------------|---------|
| Language | Java | 8 |
| Language | Kotlin | 1.6.21 |
| Build | Maven | 3+ |
| REST | Spark Java | 2.9.4 |
| HTTP | Jetty | 9.4.x |
| JSON | Gson | 2.11.0 |
| Testing | JUnit 5 | 5.9.2 |
| Testing | TestContainers | 1.17.6 |
| Logging | SLF4J + Logback | 1.2.13 |
| AWS | AWS SDK | 2.20.52 |

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

2. **Add JDBC driver** to `lib/` directory

3. **Register provider** in `ProviderManager.java`:

```java
register(new MyDbDataProvider());
```

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

## Runtime Configuration

Default settings (can be overridden):

- **Port:** 1234
- **Heap:** 4GB (`-Xmx4g`)
- **JDBC Drivers:** `lib/` directory
- **Logging:** Logback (configurable via `logback.xml`)

## Notes

- Java 8 is required (some JDBC drivers don't support newer versions)
- JDBC drivers in `lib/` are not managed by Maven (pre-built)
- Kotlin is used only for SAP HANA provider and utilities
- TestContainers tests require Docker to be running
