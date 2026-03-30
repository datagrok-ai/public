---
name: test-connectors
description: Run Java/Maven tests for the GrokConnect connectors project
argument-hint: "[all|<ClassName>|<ClassName>#<method>]"
---

# Test Connectors (GrokConnect)

Run Java/Maven tests for the GrokConnect connectors project.

## Usage

```
/test-connectors [target]
```

Where `target` is one of:
- `all` - Run all tests (default)
- `<ClassName>` - Run a specific test class (e.g., `PostgresDataProviderTest`, `SerializationTest`)
- `<ClassName>#<method>` - Run a specific test method

## Prerequisites

- **JDK 8** - Required by the project
- **Maven 3+** - Build tool
- **Docker** - Required for TestContainers-based integration tests (most provider tests)

## Instructions

When this skill is invoked, help the user run the appropriate tests for the connectors project.

### Project Location

```
connectors/
```

### Run All Tests

```bash
cd connectors && mvn test
```

### Run a Specific Test Class

```bash
cd connectors && mvn -Dtest=PostgresDataProviderTest test
```

### Run a Specific Test Method

```bash
cd connectors && mvn -Dtest=PostgresDataProviderTest#testMethodName test
```

### Build Without Tests (for quick compilation check)

```bash
cd connectors && mvn -Dmaven.test.skip=true package
```

### Run Only Serialization Tests (no Docker needed)

```bash
cd connectors && mvn -pl serialization test
```

### Run Only grok_connect Module Tests

```bash
cd connectors && mvn -pl grok_connect test
```

## Available Test Classes

### Provider Tests (require Docker for TestContainers)

Located in `grok_connect/src/test/java/grok_connect/providers/`:

| Test Class | Database |
|------------|----------|
| `PostgresDataProviderTest` | PostgreSQL |
| `MySqlDataProviderTest` | MySQL |
| `MariaDbDataProviderTest` | MariaDB |
| `MsSqlDataProviderTest` | MS SQL Server |
| `OracleDataProviderTest` | Oracle |
| `ClickHouseDataProviderTest` | ClickHouse |
| `MongoDbDataProviderTest` | MongoDB |
| `Neo4jDataProviderTest` | Neo4j |
| `CassandraDataProviderTest` | Cassandra |
| `Db2DataProviderTest` | DB2 |
| `VirtuosoDataProviderTest` | Virtuoso |
| `VerticaDataProviderTest` | Vertica |
| `Hive2DataProviderTest` | Hive2 |
| `ImpalaDataProviderTest` | Impala |
| `SnowflakeDataProviderTest` | Snowflake |
| `RedshiftDataProviderTest` | Redshift |
| `AthenaDataProviderTest` | Athena |
| `TeradataDataProviderTest` | Teradata |
| `PIDataProviderTest` | PI Data |

### Unit Tests (no Docker needed)

| Test Class | Location |
|------------|----------|
| `JdbcDataProviderTest` | `grok_connect/src/test/.../providers/` |
| `ComplexTypeConverterManagerTest` | `grok_connect/src/test/.../managers/complex_column/` |
| `TableQueryTest` | `grok_connect/src/test/.../table_query/` |
| `MsSqlTableQueryTest` | `grok_connect/src/test/.../table_query/` |
| `PostgresTableQueryTest` | `grok_connect/src/test/.../table_query/` |
| `SqlAnnotatorTest` | `grok_connect/src/test/.../utils/` |
| `SerializationTest` | `serialization/src/test/.../serialization/` |

## Behavior

1. Ask the user which tests to run if not specified
2. Check if Docker is running when integration tests are requested
3. Run the appropriate Maven test commands from `connectors/`
4. Report test results and any failures
