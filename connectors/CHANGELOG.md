# Grok Connect changelog

# 2.8.4

* GROK-20712: Fixed the main image registering zero providers — the newline `printf` writes into `providers.conf` for the default empty allowlist parsed as an empty set instead of "all providers", so every DB test failed with "Provider Postgres not found" since 2.8.3

# 2.8.3

* GROK-20712: Provider allowlist is now baked into the image at build time (`providers.conf` written by the Dockerfile `FLAVOR` layer from the `GROK_CONNECT_PROVIDERS` build arg) instead of being read from a runtime env variable, so the main/extended provider partition cannot be overridden at deploy time

# 2.8.2

* GrokConnect: ClickHouse: Fix "Execution failed" against ClickHouse >= 26

# 2.8.1

* GROK-18695: jackson 2.18.8 → 2.18.9 (pom bom + BigQuery companion lib jars; CVE-2026-59889, GHSA-mhm7-754m-9p8w)
* GROK-18695: Removed the dead Hive1 (HiveServer1) dependency tree (`hive-jdbc:0.7.1-cdh3u6` pom-type + `cassandra-thrift`) — the provider has shipped without its driver class and is probe-de-advertised; drops commons-lang 2.5 and the 2011-era transitives from the shaded jar

# 2.8.0

* GROK-18695: Split into two images built from one codebase: `datagrok/grok_connect` (main — all providers on lean or patched drivers) and `datagrok/grok_connect_extended` (opt-in; quarantines the drivers with unfixable CVE surface: Amazon Neptune 3.0.3 and Cloudera Impala). Docker `FLAVOR` build arg prunes `lib/`; runtime provider set = `GROK_CONNECT_PROVIDERS` env allowlist (empty = all) intersected with a driver-presence probe, so a provider whose driver jar is absent is never advertised on `/conn`
* ProviderManager: reflective, fault-tolerant provider registration (a broken/missing provider logs a warning instead of failing startup). HBase and Hive (HiveServer1) are no longer advertised — both have shipped without their driver classes (Phoenix query server / `org.apache.hadoop.hive.jdbc.HiveDriver`) and could never connect; Hive2 covers all modern Hive servers
* GROK-18695: Driver refresh sweep (keep providers in main image on patched drivers): hive-jdbc + hive-standalone-metastore-common 4.0.0-alpha-2 → 4.0.1 (last Java-8 build; clears hive-service CVE-2024-23945) + explicit hive-service/hive-common/hive-shims-common/hive-service-rpc/libthrift runtime set (the 4.0.1 GA jar no longer declares its runtime deps), hadoop-common 3.4.3, zookeeper pin 3.9.5, postgresql 42.7.13, mysql-connector-j 9.7.0, mariadb 3.5.10, cassandra java-driver 4.17.0, mongodb-driver 3.12.14, sqlite-jdbc 3.53.2.1, Oracle ojdbc8/xdb/xmlparserv2 23.26.3.0.0, SAP ngdbc 2.29.7, terajdbc 20.00.00.58, jackcess 4.0.11, org.json 20260719, dnsjava 3.6.5, nimbus-jose-jwt 9.48, gson 2.14.0, BouncyCastle jdk15on 1.70 → jdk18on 1.85.x
* GROK-18695: jetty 9.4.57 → 9.4.58.v20250814 (newest OSS 9.4; the remaining CVE-2026-2332 has no OSS fix on the Java-8 line and is covered by a reviewed VEX not_affected entry — grok_connect is cluster-internal with datlas as its only client, no HTTP intermediary to desynchronize)
* GROK-18695: Widened the hive-jdbc jline exclusion to `org.jline:*` (drops jline-remote-telnet 3.22.0, GHSA HIGHs)
* GROK-18695: Athena: migrated from the discontinued Simba 2.x driver (`AthenaJDBC42-2.0.35.1000.jar`, embedded jackson 2.14.0 / log4j 2.17.1) to the first-party AWS Athena JDBC v3 driver 3.8.0 (lean jar + AWS SDK v2 from Maven, so netty stays on the patched bom pin). Driver class `com.amazon.athena.jdbc.AthenaDriver`, URL `jdbc:athena://`; credentials now passed as `User`/`Password`/`SessionToken`; `S3OutputLocation`/`Schema`/`S3OutputEncOption` mapped to v3 `OutputLocation`/`Database`/`EncryptionOption`; legacy `SocketTimeout`/`UseResultsetStreaming` jdbc properties translated/dropped; CSE_KMS results fall back to the `GetQueryResults` fetcher
* GROK-18695: Security — driver updates: databricks-jdbc 2.6.40 → 2.8.3 (embeds jackson 2.21.5, clears CVE-2026-54512/13 + log4j 2.20/netty-common 4.1.86 rows), redshift-jdbc42 2.1.0.28 → 2.2.8, mssql-jdbc 12.8.2 → 12.10.2.jre8 (clears CVE-2025-59250), logback 1.2.13 → 1.3.16, netty pins 4.1.135 → 4.1.137.Final (new 2026 netty advisories)
* GROK-18695: Security — removed the unused CData DynamoDB driver jar (provider was already unregistered) and the direct `commons-lang:2.4` dependency (`NotImplementedException` → `UnsupportedOperationException`)
* Fixed `grok_connect.sh` referencing a stale hardcoded jar version
* GROK-18695: Security — added a `<dependencyManagement>` block forcing patched, Java-8-compatible versions of vulnerable transitive dependencies: netty-bom 4.1.135.Final (transitive `netty-codec-http2`/`resolver-dns`/`redis`/…), jackson-bom 2.18.8, protobuf-java 3.25.5, commons-compress 1.26.2, commons-io 2.16.1, commons-beanutils 1.11.0, json-smart 2.5.1, plexus-utils 3.6.1, hadoop-yarn-server-common 3.3.2.
* GROK-18695: Security — second round: nimbus-jose-jwt 9.37.4 (CVE-2025-53864); new pins json-io 4.14.1 (CVE-2023-34610), xmlsec 2.3.4 (CVE-2023-44483), commons-lang3 3.18.0 (CVE-2025-48924), commons-configuration2 2.15.0 (CVE-2024-29131/29133, CVE-2026-45205); excluded `org.jline:jline` (GHSA-2r2c-cx56-8933, GHSA-47qp-hqvx-6r3f) and `hadoop-shaded-protobuf_3_7` (embeds protobuf-java 3.7.1, CVE-2024-7254) from hive-jdbc; BigQuery companion jackson jars in `lib/` bumped 2.14.3 → 2.18.8 (CVE-2026-54512/54513/54514).

# 2.7.0 (2026-07-07)

* Athena: Fixed "No suitable driver found" when connecting with temporary/STS session credentials (register JDBC driver before the direct DriverManager path)

# 2.6.5 (2026-06-22)

* Set default locale to C.UTF-8 for container build

# 2.6.4 (2026-05-18)

* Snowflake: OAuth authentication support
* Generalized lazy OAuth consent for data connectors
* GROK-20093: Per-connector OAuth flavour declared on OAuthSpec (host-suffix rules + RFC 8693 token-exchange specs; removes Databricks-specific knowledge from datlas)
* OAuth scopes moved from Credentials to DataConnection
* Databricks: Added openid scope to OAuth template
* VisualQuery: Fixed BigQuery path qualification
* Skip login/password in JDBC properties when auth method doesn't use them

# 2.6.3 (2026-05-18)

* Fixes list<string> parameters support with punctuation for Snowflake and MSSQL

# 2.6.2 (2026-03-26)

* Snowflake: OAuth (JWT) authentication with per-user identity

# 2.6.1 (2026-03-13)

* Databricks: Federated SSO auth mechanism fix

# 2.6.0 (2026-03-10)

* Database catalogs support (multi-catalog browsing for Postgres, MySQL, MS SQL, Databricks, Snowflake)
* VisualQuery: Cross-catalog joins
* SqlAnnotator: Query annotation with column metadata (types, tags)
* Athena: Exposed UseResultsetStreaming property
* Support bigint in/not in patterns
* Databricks: Improved Unity Catalog detection with fallback logic
* String column encoding speed-up
* Refactored serialization, reduced code size
* Connection pool cleanup and locking improvements
* Remediated ~20 CVEs from vulnerability scan

# 2.5.6 (2025-12-18)

* Databricks: Correction in detection of environment
* Minor fixes and improvements in TableQueries

# 2.5.5 (2025-12-02)

* Athena: Possibility to set socket timeout property in Properties

# 2.5.4 (2025-11-19)

* Databricks: Use meta tables to get info about schemas and tables instead of JDBC ConnectionMeta

# 2.5.3 (2025-11-03)

* Databricks: Oauth2(M-M) scheme support, driver update, schema browsing and visual query fixes

# 2.5.2 (2025-10-02)

* Athena: Fixed STS token usage

# 2.5.1 (2025-09-30)

* Athena: Apply session token when passed

# 2.5.0 (2025-09-19)

* TableQuery: Support joining between different schemas
* BigQuery: Support token resolution
* Postgres, Oracle, Snowflake: Additional driver properties to tune
* DynamoDB: Removed connector from the list, friendly name for ssl

# 2.4.0 (2025-07-24)

* BigQuery: Improved provider support, schem browsing, service account auth
* Fixed bug with duplicated columns in schemas during browsing

## Requires

* Datagrok >= 1.24.0

# 2.3.31 (2025-07-03)

* Snowflake: RSA keys authentication
* TableQuery: Add schema name to join
* Possibility to set port using env variables
* Minor bugs fixes

# 2.3.30 (2025-03-125)

* Memory leaks fixes

# 2.3.29 (2025-02-17)

* Postgres: Set uuid as object

# 2.3.28 (2025-02-12)

* TableQuery: Order by support and having with aggregations
* Bugs fixes and improvements

## Requires

* Datagrok >= 1.24.0

# 2.2.25 (2024-12-02)

* TableQuery: Join supports

# 2.1.20 (2024-09-17)

* Foreign keys retrieval support
* Denodo provider fixes and schema browsing support
* Dependencies bump

# 2.1.19 (2024-08-01)

* Bug fixes in parameters parsing
* Bug fixes in TableQuery
* Disables redundant logging

# 2.1.18 (2024-07-17)

* Bug fixes

# 2.1.17 (2024-04-19)

* Bug fixes

# 2.1.16 (2024-03-20)

* Fix bug with missing last rows of returned DataFrame

# 2.1.15 (2024-03-19)

* Batch mode: use `--batch` for backward compatibility

## 2.1.14 (2024-02-21)

* Increase WebSocket max size of incoming string message

## 2.1.13 (2024-01-19)

### Bug fixes

* Fix empty DataFrame with no columns when no results were returned from the query
* Fix MSSQL credentials expose in logs
* Fix TableQuery. Creates a query with aliases for column names
* Fix blob, clob and nclob display for Oracle
* Fix deadlock when several drivers have been loaded in parallel

## 2.1.12

* Make logging async
* Improve queries logging and debugging

## 2.1.10 - 2.1.11

* Fix categorize()
* Fix NullPointer when update query executed
* Fix Snowflake connection building
* Fix bit type support in columns

## 2.1.9

* Update of Neptune driver version

## 2.1.8 

This release focuses on fixing bugs

### Bug fixes

* Fix connection form of Athena Provider
* Fix exposing sensitive data in logs

## 2.1.7 (2023-07-28)

This release focuses on fixing bugs

### Features

* Add possibility to set log level in command line args when grok connect is starting

### Bug fixes

* Fix connection form of Impala Provider
* Fix datetime parameter set in UTC
* Remove connection pool for Neptune due to memory leak

## 2.1.6 (2023-07-24)

This release focuses on fixing bugs

### Bug fixes

* Fix bug with socket data sending with the size of 100 when fetch size change is not supported
