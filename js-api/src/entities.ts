// noinspection JSUnusedGlobalSymbols

/**
 * Platform entity classes representing server-side objects.
 *
 * All entities extend the base {@link Entity} class and can be persisted,
 * queried, and managed through the Data API ({@link grok.dapi}).
 *
 * ## Core Entities
 *
 * - **{@link Entity}** - Base class for all database-stored objects
 * - **{@link User}** - Platform user
 * - **{@link UserSession}** - Active user session
 * - **{@link Group}** - User group for permissions
 * - **{@link Project}** - Container for tables, views, and connections
 *
 * ## Data Access
 *
 * - **{@link DataConnection}** - Database/file share connection
 * - **{@link DataQuery}** - Parameterized SQL query
 * - **{@link TableQuery}** - Query with table results
 * - **{@link TableQueryBuilder}** - Fluent API for building queries
 * - **{@link DataJob}** - Scheduled data job
 * - **{@link Credentials}** - Connection credentials
 *
 * ## Functions and Scripts
 *
 * - **{@link Func}** - Base class for executable functions
 * - **{@link Script}** - User-defined script (Python, R, JS, etc.)
 * - **{@link ScriptEnvironment}** - Script execution environment
 * - **{@link Package}** - Datagrok package
 *
 * ## Tables and Files
 *
 * - **{@link TableInfo}** - Table metadata
 * - **{@link ColumnInfo}** - Column metadata
 * - **{@link FileInfo}** - File/folder in file shares
 *
 * ## Views and Layouts
 *
 * - **{@link ViewLayout}** - Saved view layout
 * - **{@link ViewInfo}** - View metadata
 *
 * ## Logging and Reports
 *
 * - **{@link LogEvent}** / **{@link LogEventType}** - Audit log entries
 * - **{@link UserReport}** / **{@link UserReportsRule}** - User reports
 * - **{@link UserNotification}** - User notifications
 *
 * ## Other
 *
 * - **{@link Model}** - Machine learning model
 * - **{@link Notebook}** - Jupyter notebook
 * - **{@link Property}** - Entity property definition
 * - **{@link DockerContainer}** - Docker container
 * - **{@link ProgressIndicator}** - Progress tracking
 *
 * ## Connection Interfaces
 *
 * - **{@link DataConnectionProperties}** - Connection configuration
 * - **{@link DatabaseConnectionProperties}** - Database-specific config
 * - **{@link DataConnectionCacheProperties}** - Caching config
 *
 * @module entities
 */

export * from './entities/index';
