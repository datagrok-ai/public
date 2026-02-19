/**
 * Type definitions for entities module.
 * @module entities/types
 */

import {JoinType} from "../const";

export type PropertyGetter<TSource = any, TProp = any> = (a: TSource) => TProp;
export type PropertySetter<TSource = any, TProp = any> = (a: TSource, value: TProp) => void;
export type ValueValidator<T> = (value: T) => string;

/**
 * Represents basic properties of database connection.
 *
 * Note: The actual list of supported parameters varies by connector. For a complete list,
 * see {@link https://datagrok.ai/help/access/databases/connectors/}
 */
export interface DatabaseConnectionProperties {
  /** Database server hostname or IP address */
  server?: string;
  /** Database server port number */
  port?: number;
  /** Database name */
  db?: string;
  /** Full connection string (alternative to individual parameters) */
  connString?: string;
}

/**
 * Represents connection cache properties.
 * See also: {@link https://datagrok.ai/help/develop/how-to/function_results_cache}
 */
export interface DataConnectionCacheProperties {
  /** Whether to cache query results */
  cacheResults?: boolean;
  /** Whether to cache database schema information */
  cacheSchema?: boolean;
  /** Cron expression for cache invalidation schedule */
  cacheInvalidateSchedule?: string;
}

/**
 * Represents data connection properties.
 *
 * The available properties depend on the {@link DataSourceType}. Common properties are
 * defined explicitly, but each connector may support additional parameters.
 *
 * For a complete list of supported parameters:
 * - Database connectors: {@link https://datagrok.ai/help/access/databases/connectors/}
 * - File shares: {@link https://datagrok.ai/help/access/files/shares/}
 *
 * @example
 * // PostgreSQL connection
 * const props: DataConnectionProperties = {
 *   dataSource: 'Postgres',
 *   server: 'localhost',
 *   port: 5432,
 *   db: 'mydb',
 *   login: 'user',
 *   password: 'pass'
 * };
 *
 * @example
 * // S3 file share
 * const props: DataConnectionProperties = {
 *   dataSource: 'S3',
 *   accessKey: 'AKIA...',
 *   secretKey: '...',
 *   region: 'us-east-1',  // connector-specific parameter
 *   bucket: 'my-bucket'   // connector-specific parameter
 * };
 */
export interface DataConnectionProperties extends DatabaseConnectionProperties, DataConnectionCacheProperties {
  /** Data source type identifier (e.g., 'Postgres', 'MySQL', 'S3', 'Git') */
  dataSource: string;
  /** Login/username for authentication */
  login?: string;
  /** Password for authentication */
  password?: string;
  /** AWS access key (for AWS-based connectors) */
  accessKey?: string;
  /** AWS secret key (for AWS-based connectors) */
  secretKey?: string;
  /** AWS region (for AWS-based connectors like S3, Athena) */
  region?: string;
  /** SSL/TLS mode (for database connectors) */
  ssl?: boolean | string;
  /**
   * Additional connector-specific parameters.
   * The available parameters vary by {@link dataSource} type.
   */
  [x: string]: string | number | boolean | undefined;
}

export type FieldPredicate = {field: string, pattern: string, dataType: string};
export type FieldOrder = {field: string, asc?: boolean};
export type GroupAggregation = {aggType: string, colName: string, resultColName?: string, function?: string};
export type TableJoin = {leftTableName: string, rightTableName: string, rightTableAlias?: string, joinType: JoinType, leftTableKeys: string[], rightTableKeys: string[]};
export type DockerContainerStatus = 'stopped' | 'started' | 'pending change' | 'changing' | 'error' | 'checking';
