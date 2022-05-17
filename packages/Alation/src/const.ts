export const STORAGE_NAME = 'alation';
export const REFRESH_TOKEN_KEY = 'refresh_token';
export const API_TOKEN_KEY = 'api_access_token';
export const USER_ID = 'user_id';

export const URI_MAP = {
  refresh_token: 'validateRefreshToken',
  api_access_token: 'validateAPIAccessToken',
  datasource: 'datasource',
  schema: 'schema',
  table: 'table',
  column: 'column',
  result: 'result',
  create_api_token: 'createAPIAccessToken',
  account_auth: 'account/auth',
};

export const DATASOURCE_URI_MAP = {
  schemas: (ds_id: number) => `${ds_id}/available_schemas`,
};

export enum TOKEN_STATUS {
  ACTIVE = 'ACTIVE',
}

export const FORMAT = {
  csv: 'csv',
};

export const DATA_SOURCE_TYPES = ['Access', 'Athena', 'BigQuery', 'Cassandra', 'DB2', 'Excel', 'Firebird', 'HBase',
  'Hive', 'Hive2', 'MariaDB', 'MongoDB', 'MS SQL', 'MySQL', 'Neo4j', 'ODATA', 'ODBC', 'Oracle', 'PostgreSQL',
  'PostgresNet', 'Regshift', 'SQLite', 'Square', 'Socrata', 'Sparql', 'Terdata', 'Twitter', 'Vertica', 'Web', 'Dropbox',
  'GitHub', 'Git', 'FTP', 'GoogleCloud', 'S3', 'HTTP', 'Files', 'SFTP', 'AWS'];
