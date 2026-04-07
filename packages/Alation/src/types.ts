export type alationTokenType = 'refresh_token' | 'api_access_token';
export type alationDataType = {[key: string]: string | number};
export type refreshTokenResponse = {created_at: string, last_used_at: string | null, name: string,
  refresh_token: string, token_expires_at: string, token_status: string, user_id: number};
export type apiTokenResponse = {api_access_token: string, created_at: string, token_expires_at: string,
  token_status: string, user_id: number};
export type tokenResponse = refreshTokenResponse | apiTokenResponse;
export type dataSource = {host: string, port: number, deployment_setup_complete: boolean, db_username: string,
  dbname: string, supports_explain: boolean, has_aws_glue_metastore: boolean, supports_qli_daterange: boolean,
  icon: boolean, qli_aws_access_key_id: string, private: boolean, data_upload_disabled_message: string,
  obfuscate_literals: null, is_virtual: boolean, metastore_type: number, latest_extraction_successful: boolean,
  is_hidden: boolean, id: number, owner_ids: number[], exclude_additional_columns_in_qli: boolean,
  disable_auto_extraction: boolean, hive_logs_source_type: number, metastore_uri: null, qualified_name: string,
  is_hive: boolean, title: string, dbtype: string, supports_default_schema_extraction: boolean,
  supports_profiling_v2: boolean, qli_aws_region: string, is_gone: boolean, favorited_by_list: null,
  webhdfs_server: null, enabled_in_compose: boolean, supports_profiling: boolean, supports_qli_diagnostics: boolean,
  is_presto_hive: boolean, nosql_mde_sample_size: number, qli_hive_connection_source: number, cron_extraction: string,
  negative_filter_words: null, supports_compose: boolean, hive_logs_source: null, has_hdfs_based_qli: boolean,
  has_aws_s3_based_qli: boolean, can_data_upload: boolean, description: string, enable_designated_credential: boolean,
  all_schemas: string | null, deleted: boolean, supports_md_diagnostics: boolean, aws_region: null, limit_schemas: null,
  has_previewable_qli: boolean, hive_tez_logs_source: null, remove_filtered_schemas: boolean, profiling_tip: null,
  aws_access_key_id: null, has_metastore_uri: boolean, webhdfs_username: string, webhdfs_port: number,
  latest_extraction_time: null, uri: string, url: string, enable_default_schema_extraction: boolean, jdbc_driver: null,
  unresolved_mention_fingerprint_method: number, builtin_datasource: null, otype: string,
  can_toggle_ds_privacy: boolean, exclude_schemas: null, compose_oauth_enabled: boolean};
export type schema = {custom_fields: any[], db_comment: null | string, description: string, ds_id: number, id: number,
  key: string, name: string, title: string, url: string};
export type table = {base_table_key: null | string, description: string, ds_id: number, id: number, key: string,
  name: string, custom_fields: {value: {otype: string, oid: number}[], field_id: number, field_name: string}[],
  partition_columns: null | string, partition_definition: null | string, schema_id: number, schema_name: string,
  sql: null | string, table_comment: null | string, table_type: string, title: string, url: string};
export type column = {name: string};
export type baseEntity = dataSource | schema | table | query;
export type specialType = 'data-source' | 'schema' | 'table' | 'query';
export type createApiTokenResponse = {api_access_token: string, user_id: number, created_at: string,
  token_expires_at: string, token_status: string};
export type createRefreshTokenResponse = {refresh_token: string, user_id: number, created_at: string, name: string,
  token_expires_at: string, token_status: string};
export type query = {datasource_id: number, autosave_content: null | boolean, content: string, saved: boolean,
  published: boolean, url: string, id: number, title: string, description: string,
  datasource: {id: number, title: string, uri: string, url: string}, ts_last_saved: null | string,
  has_unsaved_changes: null | boolean, catalog_url: string, compose_url: string, schedules: any[]};
