/// this file was generated automatically from grok_shared classes declarations
import { toDart } from "../wrappers";
let api = (typeof window !== 'undefined' ? window : global.window) as any;

export class ServerMessageTypes {
  /// Chat message sent.
  static CHAT_MESSAGE_SENT = 'chat-message-sent';

  /// Notification sent to a user.
  static NOTIFICATION_SENT = 'notification-sent';

  /// Notebook table updated.
  static NOTEBOOK_TABLE_UPDATED = 'notebook-table-updated';

  /// Docker image built.
  static DOCKER_IMAGE_BUILT = 'docker-image-built';

  /// Docker container updated.
  static DOCKER_CONTAINER_UPDATED = 'docker-container-updated';

  /// Package repository error.
  static PACKAGE_REPOSITORY_ERROR = 'package-repository-error';

  /// Logger settings changed.
  static LOGGER_SETTINGS_CHANGED = 'logger-settings-changed';

  /// Package installed.
  static PACKAGE_INSTALLED = 'package-installed';

}
export class DataSourceType {
  static Access = 'Access';

  static Athena = 'Athena';

  static BigQuery = 'BigQuery';

  static Cassandra = 'Cassandra';

  static Databricks = 'Databricks';

  static DB2 = 'DB2';

  static Excel = 'Excel';

  static Firebird = 'Firebird';

  static HBase = 'HBase';

  static Hive = 'Hive';

  static Hive2 = 'Hive2';

  static MariaDB = 'MariaDB';

  static MlFlow = 'MLFlow';

  static MongoDB = 'MongoDB';

  static MsSql = 'MS SQL';

  static MySql = 'MySQL';

  static Neo4j = 'Neo4j';

  static Odata = 'ODATA';

  static Odbc = 'ODBC';

  static Oracle = 'Oracle';

  static PostgresDart = 'PostgresDart';

  static Postgres = 'Postgres';

  static Redshift = 'Redshift';

  static SQLite = 'SQLite';

  static Socrata = 'Socrata';

  static Sparql = 'Sparql';

  static Teradata = 'Teradata';

  static Twitter = 'Twitter';

  static Vertica = 'Vertica';

  static Web = 'Web';

  static AzureBlob = 'Azure Blob';

  static Dropbox = 'Dropbox';

  static GitHub = 'GitHub';

  static Ftp = 'FTP';

  static GoogleCloud = 'GoogleCloud';

  static S3 = 'S3';

  static EFS = 'Amazon EFS';

  static SharePoint = 'SharePoint';

  static Http = 'HTTP';

  static Files = 'Files';

  static Sftp = 'SFTP';

  static AWSSM = 'AWS';

  static CoreWeave = 'CoreWeave';

  static get GCP(): any { return api.grok_DataSourceType_Get_GCP(); };

  static fileDataSources = ['Azure Blob', 'Dropbox', 'Files', 'GitHub', 'GoogleCloud', 'S3', 'CoreWeave', 'SharePoint', 'Amazon EFS'];

  static systemDataSources = ['PostgresDart'];

  static secureSources = ['AWS', 'GCP'];

}
export class Permission {
  // ===== Admin =====
  static CREATE_USER = 'CreateUser';

  static EDIT_USER = 'EditUser';

  static EDIT_GROUP = 'EditGroup';

  static EDIT_GLOBAL_PERMISSIONS = 'EditGlobalPermissions';

  static START_ADMIN_SESSION = 'StartAdminSession';

  static EDIT_PLUGINS_SETTINGS = 'EditPluginsSettings';

  static PUBLISH_PACKAGE = 'PublishPackage';

  static DELETE_COMMENTS = 'DeleteComments';

  static ADMIN_SYSTEM_CONNECTIONS = 'AdminSystemConnections';

  static ADMIN_STICKY_META = 'AdminStickyMeta';

  static CREATE_REPOSITORY = 'CreateRepository';

  static CREATE_GROUP = 'CreateGroup';

  static CREATE_ROLE = 'CreateRole';

  // ===== Create =====
  static SAVE_ENTITY_TYPE = 'SaveEntityType';

  static CREATE_ENTITY = 'CreateEntity';

  static CREATE_SCRIPT = 'CreateScript';

  static CREATE_SECURITY_CONNECTION = 'CreateSecurityConnection';

  static CREATE_DATABASE_CONNECTION = 'CreateDatabaseConnection';

  static CREATE_FILE_CONNECTION = 'CreateFileConnection';

  static CREATE_DATA_QUERY = 'CreateDataQuery';

  static CREATE_DASHBOARD = 'CreateDashboard';

  static CREATE_SPACE = 'CreateSpace';

  // ===== General =====
  static INVITE_USER = 'InviteUser';

  static SHARE_WITH_EVERYONE = 'ShareWithEveryone';

  static SEND_EMAIL = 'SendEmail';

  // ===== Browse =====
  static BROWSE_FILE_CONNECTIONS = 'BrowseFileConnections';

  static BROWSE_DATABASE_CONNECTIONS = 'BrowseDatabaseConnections';

  static BROWSE_APPS = 'BrowseApps';

  static BROWSE_SPACES = 'BrowseSpaces';

  static BROWSE_DASHBOARDS = 'BrowseDashboards';

  static BROWSE_PLUGINS = 'BrowsePlugins';

  static BROWSE_FUNCTIONS = 'BrowseFunctions';

  static BROWSE_QUERIES = 'BrowseQueries';

  static BROWSE_SCRIPTS = 'BrowseScripts';

  static BROWSE_OPEN_API = 'BrowseOpenApi';

  static BROWSE_USERS = 'BrowseUsers';

  static BROWSE_GROUPS = 'BrowseGroups';

  static BROWSE_ROLES = 'BrowseRoles';

  static BROWSE_MODELS = 'BrowseModels';

  static BROWSE_DOCKERS = 'BrowseDockers';

  static BROWSE_LAYOUTS = 'BrowseLayouts';

  static BROWSE_SHARED_DATA = 'BrowseSharedData';

  // ===== Entity: Common =====
  static VIEW = 'View';

  static EDIT = 'Edit';

  static DELETE = 'Delete';

  static SHARE = 'Share';

  // ===== Entity: DataConnection =====
  static DATA_CONNECTION_QUERY = 'DataConnection.Query';

  static DATA_CONNECTION_GET_SCHEMA = 'DataConnection.GetSchema';

  static DATA_CONNECTION_LIST_FILES = 'DataConnection.ListFiles';

  // ===== Entity: DataQuery =====
  static DATA_QUERY_EXECUTE = 'DataQuery.Execute';

  // ===== Entity: DataQuery =====
  static SCRIPT_EXECUTE = 'Script.Execute';

  // ===== Entity: TableInfo =====
  static TABLE_READ_DATA = 'Table.ReadData';

}
export class ScriptLanguage {
  static Grok = 'grok';

  static Julia = 'julia';

  static Python = 'python';

  static Pyodide = 'pyodide';

  static R = 'r';

  static NodeJs = 'nodejs';

  static JavaScript = 'javascript';

  static Octave = 'octave';

  static PythonDocker = 'docker';

}
import {Entity} from '../entities'

export class DockerImage extends Entity {
  constructor(dart: any) { super(dart); };

  static fromJson(map: {[index: string]: any}): DockerImage { return new DockerImage(api.grok_DockerImage_fromJson(toDart(map)))};
  static STATUS_READY = 'ready';

  static STATUS_ERROR = 'error';

  static STATUS_BUILDING = 'building';

  static STATUS_PENDING_BUILD = 'pending build';

  static dbTableName = 'dockerfiles';

  get description(): string { return api.grok_DockerImage_Get_description(this.dart); };
  set description(x: string) {api.grok_DockerImage_Set_description(this.dart, toDart(x)); }
  get dockerfile(): string { return api.grok_DockerImage_Get_dockerfile(this.dart); };
  set dockerfile(x: string) {api.grok_DockerImage_Set_dockerfile(this.dart, toDart(x)); }
  get status(): string { return api.grok_DockerImage_Get_status(this.dart); };
  set status(x: string) {api.grok_DockerImage_Set_status(this.dart, toDart(x)); }
  get dockerName(): string { return api.grok_DockerImage_Get_dockerName(this.dart); };
  set dockerName(x: string) {api.grok_DockerImage_Set_dockerName(this.dart, toDart(x)); }
  get version(): string { return api.grok_DockerImage_Get_version(this.dart); };
  set version(x: string) {api.grok_DockerImage_Set_version(this.dart, toDart(x)); }
  get dockerfilePath(): string { return api.grok_DockerImage_Get_dockerfilePath(this.dart); };
  set dockerfilePath(x: string) {api.grok_DockerImage_Set_dockerfilePath(this.dart, toDart(x)); }
  get updatedBy(): string { return api.grok_DockerImage_Get_updatedBy(this.dart); };
  set updatedBy(x: string) {api.grok_DockerImage_Set_updatedBy(this.dart, toDart(x)); }
  get logs(): string { return api.grok_DockerImage_Get_logs(this.dart); };
  set logs(x: string) {api.grok_DockerImage_Set_logs(this.dart, toDart(x)); }
  get completed(): boolean { return api.grok_DockerImage_Get_completed(this.dart); };

  get iconStatus(): string { return api.grok_DockerImage_Get_iconStatus(this.dart); };

  get dockerFullName(): string { return api.grok_DockerImage_Get_dockerFullName(this.dart); };

}
