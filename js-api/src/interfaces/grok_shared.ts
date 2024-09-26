/// this file was generated automatically from grok_shared classes declarations


export interface IUserReport {
  sendEmail: boolean;

  /// Creation date.
  createdOn: any;

  /// Data source name (such as 'Oracle'). [DataConnection.dataSource] refers to it.
  /// See also [DataSourceType].
  type: string;

  jiraTicket: string;

  isAuto: boolean;

  isResolved: boolean;

  number: number;

  /** Property description. It will be shown in the UI where possible. */
  /// Model description
  /** Property description. It will be shown in the UI where possible. */
  /// Free-text description.
  /// Description (optional).
  description: string;

  reporter: any;

  assignee: any;

  errorMessage: string;

  errorStackTrace: string;

  errorStackTraceHash: string;

  labels: Array<string>;

  /// Options.
  options: {[index: string]: any};

  data: any;

  /// Used as display name in UI.
  friendlyName: any;

  /// Name;
  /** Property name. */
  /** Property name. */
  /// File name (part of [path] after the last separator).
  name: string;

  /// Identifier
  id: string;

  securityObject: any;

  projectRelations: Array<any>;

  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  /// True is this object was deleted.
  isDeleted: boolean;

  bindId: string;

  bid: string;

  namespace: string;

  isOnServer: boolean;

  isDirty: boolean;

}

