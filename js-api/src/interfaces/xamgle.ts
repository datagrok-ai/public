/// this file was generated automatically from xamgle classes declarations


export interface IDetailedLogSettings {
  type: string;

  jiraTicketNumber: string;

  assigneeId: string;

  currentTime: number;

  id: string;

  isAuto: boolean;

  reportNumber: string;

  screenshot: string;

  description: string;

  reporterId: string;

  /// details
  details: Map<string, any>;

  sendEmail: boolean;

  /// reports email
  reportEmail: string;

  clientSettings: Array<Map<string, any>>;

  serverSettings: Array<Map<string, any>>;

  /// errors
  errors: Array<Map<string, any>>;

  /// client logs
  clientLog: Array<Map<string, any>>;

  /// server logs
  serverLog: Array<Map<string, any>>;

  /// console logs
  console: Array<Map<string, any>>;

  /// query logs
  grokConnectLog: Array<Map<string, any>>;

  /// scripting logs
  scriptingLog: Array<Map<string, any>>;

  /// containers logs
  containersLog: Array<Map<string, any>>;

  /// images logs
  imagesLog: Array<Map<string, any>>;

  /// services info
  services: Array<Map<string, any>>;

  data: Array<any>;

}

