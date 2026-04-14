export interface ToolInputs {
  Read: {file_path?: string};
  Write: {file_path?: string};
  Edit: {file_path?: string};
  Bash: {command?: string};
  Glob: {pattern?: string; path?: string};
  Grep: {pattern?: string; path?: string};
  WebSearch: {query?: string};
  WebFetch: {url?: string};
  AskUserQuestion: {questions?: {question?: string}[]};
}

export interface McpInputs {
  call_function: {name?: string};
  list_functions: {filter?: string};
  get_function: {id?: string};
  list_files: {path?: string};
  download_file: {path?: string};
  upload_file: {path?: string};
  index_codebase: {path?: string};
  search_code: {query?: string};
  get_indexing_status: {path?: string};
  clear_index: {path?: string};
}

export type ToolName = keyof ToolInputs;
export type McpName = keyof McpInputs;

export interface UserMessage {
  type: 'user_message';
  sessionId: string;
  message: string;
  apiKey?: string;
  mcpServerUrl?: string;
  outputSchema?: object;
  systemPromptMode?: 'datagrok' | 'bash' | 'none';
}

export interface AbortMessage {
  type: 'abort';
  sessionId: string;
}

export interface InputResponseMessage {
  type: 'input_response';
  sessionId: string;
  value: any;
}

export interface SyncMessage {
  type: 'sync_user_files';
  apiKey: string;
  mcpServerUrl: string;
  scope?: string;
  packageName?: string;
}

export type IncomingMessage = UserMessage | AbortMessage | InputResponseMessage | SyncMessage;

export type OutgoingMessage =
  | {type: 'chunk'; sessionId: string; content: string}
  | {type: 'tool_activity'; sessionId: string; summary: string}
  | {type: 'tool_result'; sessionId: string; content: string}
  | {type: 'final'; sessionId: string; content: string; structured_output?: any}
  | {type: 'error'; sessionId: string; message: string}
  | {type: 'aborted'; sessionId: string}
  | {type: 'input_request'; sessionId: string; toolName: string; input: any}
  | {type: 'sync_status'; status: 'done' | 'error'; message?: string};
