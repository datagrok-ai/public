export interface ToolInputs {
  Read: {file_path?: string};
  Write: {file_path?: string};
  Edit: {file_path?: string};
  Bash: {command?: string};
  Glob: {pattern?: string; path?: string};
  Grep: {pattern?: string; path?: string};
  WebSearch: {query?: string};
  WebFetch: {url?: string};
}

export interface McpInputs {
  call_function: {name?: string};
  list_functions: {filter?: string};
  get_function: {id?: string};
  list_files: {path?: string};
  download_file: {path?: string};
  upload_file: {path?: string};
}

export type ToolName = keyof ToolInputs;
export type McpName = keyof McpInputs;

export interface UserMessage {
  type: 'user_message';
  sessionId: string;
  message: string;
  apiKey?: string;
  mcpServerUrl?: string;
}

export interface AbortMessage {
  type: 'abort';
  sessionId: string;
}

export type IncomingMessage = UserMessage | AbortMessage;

export type OutgoingMessage =
  | {type: 'chunk'; sessionId: string; content: string}
  | {type: 'tool_activity'; sessionId: string; summary: string}
  | {type: 'tool_result'; sessionId: string; content: string}
  | {type: 'final'; sessionId: string; content: string}
  | {type: 'error'; sessionId: string; message: string}
  | {type: 'aborted'; sessionId: string};
