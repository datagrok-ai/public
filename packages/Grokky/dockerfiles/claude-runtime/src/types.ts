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
  datagrok_exec: {code?: string};
  datagrok_verify: {assertion?: string; description?: string};
  datagrok_show_entities: {entities?: any[]};
}

export type ToolName = keyof ToolInputs;
export type McpName = keyof McpInputs;

export const ClaudeModel = {
  Haiku: 'haiku',
  Sonnet: 'sonnet',
  Opus: 'opus',
} as const;
export type ClaudeModel = typeof ClaudeModel[keyof typeof ClaudeModel];

export interface ImageAttachment {
  mediaType: 'image/jpeg' | 'image/png' | 'image/gif' | 'image/webp';
  data: string;
}

/** A tool declared by the browser for this turn (usually supplied by the current view).
 * Calls round-trip to the browser via input_request / input_response. */
export interface ClientToolDef {
  name: string;
  description: string;
  /** JSON Schema (object type) describing the arguments. */
  inputSchema?: object;
}

export interface UserMessage {
  type: 'user_message';
  sessionId: string;
  message: string;
  images?: ImageAttachment[];
  apiKey?: string;
  mcpServerUrl?: string;
  outputSchema?: object;
  systemPromptMode?: 'datagrok' | 'bash' | 'none';
  model?: ClaudeModel;
  clientTools?: ClientToolDef[];
}

export interface AbortMessage {
  type: 'abort';
  sessionId: string;
}

export interface InputResponseMessage {
  type: 'input_response';
  sessionId: string;
  requestId?: string;
  value: any;
}

export interface SyncMessage {
  type: 'sync_user_files';
  apiKey: string;
  mcpServerUrl: string;
  scope?: string;
  packageName?: string;
}

export interface AuthStartMessage {type: 'auth_start'}
export interface AuthCodeMessage {type: 'auth_code'; code: string}

export type IncomingMessage = UserMessage | AbortMessage | InputResponseMessage | SyncMessage | AuthStartMessage | AuthCodeMessage;

/** Per-turn metrics forwarded from the SDK `result` message (see docs/BENCHMARK.md). */
export interface TurnMetrics {
  inputTokens: number | null;
  outputTokens: number | null;
  cacheReadTokens: number | null;
  cacheCreationTokens: number | null;
  costUsd: number | null;
  numTurns: number | null;
  durationMs: number | null;
  durationApiMs: number | null;
}

export type OutgoingMessage =
  | {type: 'chunk'; sessionId: string; content: string}
  | {type: 'tool_activity'; sessionId: string; summary: string}
  // A gate (verifier / grounding) blocked the turn's Stop and a revision is being generated.
  // The visible answer stays; the revision streams hidden, and `final.revision` says whether it
  // replaces the original ('replaced') or the original stands ('kept').
  | {type: 'revision_start'; sessionId: string}
  | {type: 'final'; sessionId: string; content: string; structured_output?: any; unverified?: boolean; metrics?: TurnMetrics; revision?: 'kept' | 'replaced'}
  | {type: 'error'; sessionId: string; message: string}
  | {type: 'queued'; sessionId: string}
  | {type: 'aborted'; sessionId: string}
  | {type: 'input_request'; sessionId: string; requestId: string; toolName: string; input: any}
  | {type: 'sync_status'; status: 'done' | 'error'; message?: string; files?: string[]}
  | {type: 'auth_url'; url: string}
  | {type: 'auth_done'}
  | {type: 'auth_error'; message: string};
