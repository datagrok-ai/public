// Browser ↔ runtime WebSocket protocol

export type PromptMode = 'datagrok' | 'bash' | 'none';

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
  systemPromptMode?: PromptMode;
  clientTools?: ClientToolDef[];
  /** Per-turn gate switches, for A/B experiments (benchmark arms, dev harness). Both default
   * true; the gates are quality mechanisms, not security boundaries, so a client may turn its
   * own off. */
  gates?: {grounding?: boolean; verify?: boolean};
  resumeExpected?: boolean;
  /** Id of the queued admission task: the turn holds until the matching
   * `Grokky:aiChatTurnTask` celery task claims it (see tasks.ts). */
  taskId?: string;
}

export interface AbortMessage {
  type: 'abort';
  sessionId: string;
}

export interface InputResponseMessage {
  type: 'input_response';
  sessionId: string;
  requestId?: string;
  value: unknown;
}

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
  // `name` is the bare tool name (mcp prefix stripped), so a harness can assert which path a turn took.
  | {type: 'tool_activity'; sessionId: string; summary: string; name: string}
  // A gate (verifier / grounding) blocked the turn's Stop and a revision is being generated.
  // The visible answer stays; the revision streams hidden, and `final.revision` says whether it
  // replaces the original ('replaced') or the original stands ('kept').
  | {type: 'revision_start'; sessionId: string}
  | {type: 'final'; sessionId: string; content: string; structured_output?: unknown; unverified?: boolean;
      metrics?: TurnMetrics; revision?: 'kept' | 'replaced'}
  | {type: 'error'; sessionId: string; message: string}
  | {type: 'session_reset'; sessionId: string}
  | {type: 'queued'; sessionId: string}
  | {type: 'aborted'; sessionId: string}
  | {type: 'input_request'; sessionId: string; requestId: string; toolName: string; input: unknown}
  | {type: 'sync_status'; status: 'done' | 'error'; message?: string; files?: string[]};

// One agent turn

/** A quality gate (verifier, grounding), fed every finished tool call of the turn and then its stop. */
export interface TurnGate {
  afterTool(toolName: string, args: Record<string, unknown>, result: string): void;
  /** Feedback that continues the turn as a hidden revision, or undefined to let the turn end. */
  atStop(): string | undefined;
}

export interface AgentTurnRequest {
  sessionKey: string;
  text: string;
  images?: ImageAttachment[];
  systemPrompt: string;
  /** Reasoning effort; without one the agent asks for no thinking where the model allows it. */
  effort?: 'high';
  cwd: string;
  userDir?: string;
  /** The browser holds a transcript for this session: if the runtime has none, the turn resets instead. */
  resumeExpected?: boolean;
  outputSchema?: object;
  clientTools?: ClientToolDef[];
  /** The Datagrok MCP server, reached with the user's key in these headers. */
  datagrokMcp?: {url: string; headers: Record<string, string>};
  limits: {maxTurns: number; maxBudgetUsd: number};
  gates: TurnGate[];
  /** Round-trips a tool call to the browser; rejects on timeout or abort. */
  awaitBrowser: (toolName: string, input: unknown, timeoutMs: number) => Promise<unknown>;
  /** Aborts the turn: the agent is cancelled, nothing more is sent to it, and the event stream ends. */
  signal: AbortSignal;
}

export type AgentEvent =
  | {type: 'text'; text: string}
  | {type: 'tool_activity'; name: string; summary: string}
  | {type: 'activity'}
  | {type: 'result'; text: string; structuredOutput?: unknown; metrics: TurnMetrics}
  | {type: 'limit'; kind: 'turns' | 'budget'}
  /** The session's stored history could not be restored; the browser replays its transcript. */
  | {type: 'session_reset'}
  | {type: 'error'; message: string};

// Runtime ↔ agent process protocol (agent/grokky_agent), one JSON object per line

export type AgentCommand =
  | {type: 'init'; credentials: Record<string, string>;
      settings: {homes: string; skillsDir: string; compressionThresholdTokens: number; warm: boolean}}
  | ({type: 'turn';
      /** The stored Hermes session a reloaded process continues. */
      restore?: string} &
      Omit<AgentTurnRequest, 'sessionKey' | 'resumeExpected' | 'gates' | 'awaitBrowser' | 'signal'>)
  | {type: 'continue'; text: string}
  | {type: 'cancel'}
  | {type: 'browser_result'; id: string; value?: unknown; error?: string};

export type AgentProcessEvent =
  | {type: 'ready'}
  | {type: 'checked'; error: string | null}
  | {type: 'fatal'; message: string}
  | {type: 'session'; id: string}
  | {type: 'text'; text: string}
  | {type: 'activity'}
  | {type: 'tool_start'; name: string; args: Record<string, unknown>}
  | {type: 'tool_end'; name: string; args: Record<string, unknown>; result: string}
  | {type: 'browser_call'; id: string; name: string; input: unknown; timeoutMs: number}
  | {type: 'stopped'; text: string; structuredOutput?: unknown; metrics: TurnMetrics}
  | {type: 'limit'; kind: 'turns' | 'budget'}
  | {type: 'error'; message: string}
  | {type: 'cancelled'}
  | {type: 'session_lost'};
