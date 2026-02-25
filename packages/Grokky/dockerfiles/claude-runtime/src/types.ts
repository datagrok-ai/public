export interface UserMessage {
  type: 'user_message';
  sessionId: string;
  message: string;
}

export type OutgoingMessage =
  | {type: 'chunk'; sessionId: string; content: string}
  | {type: 'tool_activity'; sessionId: string; summary: string}
  | {type: 'tool_result'; sessionId: string; content: string}
  | {type: 'final'; sessionId: string; content: string}
  | {type: 'error'; sessionId: string; message: string};
