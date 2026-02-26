import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';

export type ChunkEvent = {sessionId: string, content: string};
export type ToolActivityEvent = {sessionId: string, summary: string};
export type ToolResultEvent = {sessionId: string, content: string};
export type FinalEvent = {sessionId: string, content: string};
export type ErrorEvent = {sessionId: string, message: string};

const CLAUDE_RUNTIME_WS_URL = 'ws://localhost:5353/ws';

export class ClaudeRuntimeClient {
  private static instance: ClaudeRuntimeClient | null = null;
  private ws: WebSocket | null = null;
  private containerId: string | null = null;

  public onChunk = new rxjs.Subject<ChunkEvent>();
  public onToolActivity = new rxjs.Subject<ToolActivityEvent>();
  public onToolResult = new rxjs.Subject<ToolResultEvent>();
  public onFinal = new rxjs.Subject<FinalEvent>();
  public onError = new rxjs.Subject<ErrorEvent>();
  public onClose = new rxjs.Subject<void>();

  private constructor() {}

  static getInstance(): ClaudeRuntimeClient {
    if (!ClaudeRuntimeClient.instance)
      ClaudeRuntimeClient.instance = new ClaudeRuntimeClient();
    return ClaudeRuntimeClient.instance;
  }

  get connected(): boolean {
    return this.ws !== null && this.ws.readyState === WebSocket.OPEN;
  }

  async ensureConnected(): Promise<void> {
    if (!this.connected)
      await this.connect();
  }

  async connect(): Promise<void> {
    if (this.connected)
      return;

    try {
      const containers = await grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list();
      if (containers.length > 0) {
        this.containerId = containers[0].id;
        this.ws = await grok.dapi.docker.dockerContainers.webSocketProxy(
          this.containerId, '/ws', 120000,
        ) as unknown as WebSocket;
      }
    } catch (e) {
      console.error('Failed to connect to Claude runtime:', e);
    }

    if (!this.ws)
      this.ws = new WebSocket(CLAUDE_RUNTIME_WS_URL);

    await new Promise<void>((resolve, reject) => {
      const timeout = setTimeout(() => reject(new Error('Claude runtime: connection timed out')), 10000);
      this.ws!.onopen = () => { clearTimeout(timeout); resolve(); };
      this.ws!.onerror = () => { clearTimeout(timeout); reject(new Error('Claude runtime: failed to connect')); };
    });

    this.ws!.onmessage = (event: MessageEvent) => {
      let data: any;
      try {
        data = JSON.parse(event.data);
      }
      catch {
        console.error('ClaudeRuntimeClient: failed to parse message', event.data);
        return;
      }

      switch (data.type) {
      case 'chunk':
        this.onChunk.next({sessionId: data.sessionId, content: data.content});
        break;
      case 'tool_activity':
        this.onToolActivity.next({sessionId: data.sessionId, summary: data.summary});
        break;
      case 'tool_result':
        this.onToolResult.next({sessionId: data.sessionId, content: data.content});
        break;
      case 'final':
        this.onFinal.next({
          sessionId: data.sessionId, content: data.content,
        });
        break;
      case 'error':
        this.onError.next({sessionId: data.sessionId, message: data.message});
        break;
      }
    };

    this.ws!.onclose = () => {
      this.ws = null;
      this.onClose.next();
    };

    this.ws!.onerror = () => {
      this.ws = null;
    };
  }

  send(sessionId: string, message: string): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      throw new Error('ClaudeRuntimeClient: WebSocket is not connected');
    this.ws.send(JSON.stringify({type: 'user_message', sessionId, message}));
  }

  async query(message: string, sessionId?: string): Promise<string> {
    await this.ensureConnected();
    const sid = sessionId ?? `query-${Date.now()}`;
    return new Promise<string>((resolve, reject) => {
      const subs: {unsubscribe: () => void}[] = [];
      const cleanup = () => { subs.forEach((s) => s.unsubscribe()); };
      subs.push(this.onFinal.subscribe((evt) => {
        if (evt.sessionId !== sid) return;
        cleanup();
        resolve(evt.content);
      }));
      subs.push(this.onError.subscribe((evt) => {
        if (evt.sessionId !== sid) return;
        cleanup();
        reject(new Error(evt.message));
      }));
      this.send(sid, message);
    });
  }

  disconnect(): void {
    if (this.ws) {
      this.ws.close();
      this.ws = null;
    }
    this.onChunk.complete();
    this.onToolActivity.complete();
    this.onToolResult.complete();
    this.onFinal.complete();
    this.onError.complete();
    this.onClose.complete();
    this.onChunk = new rxjs.Subject<ChunkEvent>();
    this.onToolActivity = new rxjs.Subject<ToolActivityEvent>();
    this.onToolResult = new rxjs.Subject<ToolResultEvent>();
    this.onFinal = new rxjs.Subject<FinalEvent>();
    this.onError = new rxjs.Subject<ErrorEvent>();
    this.onClose = new rxjs.Subject<void>();
  }
}
