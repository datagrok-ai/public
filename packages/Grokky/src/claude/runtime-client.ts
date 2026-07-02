import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';

export const ClaudeModel = {
  Haiku: 'haiku',
  Sonnet: 'sonnet',
  Opus: 'opus',
} as const;
export type ClaudeModel = typeof ClaudeModel[keyof typeof ClaudeModel];

export type ChunkEvent = {sessionId: string, content: string};
export type ToolActivityEvent = {sessionId: string, summary: string};
export type FinalEvent = {sessionId: string, content: string, structured_output?: any, unverified?: boolean};
export type ErrorEvent = {sessionId: string, message: string};
export type AbortedEvent = {sessionId: string};
export type InputRequestEvent = {sessionId: string, requestId: string, toolName: string, input: any};
export type AuthUrlEvent = {url: string};
export type AuthErrorEvent = {message: string};

export class ClaudeRuntimeClient {
  private static instance: ClaudeRuntimeClient | null = null;
  private ws: WebSocket | null = null;
  private containerId: string | null = null;
  private mcpServerUrl: string | null = null;
  private _connectPromise: Promise<void> | null = null;


  public onChunk = new rxjs.Subject<ChunkEvent>();
  public onToolActivity = new rxjs.Subject<ToolActivityEvent>();
  public onFinal = new rxjs.Subject<FinalEvent>();
  public onError = new rxjs.Subject<ErrorEvent>();
  public onAborted = new rxjs.Subject<AbortedEvent>();
  public onInputRequest = new rxjs.Subject<InputRequestEvent>();
  public onSyncStatus = new rxjs.Subject<{status: string; message?: string; files?: string[]}>();
  public onClose = new rxjs.Subject<void>();
  public onAuthUrl = new rxjs.Subject<AuthUrlEvent>();
  public onAuthDone = new rxjs.Subject<void>();
  public onAuthError = new rxjs.Subject<AuthErrorEvent>();
  private _skillNames: string[] | null = null;

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
    if (this.connected)
      return;
    if (!this._connectPromise) {
      this._connectPromise = this.connect()
        .finally(() => { this._connectPromise = null; });
    }
    return this._connectPromise;
  }

  async connect(): Promise<void> {
    if (this.connected)
      return;

    try {
      if (!this.containerId || !this.mcpServerUrl) {
        const [runtimeContainers, mcpContainers] = await Promise.all([
          grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list(),
          grok.dapi.docker.dockerContainers.filter('name = "grokky-mcp-server"').list(),
        ]);
        this.containerId = runtimeContainers[0]?.id ?? null;
        this.mcpServerUrl = mcpContainers[0] ?
          `${grok.dapi.root}/docker/containers/proxy/${mcpContainers[0].id}/mcp` :
          null;
      }
      if (this.containerId)
        this.ws = await grok.dapi.docker.dockerContainers.webSocketProxy(this.containerId, '/ws');
    } catch (e) {
      this.containerId = null;
      this.mcpServerUrl = null;
      console.error('Failed to connect to Claude runtime:', e);
    }

    if (!this.ws)
      throw new Error('Claude runtime container is not running');

    this.ws.onmessage = (event: MessageEvent) => {
      let data: any;
      try {
        data = JSON.parse(event.data);
      } catch {
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
      case 'final':
        this.onFinal.next({
          sessionId: data.sessionId, content: data.content,
          ...(data.structured_output ? {structured_output: data.structured_output} : {}),
          ...(data.unverified ? {unverified: true} : {}),
        });
        break;
      case 'error':
        this.onError.next({sessionId: data.sessionId, message: data.message});
        break;
      case 'aborted':
        this.onAborted.next({sessionId: data.sessionId});
        break;
      case 'input_request':
        this.onInputRequest.next({sessionId: data.sessionId, requestId: data.requestId, toolName: data.toolName, input: data.input});
        break;
      case 'sync_status':
        if (data.status === 'done' && Array.isArray(data.files))
          this._skillNames = data.files.map((f: string) => f.replace(/\.[^.]+$/, ''));
        this.onSyncStatus.next({status: data.status, message: data.message, files: data.files});
        break;
      case 'auth_url':
        this.onAuthUrl.next({url: data.url});
        break;
      case 'auth_done':
        this.onAuthDone.next();
        break;
      case 'auth_error':
        this.onAuthError.next({message: data.message});
        break;
      }
    };

    this.ws.onclose = () => {
      this.ws = null;
      this.onClose.next();
    };

    this.ws.onerror = () => {
      this.ws = null;
    };
  }

  send(sessionId: string, message: string, options?: {outputSchema?: object; systemPromptMode?: string; model?: ClaudeModel}): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      throw new Error('ClaudeRuntimeClient: WebSocket is not connected');
    this.ws.send(JSON.stringify({
      type: 'user_message', sessionId, message,
      apiKey: grok.dapi.token,
      mcpServerUrl: this.mcpServerUrl,
      ...(options?.outputSchema ? {outputSchema: options.outputSchema} : {}),
      ...(options?.systemPromptMode ? {systemPromptMode: options.systemPromptMode} : {}),
      ...(options?.model ? {model: options.model} : {}),
    }));
  }

  async syncUserFiles(
    scope: 'all' | 'user-files' | 'packages' | 'shared' = 'all',
    packageName?: string,
  ): Promise<void> {
    await this.ensureConnected();
    this.ws!.send(JSON.stringify({
      type: 'sync_user_files',
      apiKey: grok.dapi.token,
      mcpServerUrl: this.mcpServerUrl,
      scope,
      ...(packageName ? {packageName} : {}),
    }));
  }

  getSkillNames(): string[] {
    return this._skillNames ?? [];
  }

  abort(sessionId: string): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      return;
    this.ws.send(JSON.stringify({type: 'abort', sessionId}));
  }

  startAuth(): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      throw new Error('ClaudeRuntimeClient: WebSocket is not connected');
    this.ws.send(JSON.stringify({type: 'auth_start'}));
  }

  sendAuthCode(code: string): void {
    this.ws?.send(JSON.stringify({type: 'auth_code', code}));
  }

  respondToInput(sessionId: string, requestId: string, value: any): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      return;
    let payload: string;
    try {
      payload = JSON.stringify({type: 'input_response', sessionId, requestId, value});
    } catch {
      // Non-serializable executed-JS return: reply with a failure both datagrok_exec and datagrok_verify can read.
      const fallback = {success: false, passed: false, error: 'result is not serializable'};
      payload = JSON.stringify({type: 'input_response', sessionId, requestId, value: fallback});
    }
    this.ws.send(payload);
  }

  async query(message: string, options?: {sessionId?: string, outputSchema?: object, model?: ClaudeModel, systemPromptMode?: string}): Promise<any> {
    await this.ensureConnected();
    const sid = options?.sessionId ?? `query-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
    return new Promise((resolve, reject) => {
      const subs: {unsubscribe: () => void}[] = [];
      const cleanup = () => { subs.forEach((s) => s.unsubscribe()); };
      subs.push(this.onFinal.subscribe((evt) => {
        if (evt.sessionId !== sid) return;
        cleanup();
        resolve(options?.outputSchema && evt.structured_output ? evt.structured_output : evt.content);
      }));
      subs.push(this.onError.subscribe((evt) => {
        if (evt.sessionId !== sid) return;
        cleanup();
        reject(new Error(evt.message));
      }));
      this.send(sid, message, {
        ...(options?.outputSchema ? {outputSchema: options.outputSchema} : {}),
        ...(options?.model ? {model: options.model} : {}),
        ...(options?.systemPromptMode ? {systemPromptMode: options.systemPromptMode} : {}),
      });
    });
  }

  disconnect(): void {
    if (this.ws) {
      this.ws.close();
      this.ws = null;
    }
    this.onChunk.complete();
    this.onToolActivity.complete();
    this.onFinal.complete();
    this.onError.complete();
    this.onAborted.complete();
    this.onInputRequest.complete();
    this.onSyncStatus.complete();
    this.onClose.complete();
    this.onAuthUrl.complete();
    this.onAuthDone.complete();
    this.onAuthError.complete();
    this.onChunk = new rxjs.Subject<ChunkEvent>();
    this.onToolActivity = new rxjs.Subject<ToolActivityEvent>();
    this.onFinal = new rxjs.Subject<FinalEvent>();
    this.onError = new rxjs.Subject<ErrorEvent>();
    this.onAborted = new rxjs.Subject<AbortedEvent>();
    this.onInputRequest = new rxjs.Subject<InputRequestEvent>();
    this.onSyncStatus = new rxjs.Subject<{status: string; message?: string}>();
    this.onClose = new rxjs.Subject<void>();
    this.onAuthUrl = new rxjs.Subject<AuthUrlEvent>();
    this.onAuthDone = new rxjs.Subject<void>();
    this.onAuthError = new rxjs.Subject<AuthErrorEvent>();
  }
}
