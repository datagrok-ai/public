import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';

export type ChunkEvent = {sessionId: string, content: string};
export type ToolActivityEvent = {sessionId: string, summary: string};
export type ToolResultEvent = {sessionId: string, content: string};
export type FinalEvent = {sessionId: string, content: string, structured_output?: any};
export type ErrorEvent = {sessionId: string, message: string};
export type AbortedEvent = {sessionId: string};
export type InputRequestEvent = {sessionId: string, toolName: string, input: any};

export class ClaudeRuntimeClient {
  private static instance: ClaudeRuntimeClient | null = null;
  private ws: WebSocket | null = null;
  private containerId: string | null = null;
  private mcpServerUrl: string | null = null;

  public onChunk = new rxjs.Subject<ChunkEvent>();
  public onToolActivity = new rxjs.Subject<ToolActivityEvent>();
  public onToolResult = new rxjs.Subject<ToolResultEvent>();
  public onFinal = new rxjs.Subject<FinalEvent>();
  public onError = new rxjs.Subject<ErrorEvent>();
  public onAborted = new rxjs.Subject<AbortedEvent>();
  public onInputRequest = new rxjs.Subject<InputRequestEvent>();
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
      const [runtimeContainers, mcpContainers] = await Promise.all([
        grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list(),
        grok.dapi.docker.dockerContainers.filter('name = "grokky-mcp-server"').list(),
      ]);
      if (runtimeContainers.length > 0) {
        this.containerId = runtimeContainers[0].id;
        this.ws = await grok.dapi.docker.dockerContainers.webSocketProxy(this.containerId, '/ws');
      }
      if (mcpContainers.length > 0)
        this.mcpServerUrl = `${grok.dapi.root}/docker/containers/proxy/${mcpContainers[0].id}/mcp`;
    } catch (e) {
      console.error('Failed to connect to Claude runtime:', e);
    }

    if (!this.ws)
      throw new Error('Claude runtime container is not running');

    this.ws.onmessage = (event: MessageEvent) => {
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
          ...(data.structured_output ? {structured_output: data.structured_output} : {}),
        });
        break;
      case 'error':
        this.onError.next({sessionId: data.sessionId, message: data.message});
        break;
      case 'aborted':
        this.onAborted.next({sessionId: data.sessionId});
        break;
      case 'input_request':
        this.onInputRequest.next({sessionId: data.sessionId, toolName: data.toolName, input: data.input});
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

  send(sessionId: string, message: string, options?: {outputSchema?: object}): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      throw new Error('ClaudeRuntimeClient: WebSocket is not connected');
    this.ws.send(JSON.stringify({
      type: 'user_message', sessionId, message,
      apiKey: grok.dapi.token,
      mcpServerUrl: this.mcpServerUrl,
      ...(options?.outputSchema ? {outputSchema: options.outputSchema} : {}),
    }));
  }

  abort(sessionId: string): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      return;
    this.ws.send(JSON.stringify({type: 'abort', sessionId}));
  }

  respondToInput(sessionId: string, value: any): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      return;
    this.ws.send(JSON.stringify({type: 'input_response', sessionId, value}));
  }

  async query(message: string, options?: {sessionId?: string, outputSchema?: object}): Promise<any> {
    await this.ensureConnected();
    const sid = options?.sessionId ?? `query-${Date.now()}`;
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
      this.send(sid, message, options?.outputSchema ? {outputSchema: options.outputSchema} : undefined);
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
    this.onAborted.complete();
    this.onInputRequest.complete();
    this.onClose.complete();
    this.onChunk = new rxjs.Subject<ChunkEvent>();
    this.onToolActivity = new rxjs.Subject<ToolActivityEvent>();
    this.onToolResult = new rxjs.Subject<ToolResultEvent>();
    this.onFinal = new rxjs.Subject<FinalEvent>();
    this.onError = new rxjs.Subject<ErrorEvent>();
    this.onAborted = new rxjs.Subject<AbortedEvent>();
    this.onInputRequest = new rxjs.Subject<InputRequestEvent>();
    this.onClose = new rxjs.Subject<void>();
  }
}
