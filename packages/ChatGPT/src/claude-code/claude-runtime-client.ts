import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {BrowserToolExecutor} from './tool-registry';

export type ChunkEvent = {sessionId: string, content: string};
export type ToolUseEvent = {sessionId: string, tool: string, input: any, status: string};
export type FinalEvent = {sessionId: string, content: string, usage?: any, cost_usd?: number, num_turns?: number};
export type ErrorEvent = {sessionId: string, message: string};

// Direct WebSocket URL to the claude-runtime container
const CLAUDE_RUNTIME_WS_URL = 'ws://localhost:5353/ws';

export class ClaudeRuntimeClient {
  private static instance: ClaudeRuntimeClient | null = null;
  private ws: WebSocket | null = null;
  // private container: DG.DockerContainer | null = null;
  private toolExecutor: BrowserToolExecutor | null = null;

  public onChunk = new rxjs.Subject<ChunkEvent>();
  public onToolUse = new rxjs.Subject<ToolUseEvent>();
  public onFinal = new rxjs.Subject<FinalEvent>();
  public onError = new rxjs.Subject<ErrorEvent>();

  private constructor() {}

  static getInstance(): ClaudeRuntimeClient {
    if (!ClaudeRuntimeClient.instance)
      ClaudeRuntimeClient.instance = new ClaudeRuntimeClient();
    return ClaudeRuntimeClient.instance;
  }

  get connected(): boolean {
    return this.ws !== null && this.ws.readyState === WebSocket.OPEN;
  }

  async connect(tableView?: DG.TableView | null, connectionId?: string | null): Promise<void> {
    if (this.connected)
      return;

    // -- Datagrok container proxy (commented out — running standalone) --
    // const containers = await grok.dapi.docker.dockerContainers.filter('name = "claude-runtime"').list();
    // this.container = containers.length > 0 ? containers[0] : null;
    // if (!this.container)
    //   throw new Error('Claude runtime container not found. Please ensure it is deployed.');
    // this.ws = await grok.dapi.docker.dockerContainers.webSocketProxy(
    //   this.container.id, '/ws', 120000
    // ) as unknown as WebSocket;

    this.toolExecutor = new BrowserToolExecutor(tableView ?? null, connectionId ?? null);
    this.ws = new WebSocket(CLAUDE_RUNTIME_WS_URL);

    await new Promise<void>((resolve, reject) => {
      const timeout = setTimeout(() => reject(new Error('Claude runtime: connection timed out')), 10000);
      this.ws!.onopen = () => { clearTimeout(timeout); resolve(); };
      this.ws!.onerror = () => { clearTimeout(timeout); reject(new Error('Claude runtime: failed to connect')); };
    });

    this.ws.onmessage = async (event: MessageEvent) => {
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
      case 'tool_use':
        this.onToolUse.next({sessionId: data.sessionId, tool: data.tool, input: data.input, status: data.status});
        break;
      case 'tool_execute': {
        const {callId, tool, args} = data;
        let result: string;
        try {
          result = await this.toolExecutor!.execute({callId, tool, args});
        }
        catch (e: any) {
          result = `Error: ${e.message}`;
        }
        this.ws!.send(JSON.stringify({type: 'tool_result', callId, result}));
        this.onToolUse.next({sessionId: data.sessionId, tool, input: args, status: 'completed'});
        break;
      }
      case 'final':
        this.onFinal.next({
          sessionId: data.sessionId, content: data.content,
          usage: data.usage, cost_usd: data.cost_usd, num_turns: data.num_turns,
        });
        break;
      case 'error':
        this.onError.next({sessionId: data.sessionId, message: data.message});
        break;
      default:
        console.warn('ClaudeRuntimeClient: unknown message type', data.type);
      }
    };

    this.ws.onclose = () => {
      console.log('ClaudeRuntimeClient: WebSocket closed');
      this.ws = null;
    };

    this.ws.onerror = (err) => {
      console.error('ClaudeRuntimeClient: WebSocket error', err);
      this.ws = null;
    };
  }

  send(sessionId: string, message: string, context?: {
    viewType?: string,
    connectionId?: string,
    repoFiles?: string[],
    metadata?: Record<string, any>,
  }): void {
    if (!this.ws || this.ws.readyState !== WebSocket.OPEN)
      throw new Error('ClaudeRuntimeClient: WebSocket is not connected');
    this.ws.send(JSON.stringify({type: 'user_message', sessionId, message, context: context ?? {}}));
  }

  updateContext(tableView?: DG.TableView | null, connectionId?: string | null): void {
    this.toolExecutor = new BrowserToolExecutor(tableView ?? null, connectionId ?? null);
  }

  disconnect(): void {
    if (this.ws) {
      this.ws.close();
      this.ws = null;
    }
    this.onChunk.complete();
    this.onToolUse.complete();
    this.onFinal.complete();
    this.onError.complete();
    // Recreate subjects so the client can be reused after disconnect
    this.onChunk = new rxjs.Subject<ChunkEvent>();
    this.onToolUse = new rxjs.Subject<ToolUseEvent>();
    this.onFinal = new rxjs.Subject<FinalEvent>();
    this.onError = new rxjs.Subject<ErrorEvent>();
  }
}
