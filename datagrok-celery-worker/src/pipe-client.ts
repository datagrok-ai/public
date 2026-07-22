/** WebSocket client for grok_pipe implementing the param streaming protocol.
 *  Ground truth for the framing:
 *  - core/shared/grok_shared/lib/src/http_client/sockets.dart (Const, SocketSender,
 *    SocketReceiver): 'SENDING DATAFRAME <size> <tagsJson>' + binary chunks, each
 *    acknowledged with 'PART OK'; abort constant Const.ERROR = 'PART ERROR'.
 *  - datagrok-celery-task/datagrok_task.py (_get_param_grok_pipe/_send_param_grok_pipe):
 *    the python lib compares the bare string 'ERROR' instead — we accept BOTH as abort.
 *  - grok_pipe/src/server.ts: rooms are joined at ws://host:port/<callId> with headers
 *    'x-member-name' and 'authorization'; all frames are relayed to the other members. */
import WebSocket from 'ws';

import {logInfo, logWarn} from './logger';

export const Const = {
  PART_OK: 'PART OK',
  /** Dart abort constant (sockets.dart Const.ERROR). */
  PART_ERROR: 'PART ERROR',
  /** The python lib's abort comparison value (datagrok_task.py) — tolerated too. */
  ERROR: 'ERROR',
  CANCEL_PARAM: 'CANCEL PARAM',
  SENDING: 'SENDING DATAFRAME',
  PARAM: 'PARAM',
  PARAM_SENT: 'PARAM_SENT',
  CALL: 'CALL',
  PROGRESS: 'PROGRESS',
} as const;

type PipeMessage = string | Uint8Array;

interface Waiter {
  resolve: (message: PipeMessage) => void;
  reject: (error: Error) => void;
  timer: NodeJS.Timeout;
}

export interface PipeClientOptions {
  /** Per-message wait timeout, ms (DATAGROK_WS_MESSAGE_TIMEOUT). */
  messageTimeoutMs: number;
  /** Total budget for receiving one param, ms (DATAGROK_PARAM_TIMEOUT). */
  paramTimeoutMs: number;
  /** Outbound chunk size (call.aux.batchSize, default 2048000). */
  batchSize: number;
  connectAttempts?: number;
  connectRetryMs?: number;
}

export class PipeClient {
  private readonly url: string;
  private readonly headers: {[key: string]: string};
  private readonly options: PipeClientOptions;
  private ws: WebSocket | null = null;
  private inbound: PipeMessage[] = [];
  private waiter: Waiter | null = null;
  private closed = false;

  constructor(url: string, headers: {[key: string]: string}, options: PipeClientOptions) {
    this.url = url;
    this.headers = headers;
    this.options = options;
  }

  get isConnected(): boolean {
    return this.ws != null && this.ws.readyState === WebSocket.OPEN;
  }

  async connect(): Promise<void> {
    const attempts = this.options.connectAttempts ?? 5;
    const retryMs = this.options.connectRetryMs ?? 3000;
    let lastError: any = null;
    for (let attempt = 0; attempt < attempts; attempt++) {
      try {
        await this.connectOnce();
        return;
      }
      catch (e: any) {
        lastError = e;
        logWarn(`grok_pipe connection failed (attempt ${attempt + 1}/${attempts}): ${e?.message ?? e}`);
        if (attempt < attempts - 1)
          await new Promise((resolve) => setTimeout(resolve, retryMs));
      }
    }
    throw new Error(`Failed to connect to grok_pipe at ${this.url} after ${attempts} retries: ${lastError?.message ?? lastError}`);
  }

  private connectOnce(): Promise<void> {
    return new Promise((resolve, reject) => {
      const ws = new WebSocket(this.url, {headers: this.headers});
      let settled = false;
      ws.on('open', () => {
        settled = true;
        this.ws = ws;
        this.attach(ws);
        resolve();
      });
      ws.on('error', (e: Error) => {
        if (!settled) {
          settled = true;
          reject(e);
        }
      });
    });
  }

  private attach(ws: WebSocket): void {
    ws.on('message', (data: WebSocket.RawData, isBinary: boolean) => {
      const buffer = Array.isArray(data) ? Buffer.concat(data) :
        data instanceof ArrayBuffer ? Buffer.from(data) : data;
      this.push(isBinary ? new Uint8Array(buffer) : buffer.toString('utf8'));
    });
    ws.on('error', (e: Error) => logWarn(`grok_pipe socket error: ${e.message}`));
    ws.on('close', () => {
      if (this.ws === ws)
        this.ws = null;
      const waiter = this.waiter;
      if (waiter != null) {
        this.waiter = null;
        clearTimeout(waiter.timer);
        waiter.reject(new Error('grok_pipe connection was closed'));
      }
    });
  }

  private push(message: PipeMessage): void {
    const waiter = this.waiter;
    if (waiter != null) {
      this.waiter = null;
      clearTimeout(waiter.timer);
      waiter.resolve(message);
    }
    else
      this.inbound.push(message);
  }

  /** Next inbound message from the serialized queue, bounded by [timeoutMs]. */
  private nextMessage(timeoutMs: number): Promise<PipeMessage> {
    if (this.inbound.length > 0)
      return Promise.resolve(this.inbound.shift()!);
    if (!this.isConnected)
      return Promise.reject(new Error('grok_pipe connection is closed'));
    return new Promise((resolve, reject) => {
      const timer = setTimeout(() => {
        this.waiter = null;
        reject(new Error(`Timeout (${timeoutMs} ms) waiting for a message from grok_pipe`));
      }, timeoutMs);
      this.waiter = {resolve: resolve, reject: reject, timer: timer};
    });
  }

  sendText(message: string): void {
    if (!this.isConnected)
      throw new Error('grok_pipe connection is closed');
    this.ws!.send(message);
  }

  private sendBinary(data: Uint8Array): void {
    if (!this.isConnected)
      throw new Error('grok_pipe connection is closed');
    this.ws!.send(data, {binary: true});
  }

  /** Requests an input param ('PARAM <name>') and collects the peer's
   *  'SENDING DATAFRAME <size> <tags>' + binary chunks (each acked with 'PART OK'),
   *  terminated by 'PARAM_SENT <name>'. Returns null when no data was sent. */
  async receiveParam(name: string): Promise<Uint8Array | null> {
    logInfo(`Receiving param ${name}`);
    const start = Date.now();
    this.sendText(`${Const.PARAM} ${name}`);
    const chunks: Uint8Array[] = [];
    let expectedSize: number | null = null;
    let received = 0;
    for (;;) {
      if (Date.now() - start > this.options.paramTimeoutMs)
        throw new Error(`Timeout receiving param ${name}`);
      const message = await this.nextMessage(this.options.messageTimeoutMs);
      if (typeof message === 'string') {
        if (message.startsWith(`${Const.PARAM_SENT} ${name}`))
          break;
        else if (message.startsWith('SENDING'))
          expectedSize = PipeClient.parseExpectedSize(message);
        else if (message === Const.PART_ERROR || message === Const.ERROR)
          throw new Error(`Peer reported an error while sending param ${name}`);
        // unrelated frames (LOG/PROGRESS/relay echoes) are tolerated and skipped
      }
      else {
        if (expectedSize == null)
          throw new Error('Protocol error: binary data received before the SENDING header');
        chunks.push(message);
        received += message.length;
        this.sendText(Const.PART_OK);
        if (received > expectedSize)
          throw new Error(`Received more binary data than expected for param ${name}`);
      }
    }
    if (chunks.length === 0)
      return null;
    logInfo(`Received param ${name} (${received} bytes in ${chunks.length} chunks)`);
    return chunks.length === 1 ? chunks[0] : PipeClient.concat(chunks, received);
  }

  /** Sends an output param: 'SENDING DATAFRAME <size> <tagsJson>' followed by binary
   *  chunks, waiting for 'PART OK' after each one (Dart SocketSender._send parity;
   *  'CANCEL PARAM' stops sending, 'PART ERROR'/'ERROR' aborts). */
  async sendParam(data: Uint8Array, tags: {[key: string]: string}): Promise<void> {
    this.sendText(`${Const.SENDING} ${data.length} ${JSON.stringify(tags)}`);
    const batchSize = this.options.batchSize;
    for (let offset = 0; offset < data.length; offset += batchSize) {
      this.sendBinary(data.subarray(offset, Math.min(offset + batchSize, data.length)));
      for (;;) {
        const response = await this.nextMessage(this.options.messageTimeoutMs);
        if (typeof response !== 'string')
          continue;
        if (response === Const.PART_OK)
          break;
        if (response === Const.CANCEL_PARAM)
          return;
        if (response === Const.PART_ERROR || response === Const.ERROR)
          throw new Error('Peer reported an error while receiving the parameter');
        // skip unrelated text frames, keep waiting for the ack
      }
    }
  }

  /** Safe to call twice. */
  close(): void {
    if (this.closed)
      return;
    this.closed = true;
    const ws = this.ws;
    this.ws = null;
    if (ws != null) {
      try {
        ws.close();
      }
      catch (_) {}
    }
    const waiter = this.waiter;
    if (waiter != null) {
      this.waiter = null;
      clearTimeout(waiter.timer);
      waiter.reject(new Error('grok_pipe client was closed'));
    }
  }

  private static parseExpectedSize(message: string): number {
    const size = parseInt(message.split(' ')[2], 10);
    if (!Number.isFinite(size))
      throw new Error(`Could not parse size from SENDING message: ${message.substring(0, 50)}`);
    return size;
  }

  private static concat(chunks: Uint8Array[], totalLength: number): Uint8Array {
    const result = new Uint8Array(totalLength);
    let offset = 0;
    for (const chunk of chunks) {
      result.set(chunk, offset);
      offset += chunk.length;
    }
    return result;
  }
}
