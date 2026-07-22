/** Task queue consumer. The worker OWNS the task queue declaration ({durable: true});
 *  datlas publishes into it with passive: true (call_queue_service.dart _putTask) and
 *  expects an 'accepted' fanout reply within taskPickupTimeout (60s), so ACCEPTED is
 *  published immediately on message receipt — datlas' own worker mode does the same
 *  (ack + ACCEPTED before parsing, _handleTasksQueueMessage). */
import * as amqplib from 'amqplib';

import {FanoutPublisher} from './fanout-publisher';
import {FuncCall, FuncCallStatus} from './func-call';
import {RevokedSet} from './pidbox';
import {Settings} from './settings';
import {logError, logInfo, logWarn} from './logger';

export class CeleryConsumer {
  static MAX_RECONNECT_WAIT_MS = 30000;

  private readonly settings: Settings;
  private readonly publisher: FanoutPublisher;
  private readonly revoked: RevokedSet;
  private readonly runTask: (call: FuncCall) => Promise<void>;
  private readonly connectFn: (url: string) => Promise<any>;
  private conn: any = null;
  private channel: any = null;
  private consumerTag: string | null = null;
  private readonly fifo: FuncCall[] = [];
  private draining = false;
  private stopped = false;
  private reconnectWaitMs = 0;

  constructor(settings: Settings, publisher: FanoutPublisher, revoked: RevokedSet,
    runTask: (call: FuncCall) => Promise<void>, connectFn?: (url: string) => Promise<any>) {
    this.settings = settings;
    this.publisher = publisher;
    this.revoked = revoked;
    this.runTask = runTask;
    this.connectFn = connectFn ?? ((url: string) => amqplib.connect(url));
  }

  get isBusy(): boolean {
    return this.draining;
  }

  async start(): Promise<void> {
    try {
      await this.connect();
      this.reconnectWaitMs = 0;
      logInfo(`Consuming task queue '${this.settings.taskQueueName}'`);
    }
    catch (e: any) {
      logError(`Task queue connection failed: ${e?.message ?? e}`);
      this.scheduleReconnect();
    }
  }

  private async connect(): Promise<void> {
    const conn = await this.connectFn(`${this.settings.brokerUrl}?heartbeat=30`);
    this.conn = conn;
    conn.on?.('error', (e: Error) => logError(`Task queue connection error: ${e.message}`));
    conn.on?.('close', () => {
      this.conn = null;
      this.channel = null;
      if (!this.stopped)
        this.scheduleReconnect();
    });
    const channel = await conn.createChannel();
    channel.on?.('error', (e: Error) => logWarn(`Task queue channel error: ${e.message}`));
    await channel.assertQueue(this.settings.taskQueueName, {durable: true});
    const consumed = await channel.consume(this.settings.taskQueueName,
      (message: any) => this.onMessage(message), {noAck: false});
    this.channel = channel;
    this.consumerTag = consumed.consumerTag;
  }

  private scheduleReconnect(): void {
    this.reconnectWaitMs = Math.min(this.reconnectWaitMs === 0 ? 1000 : this.reconnectWaitMs * 2,
      CeleryConsumer.MAX_RECONNECT_WAIT_MS);
    logInfo(`Task queue consumer reconnecting in ${this.reconnectWaitMs} ms`);
    setTimeout(() => {
      if (!this.stopped)
        void this.start();
    }, this.reconnectWaitMs);
  }

  /** Handles one task message; visible for tests (pass the channel-delivered message). */
  onMessage(message: any): void {
    if (message == null)
      return; // consumer was cancelled by the broker
    try {
      this.channel?.ack(message);
      const headers = message.properties?.headers ?? {};
      if (headers['lang'] != null && headers['lang'] !== 'js')
        logWarn(`Task '${headers['task']}' has lang='${headers['lang']}', expected 'js'`);
      const call = FuncCall.fromCeleryBody(JSON.parse(message.content.toString('utf8')));
      // ACCEPTED must reach datlas within taskPickupTimeout
      void this.publisher.publish({}, call.id, 'accepted');
      if (this.revoked.has(call.id)) {
        this.publishCanceled(call);
        return;
      }
      this.fifo.push(call);
      void this.drain();
    }
    catch (e: any) {
      logError(`Failed to handle a task message: ${e?.message ?? e}`);
      // Poison message: fail fast with an Error CALL instead of letting datlas wait
      // out the full taskPickupTimeout (60s) on a missing ACCEPTED
      const correlationId = message.properties?.correlationId;
      if (correlationId != null) {
        void this.publisher.publish({}, correlationId, 'accepted');
        void this.publisher.publish({
          id: correlationId,
          status: FuncCallStatus.ERROR,
          errorMessage: `Worker failed to parse the task message: ${e?.message ?? e}`,
        }, correlationId, 'call');
      }
    }
  }

  private publishCanceled(call: FuncCall): void {
    logInfo('Task was revoked before pickup, publishing Canceled CALL', call.id);
    call.status = FuncCallStatus.CANCELED;
    void this.publisher.publish(call.toJson(), call.id, 'call');
  }

  /** Single-in-flight executor loop draining the in-memory FIFO. */
  private async drain(): Promise<void> {
    if (this.draining)
      return;
    this.draining = true;
    try {
      while (this.fifo.length > 0) {
        const call = this.fifo.shift()!;
        if (this.revoked.has(call.id)) {
          this.publishCanceled(call);
          continue;
        }
        try {
          await this.runTask(call);
        }
        catch (e: any) {
          logError(`Task runner threw: ${e?.message ?? e}`, call.id);
        }
      }
    }
    finally {
      this.draining = false;
    }
  }

  /** Stops receiving new tasks (queued messages stay in the broker for other workers). */
  async stopConsuming(): Promise<void> {
    this.stopped = true;
    if (this.channel != null && this.consumerTag != null) {
      try {
        await this.channel.cancel(this.consumerTag);
      }
      catch (e: any) {
        logWarn(`Consumer cancel failed: ${e?.message ?? e}`);
      }
      this.consumerTag = null;
    }
  }

  /** Waits until the in-flight task (and the FIFO) is done, up to [timeoutMs]. */
  async waitForIdle(timeoutMs: number): Promise<boolean> {
    const deadline = Date.now() + timeoutMs;
    while (this.draining || this.fifo.length > 0) {
      if (Date.now() >= deadline)
        return false;
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    return true;
  }

  async close(): Promise<void> {
    this.stopped = true;
    const conn = this.conn;
    this.conn = null;
    this.channel = null;
    if (conn != null) {
      try {
        await conn.close();
      }
      catch (_) {}
    }
  }
}
