/** Celery pidbox (control) consumer for cancellation.
 *
 *  Datlas publishes {method: 'revoke', arguments: {task_id, terminate, signal}} to the
 *  'celery.pidbox' exchange with routing key '<control_queue>' and contentType
 *  application/json (server_docker_func.dart cancelImpl; the exchange is declared
 *  passive there, so the WORKER owns the declaration). kombu's pidbox convention is a
 *  fanout, non-durable exchange with a per-worker auto-delete queue named
 *  'celery@<hostname>.celery.pidbox'. Note: server_docker_func.dart passes
 *  ExchangeType.DIRECT, but with passive: true the broker ignores the type — the
 *  worker-asserted fanout wins and the routing key becomes irrelevant. */
import * as amqplib from 'amqplib';

import {FanoutPublisher} from './fanout-publisher';
import {FuncCall, FuncCallStatus, isCompletedStatus} from './func-call';
import {Settings} from './settings';
import {TaskContext} from './progress';
import {logError, logInfo, logWarn} from './logger';

export const PIDBOX_EXCHANGE = 'celery.pidbox';

/** Recently revoked task ids with a TTL, so a revoke that arrives before the task is
 *  picked up from the queue still cancels it. */
export class RevokedSet {
  static DEFAULT_TTL_MS = 600000;

  private readonly ttlMs: number;
  private readonly entries: Map<string, number> = new Map();

  constructor(ttlMs: number = RevokedSet.DEFAULT_TTL_MS) {
    this.ttlMs = ttlMs;
  }

  add(taskId: string): void {
    this.prune();
    this.entries.set(taskId, Date.now() + this.ttlMs);
  }

  has(taskId: string): boolean {
    this.prune();
    return this.entries.has(taskId);
  }

  private prune(): void {
    const now = Date.now();
    for (const [id, expires] of this.entries)
      if (expires <= now)
        this.entries.delete(id);
  }
}

export interface RunningTask {
  call: FuncCall;
  context: TaskContext;
}

export class PidboxConsumer {
  static MAX_RECONNECT_WAIT_MS = 30000;

  private readonly settings: Settings;
  private readonly publisher: FanoutPublisher;
  private readonly revoked: RevokedSet;
  private readonly getCurrent: () => RunningTask | null;
  private readonly connectFn: (url: string) => Promise<any>;
  private conn: any = null;
  private stopped = false;
  private reconnectWaitMs = 0;

  constructor(settings: Settings, publisher: FanoutPublisher, revoked: RevokedSet,
    getCurrent: () => RunningTask | null, connectFn?: (url: string) => Promise<any>) {
    this.settings = settings;
    this.publisher = publisher;
    this.revoked = revoked;
    this.getCurrent = getCurrent;
    this.connectFn = connectFn ?? ((url: string) => amqplib.connect(url));
  }

  async start(): Promise<void> {
    try {
      await this.connect();
      this.reconnectWaitMs = 0;
      logInfo('Pidbox consumer connected');
    }
    catch (e: any) {
      logError(`Pidbox consumer connection failed: ${e?.message ?? e}`);
      this.scheduleReconnect();
    }
  }

  private async connect(): Promise<void> {
    const conn = await this.connectFn(`${this.settings.brokerUrl}?heartbeat=30`);
    this.conn = conn;
    conn.on?.('error', (e: Error) => logError(`Pidbox connection error: ${e.message}`));
    conn.on?.('close', () => {
      this.conn = null;
      if (!this.stopped)
        this.scheduleReconnect();
    });
    const channel = await conn.createChannel();
    channel.on?.('error', (e: Error) => logWarn(`Pidbox channel error: ${e.message}`));
    await channel.assertExchange(PIDBOX_EXCHANGE, 'fanout', {durable: false});
    const queueName = `celery@${this.settings.celeryHostname}.celery.pidbox`;
    await channel.assertQueue(queueName, {durable: false, autoDelete: true, arguments: {'x-expires': 10000}});
    await channel.bindQueue(queueName, PIDBOX_EXCHANGE, '');
    await channel.consume(queueName, (message: any) => this.onMessage(message), {noAck: true});
  }

  private scheduleReconnect(): void {
    this.reconnectWaitMs = Math.min(this.reconnectWaitMs === 0 ? 1000 : this.reconnectWaitMs * 2,
      PidboxConsumer.MAX_RECONNECT_WAIT_MS);
    logInfo(`Pidbox consumer reconnecting in ${this.reconnectWaitMs} ms`);
    setTimeout(() => {
      if (!this.stopped)
        void this.start();
    }, this.reconnectWaitMs);
  }

  private onMessage(message: any): void {
    if (message == null)
      return;
    try {
      // both datlas (contentType application/json) and kombu (json serializer)
      // put a JSON document in the message body
      this.handleControl(JSON.parse(message.content.toString('utf8')));
    }
    catch (e: any) {
      logWarn(`Ignoring malformed pidbox message: ${e?.message ?? e}`);
    }
  }

  /** Parses a control message body; visible for tests. */
  handleControl(body: any): void {
    if (body?.['method'] !== 'revoke')
      return;
    const args = body['arguments'] ?? {};
    // datlas sends a single arguments.task_id; celery broadcast revokes may carry a list
    // (task_id or task_ids)
    const ids: string[] = ([] as any[]).concat(args['task_id'] ?? args['task_ids'] ?? [])
      .filter((x) => typeof x === 'string');
    for (const id of ids)
      this.revoke(id);
  }

  private revoke(taskId: string): void {
    logInfo('Revoke received', taskId);
    this.revoked.add(taskId);
    const current = this.getCurrent();
    if (current == null || current.call.id !== taskId)
      return;
    current.context.cancelled = true;
    // python send_canceled_call parity: publish the Canceled CALL immediately and
    // ALWAYS via AMQP, even for pipe calls; the runner then suppresses the duplicate
    if (!current.context.cancelPublished && !isCompletedStatus(current.call.status)) {
      current.call.status = FuncCallStatus.CANCELED;
      current.context.cancelPublished = true;
      void this.publisher.publish(current.call.toJson(), current.call.id, 'call');
    }
  }

  async close(): Promise<void> {
    this.stopped = true;
    const conn = this.conn;
    this.conn = null;
    if (conn != null) {
      try {
        await conn.close();
      }
      catch (_) {}
    }
  }
}
