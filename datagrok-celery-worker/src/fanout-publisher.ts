/** Publisher for the 'calls_fanout' exchange (mirrors amqp_publisher.py).
 *  Datlas consumes this exchange in call_queue_service.dart and expects message
 *  properties correlationId = call id and type in 'accepted'|'call'|'progress'|'log'. */
import * as amqplib from 'amqplib';

import {logError, logWarn} from './logger';

export type FanoutType = 'call' | 'log' | 'accepted' | 'progress';

export const CALLS_FANOUT = 'calls_fanout';

function sleep(ms: number): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

export class FanoutPublisher {
  static MAX_ATTEMPTS = 3;

  private readonly brokerUrl: string;
  private readonly connectFn: (url: string) => Promise<any>;
  private conn: any = null;
  private channelPromise: Promise<any> | null = null;

  constructor(brokerUrl: string, connectFn?: (url: string) => Promise<any>) {
    this.brokerUrl = brokerUrl;
    this.connectFn = connectFn ?? ((url: string) => amqplib.connect(url));
  }

  private ensureChannel(): Promise<any> {
    if (this.channelPromise == null) {
      const promise = this.openChannel();
      promise.catch(() => {
        if (this.channelPromise === promise)
          this.channelPromise = null;
      });
      this.channelPromise = promise;
    }
    return this.channelPromise;
  }

  private async openChannel(): Promise<any> {
    const conn = await this.connectFn(this.brokerUrl);
    this.conn = conn;
    conn.on?.('error', (e: Error) => {
      logWarn(`Fanout publisher connection error: ${e.message}`);
      this.reset(conn);
    });
    conn.on?.('close', () => this.reset(conn));
    const channel = await conn.createChannel();
    channel.on?.('error', (e: Error) => logWarn(`Fanout publisher channel error: ${e.message}`));
    // datlas binds a private queue to this exchange with default (non-durable) settings
    await channel.assertExchange(CALLS_FANOUT, 'fanout', {durable: false});
    return channel;
  }

  private reset(conn?: any): void {
    if (conn != null && conn !== this.conn)
      return;
    this.channelPromise = null;
    const current = this.conn;
    this.conn = null;
    if (current != null) {
      try {
        current.close().catch(() => {});
      }
      catch (_) {}
    }
  }

  /** Publishes a JSON payload; 3 attempts with 2^attempt seconds backoff, reconnecting
   *  between attempts. Never throws — returns false when all attempts failed. */
  async publish(payload: object, correlationId: string, type: FanoutType): Promise<boolean> {
    for (let attempt = 0; attempt < FanoutPublisher.MAX_ATTEMPTS; attempt++) {
      try {
        const channel = await this.ensureChannel();
        channel.publish(CALLS_FANOUT, '', Buffer.from(JSON.stringify(payload)), {
          contentType: 'application/json',
          correlationId: correlationId,
          type: type,
        });
        return true;
      }
      catch (e: any) {
        this.reset();
        logWarn(`Fanout publish '${type}' failed (attempt ${attempt + 1}): ${e?.message ?? e}`, correlationId);
        if (attempt < FanoutPublisher.MAX_ATTEMPTS - 1)
          await sleep(Math.pow(2, attempt + 1) * 1000);
      }
    }
    logError(`Fanout publish '${type}' failed after ${FanoutPublisher.MAX_ATTEMPTS} attempts`, correlationId);
    return false;
  }

  async close(): Promise<void> {
    const conn = this.conn;
    this.channelPromise = null;
    this.conn = null;
    if (conn != null) {
      try {
        await conn.close();
      }
      catch (_) {}
    }
  }
}
