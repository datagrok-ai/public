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
  private channel: any = null;

  constructor(brokerUrl: string, connectFn?: (url: string) => Promise<any>) {
    this.brokerUrl = brokerUrl;
    this.connectFn = connectFn ?? ((url: string) => amqplib.connect(url));
  }

  private async ensureChannel(): Promise<any> {
    if (this.channel != null)
      return this.channel;
    this.conn = await this.connectFn(this.brokerUrl);
    this.conn.on?.('error', (e: Error) => {
      logWarn(`Fanout publisher connection error: ${e.message}`);
      this.reset();
    });
    this.conn.on?.('close', () => this.reset());
    const channel = await this.conn.createChannel();
    channel.on?.('error', (e: Error) => logWarn(`Fanout publisher channel error: ${e.message}`));
    // datlas binds a private queue to this exchange with default (non-durable) settings
    await channel.assertExchange(CALLS_FANOUT, 'fanout', {durable: false});
    this.channel = channel;
    return channel;
  }

  private reset(): void {
    this.channel = null;
    const conn = this.conn;
    this.conn = null;
    if (conn != null) {
      try {
        conn.close().catch(() => {});
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
    this.channel = null;
    this.conn = null;
    if (conn != null) {
      try {
        await conn.close();
      }
      catch (_) {}
    }
  }
}
