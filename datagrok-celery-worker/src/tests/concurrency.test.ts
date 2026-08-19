import {CeleryConsumer} from '../celery-consumer';
import {FuncCall} from '../func-call';
import {RevokedSet} from '../pidbox';
import {Settings} from '../settings';
import {TaskRunner} from '../task-runner';

function testSettings(maxConcurrentTasks?: number): Settings {
  return new Settings({
    taskQueueName: 'test_queue',
    celeryHostname: 'test-host',
    celeryName: 'datagrok-celery',
    packageName: 'TestPkg',
    amqpHost: 'localhost', amqpPort: 5672, amqpUser: 'guest', amqpPassword: 'guest', amqpTls: false,
    pipeHost: 'localhost', pipePort: 3000, pipeKey: '',
    paramTimeoutMinutes: 5, wsMessageTimeoutSeconds: 30, healthPort: 8000,
    maxConcurrentTasks: maxConcurrentTasks,
  });
}

const mockPublisher = (): any => ({publish: jest.fn().mockResolvedValue(true)});
const callJson = (id: string): any => ({'id': id, 'func': {'name': 'f', 'params': []}, 'aux': {}});
const messageFor = (id: string): any => ({
  content: Buffer.from(JSON.stringify({'args': [callJson(id)], 'kwargs': {}}), 'utf8'),
  properties: {headers: {}},
});
const tick = (): Promise<void> => new Promise((resolve) => setImmediate(resolve));

class FakeChannel {
  prefetchCount = Infinity;
  readonly unacked = new Set<any>();
  readonly acked: any[] = [];
  private handler: ((m: any) => void) | null = null;
  private readonly pending: any[] = [];

  on(): void {}
  async assertQueue(): Promise<void> {}
  async prefetch(n: number): Promise<void> {
    this.prefetchCount = n;
  }
  async consume(_queue: string, handler: (m: any) => void): Promise<{consumerTag: string}> {
    this.handler = handler;
    this.pump();
    return {consumerTag: 'ct'};
  }
  deliver(message: any): void {
    this.pending.push(message);
    this.pump();
  }
  ack(message: any): void {
    this.unacked.delete(message);
    this.acked.push(message);
    this.pump();
  }
  private pump(): void {
    while (this.handler != null && this.pending.length > 0 && this.unacked.size < this.prefetchCount) {
      const message = this.pending.shift();
      this.unacked.add(message);
      this.handler(message);
    }
  }
}

async function startedConsumer(max: number, runTask: jest.Mock):
  Promise<{channel: FakeChannel, consumer: CeleryConsumer}> {
  const channel = new FakeChannel();
  const connection = {on: () => {}, createChannel: async () => channel, close: async () => {}};
  const consumer = new CeleryConsumer(testSettings(max), mockPublisher(), new RevokedSet(),
    runTask, async () => connection);
  await consumer.start();
  return {channel: channel, consumer: consumer};
}

function gatedRunTask(started: string[], finish: Map<string, () => void>): jest.Mock {
  return jest.fn().mockImplementation((call: FuncCall) => {
    started.push(call.id);
    return new Promise<void>((resolve) => finish.set(call.id, resolve));
  });
}

describe('CeleryConsumer concurrency', () => {
  test('runs at most maxConcurrentTasks; acks follow completion', async () => {
    const started: string[] = [];
    const finish = new Map<string, () => void>();
    const {channel, consumer} = await startedConsumer(2, gatedRunTask(started, finish));
    for (const id of ['t1', 't2', 't3'])
      channel.deliver(messageFor(id));
    await tick();

    expect(started).toEqual(['t1', 't2']);
    expect(channel.acked.length).toBe(0);

    finish.get('t1')?.();
    await tick();
    expect(started).toEqual(['t1', 't2', 't3']);
    expect(channel.acked.length).toBe(1);

    finish.get('t2')?.();
    finish.get('t3')?.();
    expect(await consumer.waitForIdle(1000)).toBe(true);
    expect(channel.acked.length).toBe(3);
    expect(channel.unacked.size).toBe(0);
  });

  test('after a reconnect: execution stays capped, an in-flight redelivery runs once and acks late', async () => {
    const started: string[] = [];
    const finish = new Map<string, () => void>();
    const {channel, consumer} = await startedConsumer(2, gatedRunTask(started, finish));
    channel.deliver(messageFor('t1'));
    channel.deliver(messageFor('t2'));
    await tick();
    const newChannel = {ack: jest.fn()};
    consumer.onMessage(messageFor('t1'), newChannel);
    consumer.onMessage(messageFor('t3'), newChannel);
    await tick();

    expect(started).toEqual(['t1', 't2']);
    expect(newChannel.ack).not.toHaveBeenCalled();

    finish.get('t1')?.();
    await tick();
    await tick();
    expect(started).toEqual(['t1', 't2', 't3']);
    expect(newChannel.ack).toHaveBeenCalledTimes(1);

    finish.get('t2')?.();
    finish.get('t3')?.();
    expect(await consumer.waitForIdle(1000)).toBe(true);
    expect(newChannel.ack).toHaveBeenCalledTimes(2);
  });
});

describe('TaskRunner concurrency', () => {
  test('concurrent tasks report progress to their own call and clean up tracking', async () => {
    const publisher = mockPublisher();
    const gates: (() => void)[] = [];
    const impl = async () => {
      await new Promise<void>((resolve) => gates.push(resolve));
      (globalThis as any).DG_TASK_PROGRESS(50, 'halfway');
    };
    const host: any = {resolve: async () => impl, dg: null, runWithToken: (_t: any, fn: any) => fn()};
    const runner = new TaskRunner(testSettings(2), publisher, host);
    const runs = [runner.run(new FuncCall(callJson('call-a'))), runner.run(new FuncCall(callJson('call-b')))];
    await tick();
    expect(runner.runningCount).toBe(2);
    gates.forEach((open) => open());
    await Promise.all(runs);

    const progressCalls = publisher.publish.mock.calls.filter((c: any[]) => c[2] === 'progress');
    expect(progressCalls.map((c: any[]) => c[1]).sort()).toEqual(['call-a', 'call-b']);
    const resultCalls = publisher.publish.mock.calls.filter((c: any[]) => c[2] === 'call');
    expect(resultCalls.map((c: any[]) => c[0]['status'])).toEqual(['Completed', 'Completed']);
    expect(runner.runningCount).toBe(0);
  });
});
