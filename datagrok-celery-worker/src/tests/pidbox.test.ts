import {CeleryConsumer} from '../celery-consumer';
import {FuncCall, FuncCallStatus} from '../func-call';
import {PidboxConsumer, RevokedSet, RunningTask} from '../pidbox';
import {Settings} from '../settings';
import {TaskContext} from '../progress';

function testSettings(): Settings {
  return new Settings({
    taskQueueName: 'test_queue',
    celeryHostname: 'test-host',
    celeryName: 'datagrok-celery',
    packageName: 'TestPkg',
    amqpHost: 'localhost', amqpPort: 5672, amqpUser: 'guest', amqpPassword: 'guest', amqpTls: false,
    pipeHost: 'localhost', pipePort: 3000, pipeKey: '',
    paramTimeoutMinutes: 5, wsMessageTimeoutSeconds: 30, healthPort: 8000,
  });
}

function mockPublisher(): any {
  return {publish: jest.fn().mockResolvedValue(true)};
}

function callJson(id: string): any {
  return {'id': id, 'func': {'name': 'f', 'params': []}, 'aux': {}};
}

describe('RevokedSet', () => {
  test('entries expire after the TTL', () => {
    jest.useFakeTimers();
    try {
      const revoked = new RevokedSet(1000);
      revoked.add('t1');
      expect(revoked.has('t1')).toBe(true);
      jest.advanceTimersByTime(1001);
      expect(revoked.has('t1')).toBe(false);
    }
    finally {
      jest.useRealTimers();
    }
  });
});

describe('PidboxConsumer', () => {
  test('parses the datlas revoke payload and adds to the revoked set', () => {
    const revoked = new RevokedSet();
    const pidbox = new PidboxConsumer(testSettings(), mockPublisher(), revoked, () => null);
    // exact shape from server_docker_func.dart cancelImpl
    pidbox.handleControl({'method': 'revoke', 'arguments': {'task_id': 'task-9', 'terminate': true, 'signal': 'SIGTERM'}});
    expect(revoked.has('task-9')).toBe(true);
  });

  test('accepts kombu-style task_ids lists and ignores other methods', () => {
    const revoked = new RevokedSet();
    const pidbox = new PidboxConsumer(testSettings(), mockPublisher(), revoked, () => null);
    pidbox.handleControl({'method': 'revoke', 'arguments': {'task_ids': ['a', 'b']}});
    expect(revoked.has('a')).toBe(true);
    expect(revoked.has('b')).toBe(true);
    pidbox.handleControl({'method': 'ping', 'arguments': {}});
    pidbox.handleControl(null);
    expect(revoked.has('undefined')).toBe(false);
  });

  test('revoking the running call sets its cancel flag and publishes the Canceled CALL once', () => {
    const publisher = mockPublisher();
    const call = new FuncCall(callJson('running-1'));
    const context = new TaskContext(call.id, publisher);
    const current: RunningTask = {call: call, context: context};
    const pidbox = new PidboxConsumer(testSettings(), publisher, new RevokedSet(), () => current);

    pidbox.handleControl({'method': 'revoke', 'arguments': {'task_id': 'running-1'}});
    expect(context.cancelled).toBe(true);
    expect(context.cancelPublished).toBe(true);
    expect(call.status).toBe(FuncCallStatus.CANCELED);
    expect(publisher.publish).toHaveBeenCalledTimes(1);
    const [payload, correlationId, type] = publisher.publish.mock.calls[0];
    expect(payload['status']).toBe('Canceled');
    expect(correlationId).toBe('running-1');
    expect(type).toBe('call');

    // a duplicate revoke does not publish again
    pidbox.handleControl({'method': 'revoke', 'arguments': {'task_id': 'running-1'}});
    expect(publisher.publish).toHaveBeenCalledTimes(1);
  });
});

describe('CeleryConsumer cancel-before-pickup', () => {
  test('a task revoked before pickup gets a Canceled CALL instead of running', async () => {
    const publisher = mockPublisher();
    const revoked = new RevokedSet();
    revoked.add('task-1');
    const runTask = jest.fn().mockResolvedValue(undefined);
    const consumer = new CeleryConsumer(testSettings(), publisher, revoked, runTask);

    consumer.onMessage({
      content: Buffer.from(JSON.stringify({'args': [callJson('task-1')], 'kwargs': {}}), 'utf8'),
      properties: {headers: {'lang': 'js', 'task': 'f', 'id': 'task-1'}},
    });
    await new Promise((resolve) => setImmediate(resolve));

    expect(runTask).not.toHaveBeenCalled();
    // ACCEPTED first, then the Canceled CALL
    expect(publisher.publish).toHaveBeenCalledTimes(2);
    expect(publisher.publish.mock.calls[0]).toEqual([{}, 'task-1', 'accepted']);
    const [payload, correlationId, type] = publisher.publish.mock.calls[1];
    expect(payload['status']).toBe('Canceled');
    expect(correlationId).toBe('task-1');
    expect(type).toBe('call');
  });

  test('a normal task is accepted and executed single-in-flight in FIFO order', async () => {
    const publisher = mockPublisher();
    const order: string[] = [];
    const runTask = jest.fn().mockImplementation(async (call: FuncCall) => {
      order.push(`start-${call.id}`);
      await new Promise((resolve) => setImmediate(resolve));
      order.push(`end-${call.id}`);
    });
    const consumer = new CeleryConsumer(testSettings(), publisher, new RevokedSet(), runTask);
    const messageFor = (id: string): any => ({
      content: Buffer.from(JSON.stringify({'args': [callJson(id)], 'kwargs': {}}), 'utf8'),
      properties: {headers: {}},
    });
    consumer.onMessage(messageFor('t1'));
    consumer.onMessage(messageFor('t2'));
    await consumer.waitForIdle(1000);

    expect(publisher.publish).toHaveBeenCalledWith({}, 't1', 'accepted');
    expect(publisher.publish).toHaveBeenCalledWith({}, 't2', 'accepted');
    expect(order).toEqual(['start-t1', 'end-t1', 'start-t2', 'end-t2']);
  });
});
