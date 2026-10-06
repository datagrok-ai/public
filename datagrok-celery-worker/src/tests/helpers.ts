import {Settings} from '../settings';

export function testSettings(maxConcurrentTasks?: number): Settings {
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

export const mockPublisher = (): any => ({publish: jest.fn().mockResolvedValue(true)});

export const callJson = (id: string): any => ({'id': id, 'func': {'name': 'f', 'params': []}, 'aux': {}});
