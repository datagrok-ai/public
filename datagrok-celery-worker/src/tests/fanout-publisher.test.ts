import {CALLS_FANOUT, FanoutPublisher} from '../fanout-publisher';

function mockConnection(publishFn: jest.Mock): any {
  return {
    on: jest.fn(),
    close: jest.fn().mockResolvedValue(undefined),
    createChannel: jest.fn().mockResolvedValue({
      on: jest.fn(),
      assertExchange: jest.fn().mockResolvedValue(undefined),
      publish: publishFn,
    }),
  };
}

describe('FanoutPublisher', () => {
  beforeEach(() => jest.useFakeTimers());
  afterEach(() => jest.useRealTimers());

  test('retries with 2^attempt backoff and succeeds on the third attempt', async () => {
    const publishFn = jest.fn().mockReturnValue(true);
    const connectFn = jest.fn()
      .mockRejectedValueOnce(new Error('conn refused 1'))
      .mockRejectedValueOnce(new Error('conn refused 2'))
      .mockResolvedValue(mockConnection(publishFn));
    const publisher = new FanoutPublisher('amqp://guest:guest@localhost:5672', connectFn);

    const resultPromise = publisher.publish({'percent': 50, 'description': 'half'}, 'call-1', 'progress');
    // attempt 1 fails -> 2s backoff, attempt 2 fails -> 4s backoff, attempt 3 succeeds
    await jest.advanceTimersByTimeAsync(2000);
    expect(connectFn).toHaveBeenCalledTimes(2);
    await jest.advanceTimersByTimeAsync(4000);
    const result = await resultPromise;

    expect(result).toBe(true);
    expect(connectFn).toHaveBeenCalledTimes(3);
    expect(publishFn).toHaveBeenCalledTimes(1);
    const [exchange, routingKey, content, options] = publishFn.mock.calls[0];
    expect(exchange).toBe(CALLS_FANOUT);
    expect(routingKey).toBe('');
    expect(JSON.parse(content.toString('utf8'))).toEqual({'percent': 50, 'description': 'half'});
    expect(options).toEqual({contentType: 'application/json', correlationId: 'call-1', type: 'progress'});
  });

  test('returns false (never throws) after all attempts fail', async () => {
    const connectFn = jest.fn().mockRejectedValue(new Error('down'));
    const publisher = new FanoutPublisher('amqp://guest:guest@localhost:5672', connectFn);
    const resultPromise = publisher.publish({}, 'call-2', 'accepted');
    await jest.advanceTimersByTimeAsync(2000);
    await jest.advanceTimersByTimeAsync(4000);
    await expect(resultPromise).resolves.toBe(false);
    expect(connectFn).toHaveBeenCalledTimes(3);
  });

  test('reuses the channel across publishes and asserts the exchange non-durable', async () => {
    const publishFn = jest.fn().mockReturnValue(true);
    const conn = mockConnection(publishFn);
    const connectFn = jest.fn().mockResolvedValue(conn);
    const publisher = new FanoutPublisher('amqp://guest:guest@localhost:5672', connectFn);
    await publisher.publish({}, 'c1', 'accepted');
    await publisher.publish({}, 'c1', 'call');
    expect(connectFn).toHaveBeenCalledTimes(1);
    const channel = await conn.createChannel.mock.results[0].value;
    expect(channel.assertExchange).toHaveBeenCalledWith(CALLS_FANOUT, 'fanout', {durable: false});
    expect(publishFn).toHaveBeenCalledTimes(2);
  });
});
