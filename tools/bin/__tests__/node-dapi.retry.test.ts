import {describe, it, expect, afterEach, beforeEach, vi} from 'vitest';
import {NodeApiClient} from '../utils/node-dapi';

const URL = 'http://stand/api';
const ok = (body: any) => new Response(JSON.stringify(body), {status: 200, headers: {'content-type': 'application/json'}});
const shed = (status: number) => new Response('<html>503</html>', {status, headers: {'content-type': 'text/html'}});

function stub(...responses: (() => Response | Promise<Response>)[]) {
  let call = 0;
  const spy = vi.fn(async () => responses[Math.min(call++, responses.length - 1)]());
  vi.stubGlobal('fetch', spy);
  return spy;
}

// The backoff is read per call, so the suite does not have to sleep through it.
beforeEach(() => { process.env['GROK_HTTP_BACKOFF'] = '0'; });
afterEach(() => { delete process.env['GROK_HTTP_BACKOFF']; vi.unstubAllGlobals(); });

describe('fetchOrRetry', () => {
  it('retries a load-shedding status and returns the answer that follows', async () => {
    const spy = stub(() => shed(503), () => shed(503), () => ok({id: 'x'}));
    const client = new NodeApiClient(URL, 'token');
    expect(await client.get('/projects/x')).toEqual({id: 'x'});
    expect(spy).toHaveBeenCalledTimes(3);
  });

  it('gives up after the retry budget and reports the status, not the HTML body', async () => {
    stub(() => shed(503));
    const client = new NodeApiClient(URL, 'token');
    await expect(client.get('/projects/x')).rejects.toThrow('HTTP 503');
  });

  it('does not retry a status the server meant', async () => {
    const spy = stub(() => new Response('{"message":"nope"}', {status: 403, headers: {'content-type': 'application/json'}}));
    const client = new NodeApiClient(URL, 'token');
    await expect(client.get('/projects/x')).rejects.toThrow('nope');
    expect(spy).toHaveBeenCalledTimes(1);
  });

  it('does not retry a write — a POST the server may have applied is not repeated', async () => {
    const spy = stub(() => shed(503));
    const client = new NodeApiClient(URL, 'token');
    await expect(client.post('/projects', {id: 'x'})).rejects.toThrow('HTTP 503');
    expect(spy).toHaveBeenCalledTimes(1);
  });

  it('names the route when a request never answers', async () => {
    stub(() => { throw Object.assign(new Error('timed out'), {name: 'TimeoutError'}); });
    const client = new NodeApiClient(URL, 'token');
    await expect(client.get('/projects/x')).rejects.toThrow('GET http://stand/api/projects/x: no answer in');
  });
});
