import {AddressInfo} from 'node:net';
import {WebSocketServer, WebSocket} from 'ws';

import {Const, PipeClient} from '../pipe-client';

interface PeerFrame {
  data: string | Uint8Array;
  isBinary: boolean;
}

/** Scripted grok_pipe stand-in: records inbound frames and lets tests await them. */
class ScriptedPeer {
  wss!: WebSocketServer;
  socket: WebSocket | null = null;
  headers: any = null;
  frames: PeerFrame[] = [];
  private waiters: ((frame: PeerFrame) => void)[] = [];
  private processed = 0;

  static async start(): Promise<ScriptedPeer> {
    const peer = new ScriptedPeer();
    peer.wss = new WebSocketServer({port: 0});
    await new Promise<void>((resolve) => peer.wss.on('listening', () => resolve()));
    peer.wss.on('connection', (ws, req) => {
      peer.socket = ws;
      peer.headers = req.headers;
      ws.on('message', (data: any, isBinary: boolean) => {
        const frame: PeerFrame = {
          data: isBinary ? new Uint8Array(data as Buffer) : (data as Buffer).toString('utf8'),
          isBinary: isBinary,
        };
        peer.frames.push(frame);
        const waiter = peer.waiters.shift();
        if (waiter != null)
          waiter(frame);
      });
    });
    return peer;
  }

  get port(): number {
    return (this.wss.address() as AddressInfo).port;
  }

  nextFrame(): Promise<PeerFrame> {
    if (this.frames.length > this.processed)
      return Promise.resolve(this.frames[this.processed++]);
    return new Promise((resolve) => this.waiters.push((frame) => {
      this.processed++;
      resolve(frame);
    }));
  }

  send(message: string | Uint8Array): void {
    this.socket!.send(message, {binary: typeof message !== 'string'});
  }

  async stop(): Promise<void> {
    for (const ws of this.wss.clients)
      ws.terminate();
    await new Promise<void>((resolve) => this.wss.close(() => resolve()));
  }
}

function makeClient(peer: ScriptedPeer, batchSize: number = 4): PipeClient {
  return new PipeClient(`ws://127.0.0.1:${peer.port}/call-1`, {
    'x-member-name': 'celery-datagrok-celery',
    'authorization': 'test-key',
  }, {messageTimeoutMs: 2000, paramTimeoutMs: 10000, batchSize: batchSize, connectAttempts: 1});
}

describe('PipeClient', () => {
  let peer: ScriptedPeer;
  let client: PipeClient;

  beforeEach(async () => {
    peer = await ScriptedPeer.start();
    client = makeClient(peer);
    await client.connect();
  });

  afterEach(async () => {
    client.close();
    client.close(); // safe to call twice
    await peer.stop();
  });

  test('sends the room-join headers', () => {
    expect(peer.headers['x-member-name']).toBe('celery-datagrok-celery');
    expect(peer.headers['authorization']).toBe('test-key');
  });

  test('receiveParam: PARAM -> SENDING/chunks/PARAM_SENT with PART OK cadence, unrelated frames skipped', async () => {
    const chunk1 = new Uint8Array([1, 2, 3]);
    const chunk2 = new Uint8Array([4, 5]);
    const resultPromise = client.receiveParam('df');

    expect((await peer.nextFrame()).data).toBe('PARAM df');
    peer.send('LOG {"message":"unrelated"}'); // must be tolerated and skipped
    peer.send(`${Const.SENDING} 5 {".id":"x",".type":"csv"}`);
    peer.send(chunk1);
    expect((await peer.nextFrame()).data).toBe(Const.PART_OK);
    peer.send(chunk2);
    expect((await peer.nextFrame()).data).toBe(Const.PART_OK);
    peer.send('PARAM_SENT df');

    const result = await resultPromise;
    expect(Array.from(result!)).toEqual([1, 2, 3, 4, 5]);
  });

  test('receiveParam returns null when the peer sends no data', async () => {
    const resultPromise = client.receiveParam('df');
    expect((await peer.nextFrame()).data).toBe('PARAM df');
    peer.send('PARAM_SENT df');
    await expect(resultPromise).resolves.toBeNull();
  });

  test('receiveParam aborts on the Dart PART ERROR constant', async () => {
    const resultPromise = client.receiveParam('df');
    await peer.nextFrame();
    peer.send(Const.PART_ERROR);
    await expect(resultPromise).rejects.toThrow(/error while sending param df/);
  });

  test('receiveParam aborts on the bare python ERROR constant', async () => {
    const resultPromise = client.receiveParam('df');
    await peer.nextFrame();
    peer.send(Const.ERROR);
    await expect(resultPromise).rejects.toThrow(/error while sending param df/);
  });

  test('sendParam: chunks by batchSize and waits for PART OK after each chunk', async () => {
    const data = new Uint8Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]); // batchSize 4 -> 4+4+2
    const sendPromise = client.sendParam(data, {'.id': 'out', '.type': 'csv'});

    const header = await peer.nextFrame();
    expect(header.data).toBe(`${Const.SENDING} 10 {".id":"out",".type":"csv"}`);

    const received: number[] = [];
    for (let i = 0; i < 3; i++) {
      const frame = await peer.nextFrame();
      expect(frame.isBinary).toBe(true);
      received.push(...(frame.data as Uint8Array));
      if (i === 0)
        peer.send('PING REPLY'); // unrelated frame while waiting for the ack — skipped
      peer.send(Const.PART_OK);
    }
    await sendPromise;
    expect(received).toEqual(Array.from(data));
    expect(peer.frames.filter((f) => f.isBinary).map((f) => (f.data as Uint8Array).length)).toEqual([4, 4, 2]);
  });

  test('sendParam aborts on PART ERROR and on bare ERROR', async () => {
    const sendPromise = client.sendParam(new Uint8Array([1, 2, 3]), {'.id': 'out', '.type': 'blob'});
    await peer.nextFrame(); // SENDING header
    await peer.nextFrame(); // binary chunk
    peer.send(Const.PART_ERROR);
    await expect(sendPromise).rejects.toThrow(/error while receiving/);

    const sendPromise2 = client.sendParam(new Uint8Array([1]), {'.id': 'out', '.type': 'blob'});
    await peer.nextFrame();
    await peer.nextFrame();
    peer.send(Const.ERROR);
    await expect(sendPromise2).rejects.toThrow(/error while receiving/);
  });

  test('sendText is used for CALL/PROGRESS frames', async () => {
    client.sendText('CALL {"id":"call-1","status":"Completed"}');
    expect((await peer.nextFrame()).data).toBe('CALL {"id":"call-1","status":"Completed"}');
  });
});
