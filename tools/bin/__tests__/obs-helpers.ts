/// Test-only: a mocked server per host and captured console output for the `grok s` observability commands.
import {vi} from 'vitest';
import {NodeDapi} from '../utils/node-dapi';
import {Connect} from '../utils/server-client';

export interface Call {host: string; method: string; path: string; body?: any}

export type Responder = (method: string, path: string, body: any, host: string) => any;

/** A [Connect] whose clients answer through [responder] and record every call. */
export function mockConnect(responder: Responder): {connect: Connect; calls: Call[]} {
  const calls: Call[] = [];
  const connect: Connect = async (host?: string) => {
    const h = host ?? '';
    const client: any = {
      async request(method: string, path: string, body?: any) {
        calls.push({host: h, method, path, body});
        return responder(method, path, body, h);
      },
      get(path: string) { return this.request('GET', path); },
      post(path: string, body?: any) { return this.request('POST', path, body); },
      del(path: string) { return this.request('DELETE', path); },
    };
    return new NodeDapi(client);
  };
  return {connect, calls};
}

export function apiError(status: number, message: string): Error {
  return Object.assign(new Error(message), {apiError: {error: message, errorCode: status}});
}

/** Lines written to stdout (console.log and process.stdout.write) and stderr during [fn]. */
export async function captureOutput(fn: () => Promise<any>): Promise<{out: string[]; err: string[]; result: any; exitCode: any}> {
  const out: string[] = [];
  const err: string[] = [];
  const log = vi.spyOn(console, 'log').mockImplementation((...a: any[]) => { out.push(a.join(' ')); });
  const write = vi.spyOn(process.stdout, 'write').mockImplementation((chunk: any) => { out.push(String(chunk)); return true; });
  const errWrite = vi.spyOn(process.stderr, 'write').mockImplementation((chunk: any) => { err.push(String(chunk).replace(/\n$/, '')); return true; });
  const before = process.exitCode;
  process.exitCode = undefined;
  try {
    const result = await fn();
    return {out, err, result, exitCode: process.exitCode};
  }
  finally {
    log.mockRestore();
    write.mockRestore();
    errWrite.mockRestore();
    process.exitCode = before;
  }
}

/** A local time today (or [daysAgo] days back) as ISO, so time formatting is independent of the zone. */
export function localIso(h: number, m: number, daysAgo: number = 0): string {
  const d = new Date();
  d.setDate(d.getDate() - daysAgo);
  d.setHours(h, m, 0, 0);
  return d.toISOString();
}
