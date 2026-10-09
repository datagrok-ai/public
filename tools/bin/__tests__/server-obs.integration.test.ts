/**
 * Live smoke of `grok s observe`: one run of each command through the real CLI against a server.
 *
 * Environment variables:
 *   GROK_IT_HOST   Server alias (required; the suite is skipped without it). The logger check adds
 *                  an entry for the admin's own group that expires in 5 minutes and removes it again.
 *
 * Run: GROK_IT_HOST=wt-GROK-20884-obs npx vitest run --project integration server-obs
 */

import {describe, expect, it} from 'vitest';
import {spawnSync} from 'child_process';
import * as path from 'path';

const HOST = process.env['GROK_IT_HOST'] ?? '';
const GROK = path.resolve(__dirname, '..', 'grok.js');

function observe(...args: string[]): any {
  const r = spawnSync(process.execPath, [GROK, 's', 'o', ...args, '--host', HOST, '--output', 'json'],
    {encoding: 'utf8', timeout: 60_000});
  if (r.status !== 0)
    throw new Error(`grok s o ${args.join(' ')} exited ${r.status}\n${r.stdout}\n${r.stderr}`);
  return r.stdout.trim().startsWith('[') || r.stdout.trim().startsWith('{') ? JSON.parse(r.stdout) : r.stdout;
}

describe.skipIf(!HOST)('grok s observe', () => {
  it('problems list, get and history', () => {
    const problems: any[] = observe('problems', 'list', '--status', 'all', '--limit', '5');
    expect(Array.isArray(problems)).toBe(true);
    if (!problems.length)
      return;
    expect(observe('problems', 'get', problems[0].id).id).toBe(problems[0].id);
    expect(Array.isArray(observe('problems', 'history', `${problems[0].kind}:${problems[0].key}`, '--limit', '3'))).toBe(true);
  });

  it('rules get and test', () => {
    expect(Array.isArray(observe('rules', 'get'))).toBe(true);
    expect(typeof observe('rules', 'test')).toBe('string');
  });

  it('logger get, set an expiring entry, history', () => {
    const me = spawnSync(process.execPath, [GROK, 's', 'raw', 'GET', '/users/current', '--host', HOST, '--output', 'json'],
      {encoding: 'utf8'});
    const login = JSON.parse(me.stdout).login;
    const before = observe('logger', 'get');
    expect(before.userGroupSettings).toBeDefined();
    try {
      observe('logger', 'set', '--user', login, '--set', 'saveHttpRequests=true', '--for', '5m', '--reason', 'grok s smoke');
      const after = observe('logger', 'get', 'userGroupSettings');
      expect(Object.keys(after).find((id) => after[id].reason === 'grok s smoke')).toBeDefined();
      const history: any[] = observe('logger', 'history', '--limit', '1');
      expect(history[0].diff.some((d: any) => String(d.path).includes('saveHttpRequests'))).toBe(true);
    }
    finally {
      observe('logger', 'set', '--set', `userGroupSettings=${JSON.stringify(before.userGroupSettings)}`);
    }
  });

  it('errors top', () => {
    const errors = observe('errors', 'top', '--since', '7d', '--by', 'signature', '--limit', '3');
    expect(errors.now).toBeDefined();
    expect(errors.top.length).toBeLessThanOrEqual(3 * Math.max(1, Object.keys(errors.bySource).length));
  });

  it('timeline', () => {
    const events: any[] = JSON.parse(spawnSync(process.execPath, [GROK, 's', 'raw', 'GET', '/log?limit=50&page=1',
      '--host', HOST, '--output', 'json'], {encoding: 'utf8'}).stdout);
    const request = events.find((e) => e.requestId)?.requestId;
    if (!request)
      return;
    const rows: any[] = observe('timeline', '--request', request, '--limit', '5');
    expect(rows.some((r) => r.requestId === request)).toBe(true);
  });
});
