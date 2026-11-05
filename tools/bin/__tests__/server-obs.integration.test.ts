/**
 * Live smoke of `grok s observe`: one run of each command through the real CLI against a server.
 *
 * Environment variables:
 *   GROK_IT_HOST   Server alias (required; the suite is skipped without it). The logger and capture checks add
 *                  an entry for the admin's own group that expires in 5 minutes and remove it again; they skip
 *                  when that group already has an entry.
 *
 * Run: GROK_IT_HOST=wt-GROK-20884-obs npx vitest run --project integration server-obs
 */

import {describe, expect, it} from 'vitest';
import {spawnSync} from 'child_process';
import * as path from 'path';

const HOST = process.env['GROK_IT_HOST'] ?? '';
const GROK = path.resolve(__dirname, '..', 'grok.js');

function grok(...args: string[]): any {
  const r = spawnSync(process.execPath, [GROK, 's', ...args, '--host', HOST, '--output', 'json'],
    {encoding: 'utf8', timeout: 60_000});
  if (r.status !== 0)
    throw new Error(`grok s ${args.join(' ')} exited ${r.status}\n${r.stdout}\n${r.stderr}`);
  const out = r.stdout.trim();
  return out.startsWith('[') || out.startsWith('{') ? JSON.parse(out) : out;
}

const observe = (...args: string[]): any => grok('o', ...args);

const sleep = (ms: number) => Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, ms);

/** The admin's login and personal group, or the reason the logger checks cannot run. */
function me(): {login: string, group: string, skip?: string} {
  const user = grok('raw', 'GET', '/users/current');
  const group = user.group?.id;
  const settings = observe('logger', 'get', 'userGroupSettings') ?? {};
  return {login: user.login, group, skip: !group ? '/users/current has no group id'
    : settings[group] ? `${user.login} already has a logger entry; not touching it` : undefined};
}

describe.skipIf(!HOST)('grok s observe', () => {
  it('problems list, get and history', (ctx) => {
    const problems: any[] = observe('problems', 'list', '--status', 'all', '--limit', '5');
    expect(Array.isArray(problems)).toBe(true);
    if (!problems.length)
      ctx.skip('The server has no problems to read');
    expect(observe('problems', 'get', problems[0].id).id).toBe(problems[0].id);
    expect(Array.isArray(observe('problems', 'history', `${problems[0].kind}:${problems[0].key}`, '--limit', '3'))).toBe(true);
  });

  it('rules get and test', () => {
    expect(Array.isArray(observe('rules', 'get'))).toBe(true);
    expect(observe('rules', 'test')).toMatch(/^(No rules|Rule )/);
  });

  it('logger get, set an expiring entry, history', (ctx) => {
    const {login, group, skip} = me();
    if (skip)
      ctx.skip(skip);
    try {
      observe('logger', 'set', '--user', login, '--set', 'saveHttpRequests=true', '--for', '5m', '--reason', 'grok s smoke');
      expect(observe('logger', 'get', `userGroupSettings.${group}`).reason).toBe('grok s smoke');
      let ours: any;
      for (let i = 0; i < 10 && !ours; i++) {
        if (i)
          sleep(1000);
        const history: any[] = observe('logger', 'history', '--limit', '20');
        ours = history.find((h) => Array.isArray(h.diff) &&
          h.diff.some((d: any) => d.path === `userGroupSettings.${group}.reason` && d.after === 'grok s smoke'));
      }
      expect(ours?.diff.some((d: any) => d.path === `userGroupSettings.${group}.saveHttpRequests`)).toBe(true);
    }
    finally {
      observe('logger', 'set', '--user', login, '--remove');
    }
    expect(observe('logger', 'get', `userGroupSettings.${group}`)).toBe('(empty)');
  });

  it('logs capture and capture stop', (ctx) => {
    const {login, group, skip} = me();
    if (skip)
      ctx.skip(skip);
    try {
      observe('logs', 'capture', login, '--for', '5m', '--reason', 'grok s smoke capture');
      const captures: any[] = observe('logs', 'capture');
      const ours = captures.find((c) => c.id === group);
      expect(ours?.reason).toBe('grok s smoke capture');
      expect(Date.parse(ours.expiresAt) - Date.now()).toBeGreaterThan(4 * 60000);
      observe('logs', 'capture', 'stop', login);
      expect(observe('logger', 'get', `userGroupSettings.${group}`)).toBe('(empty)');
    }
    finally {
      if (observe('logger', 'get', `userGroupSettings.${group}`) !== '(empty)')
        observe('logger', 'set', '--user', login, '--remove');
    }
  });

  it('errors top', () => {
    const errors = observe('errors', 'top', '--since', '7d', '--by', 'signature', '--limit', '3');
    expect(errors.now).toBeDefined();
    expect(errors.top.length).toBeLessThanOrEqual(3 * Math.max(1, Object.keys(errors.bySource).length));
  });

  it('timeline', (ctx) => {
    const events: any[] = grok('raw', 'GET', '/log?limit=50&page=1');
    const event = events.find((e) => e.session?.id && e.eventTime);
    if (!event)
      ctx.skip('No recent event with a session');
    const rows: any[] = observe('timeline', '--session', event.session.id, `--from=${event.eventTime}`, '--limit', '5');
    expect(Array.isArray(rows)).toBe(true);
  });
});
