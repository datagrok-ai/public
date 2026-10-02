import {describe, it, expect} from 'vitest';
import {handleProblems, problemRow} from '../commands/server-problems';
import {mockConnect, captureOutput, utcIso} from './obs-helpers';

const PROBLEM = {id: 'b52e0d00-0000-4000-8000-000000000001', kind: 'error-incident', key: 'a41f9c', name: 'DatagrokErrorIncident',
  severity: 'warning', status: 'active', state: 'ongoing', episodes: 3, occurrences: 41, lastSeen: utcIso(10, 5),
  summary: 'TypeError: Cannot read properties of undefined (reading \'molfile\') in Chem, 9 users in 15 min'};

describe('problemRow', () => {
  it('prints the fixed problem columns, summary cut at 60', () => {
    const row = problemRow(PROBLEM);
    expect(Object.keys(row)).toEqual(['KIND', 'KEY', 'SEV', 'STATUS', 'STATE', 'EPISODES', 'LAST SEEN', 'SUMMARY']);
    expect(row).toMatchObject({STATUS: 'active', STATE: 'ongoing', EPISODES: 3, 'LAST SEEN': '10:05Z'});
    expect(row.SUMMARY.length).toBe(60);
  });
});

describe('handleProblems', () => {
  it('lists with status and state filters on several hosts', async () => {
    const {connect, calls} = mockConnect(() => [PROBLEM]);
    const {out} = await captureOutput(() => handleProblems(connect, 'list', [], {status: 'muted,fixed', state: 'ongoing',
      host: ['prod', 'val']}, 'table'));
    expect(calls.map((c) => [c.host, c.path])).toEqual([['prod', '/problems?status=muted%2Cfixed&state=ongoing'],
      ['val', '/problems?status=muted%2Cfixed&state=ongoing']]);
    expect(out[0].split(/\s+/).slice(0, 3)).toEqual(['HOST', 'KIND', 'KEY']);
  });

  it('looks kind:key up whatever the status, then sets the status by id', async () => {
    const {connect, calls} = mockConnect((m) => m === 'GET' ? [PROBLEM] : {...PROBLEM, status: 'not-a-problem'});
    const {out} = await captureOutput(() => handleProblems(connect, 'dismiss', ['error-incident:a41f9c'],
      {reason: 'expected when a token expires'}, 'table'));
    expect(calls.map((c) => [c.method, c.path, c.body])).toEqual([
      ['GET', '/problems?kind=error-incident&key=a41f9c&status=all', undefined],
      ['POST', `/problems/${PROBLEM.id}/status`, {status: 'not-a-problem', reason: 'expected when a token expires'}],
    ]);
    expect(out).toEqual(['dismissed as not a problem: error-incident:a41f9c — expected when a token expires']);
  });

  it('mutes with an end, or until lifted', async () => {
    const {connect, calls} = mockConnect(() => ({...PROBLEM, status: 'muted'}));
    const {out} = await captureOutput(() => handleProblems(connect, 'mute', [PROBLEM.id], {reason: 'fixed in 1.14.3',
      'until-version': '1.14.3'}, 'table'));
    expect(calls[0].body).toEqual({status: 'muted', reason: 'fixed in 1.14.3', untilVersion: '1.14.3'});
    expect(out).toEqual(['muted error-incident:a41f9c until 1.14.3 — fixed in 1.14.3']);
    await captureOutput(() => handleProblems(connect, 'mute', [PROBLEM.id], {reason: 'noise'}, 'table'));
    expect(calls[1].body).toEqual({status: 'muted', reason: 'noise'});
  });

  it('fixes and activates with an optional reason; dismiss needs one', async () => {
    const {connect, calls} = mockConnect((_m, _p, body) => ({...PROBLEM, status: body.status}));
    await captureOutput(() => handleProblems(connect, 'fix', [PROBLEM.id], {}, 'table'));
    await captureOutput(() => handleProblems(connect, 'activate', [PROBLEM.id], {reason: 'back'}, 'table'));
    expect(calls.map((c) => c.body)).toEqual([{status: 'fixed', reason: undefined}, {status: 'active', reason: 'back'}]);
    const {result} = await captureOutput(() => handleProblems(connect, 'dismiss', [PROBLEM.id], {}, 'table'));
    expect(result).toBe(false);
    expect(calls.length).toBe(2);
  });

  it('lists the alerts of a problem, every status by default', async () => {
    const {connect, calls} = mockConnect(() => [{id: 'A1', kind: 'error-incident', key: 'a41f9c', status: 'resolved'}]);
    await captureOutput(() => handleProblems(connect, 'alerts', [PROBLEM.id], {}, 'table'));
    expect(calls[0].path).toBe(`/alerts?problem=${PROBLEM.id}&status=all`);
  });

  it('refuses a kind:key that names no problem', async () => {
    const {connect} = mockConnect(() => []);
    await expect(handleProblems(connect, 'get', ['health:Nope'], {}, 'table')).rejects.toThrow('No problem health:Nope');
  });

  it('answers an unknown verb with the usage', async () => {
    const {connect} = mockConnect(() => ({}));
    const {err, result} = await captureOutput(() => handleProblems(connect, 'snooze', [], {}, 'table'));
    expect(result).toBe(false);
    expect(err.join('\n')).toMatch(/Usage: grok s o problems/);
  });
});
