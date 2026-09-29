import {describe, it, expect, beforeEach, afterEach, vi} from 'vitest';
import {handleCapture, handleTimeline, parseCapture, captureBody, ruleRow, ruleSummary, timelineQuery, timelineRow} from '../commands/server-capture';
import {mockConnect, captureOutput, utcIso} from './obs-helpers';

const NOW = new Date(Date.UTC(2026, 8, 28, 10, 20));

const RULE = {id: 'r-17', number: 17, author: 'b.ivanov', subject: {type: 'user', value: 'alice.mendel'},
  scope: {type: 'view', value: 'Hit Triage'}, capture: {clicks: true, inputs: true, requests: true, calls: true, errors: true,
    serverLevel: 'debug', debugFlags: ['query', 'storage']}, reason: 'GROK-21044: campaign loses filters', status: 'active',
  createdAt: new Date(Date.UTC(2026, 8, 28, 10, 20)).toISOString(), expiresAt: new Date(Date.UTC(2026, 8, 30, 10, 20)).toISOString(),
  events: 0, maxEvents: 2000, windowMinutes: 10, maxSessions: 20};

describe('parseCapture', () => {
  it('reads plain items and a trailing server level with flags', () => {
    expect(parseCapture('clicks,inputs,requests,calls,errors,server:debug=queries,files')).toEqual({
      clicks: true, inputs: true, requests: true, calls: true, errors: true, serverLevel: 'debug', debugFlags: ['query', 'storage']});
    expect(parseCapture('server:info')).toMatchObject({clicks: false, serverLevel: 'info', debugFlags: []});
  });

  it('refuses unknown items and the credentials flag; everything after server: is flags', () => {
    expect(() => parseCapture('clicks,screens')).toThrow(/Unknown capture item 'screens'/);
    expect(parseCapture('server:debug=queries,clicks').debugFlags).toEqual(['query', 'clicks']);
    expect(() => parseCapture('server:debug=credentials')).toThrow(/credentials/);
    expect(() => parseCapture(undefined)).toThrow(/needs items/);
  });
});

describe('captureBody', () => {
  const base = {user: 'alice.mendel', view: 'Hit Triage', capture: 'clicks,errors', for: '2d', reason: 'GROK-21044'};

  it('builds the rule for example 6', () => {
    expect(captureBody({...base, capture: 'clicks,inputs,requests,calls,errors,server:debug=queries,files', limit: 2000}, NOW)).toEqual({
      name: undefined, subject: {type: 'user', value: 'alice.mendel'}, scope: {type: 'view', value: 'Hit Triage'},
      capture: {clicks: true, inputs: true, requests: true, calls: true, errors: true, serverLevel: 'debug', debugFlags: ['query', 'storage']},
      anonymous: false, reason: 'GROK-21044', maxEvents: 2000, forMinutes: 2880});
  });

  it('builds an anonymous group rule with an absolute expiry', () => {
    const body = captureBody({group: 'Chemists', view: 'Hit Triage', capture: 'clicks,requests,errors', until: '2026-10-05T10:00',
      anonymous: true, reason: 'submit drop-off', window: '5m', 'max-sessions': 50}, NOW);
    expect(body).toMatchObject({subject: {type: 'group', value: 'Chemists'}, anonymous: true, windowMinutes: 5, maxSessions: 50,
      expiresAt: new Date(Date.UTC(2026, 9, 5, 10, 0)).toISOString()});
    expect(body.forMinutes).toBeUndefined();
  });

  it('refuses what the platform refuses', () => {
    expect(() => captureBody({...base, reason: undefined}, NOW)).toThrow(/--reason/);
    expect(() => captureBody({...base, for: undefined}, NOW)).toThrow(/--for <duration> or --until/);
    expect(() => captureBody({...base, anonymous: true}, NOW)).toThrow(/--anonymous applies to --group and --everyone/);
    expect(() => captureBody({...base, group: 'Chemists'}, NOW)).toThrow(/exactly one subject/);
    expect(() => captureBody({...base, user: undefined}, NOW)).toThrow(/exactly one subject/);
    expect(() => captureBody({...base, element: 'x'}, NOW)).toThrow(/at most one scope/);
    expect(() => captureBody({...base, user: undefined, everyone: true, view: undefined}, NOW)).toThrow(/needs a scope/);
    expect(() => captureBody({...base, for: undefined, until: '2026-09-01T00:00'}, NOW)).toThrow(/in the past/);
  });
});

describe('rows', () => {
  it('prints the one-line summary of a new rule', () => {
    expect(ruleSummary(RULE)).toBe('rule cap-17  active until 2026-09-30 10:20Z · 1 user · 1 view · 0/2000 events');
    expect(ruleSummary({...RULE, subject: {type: 'everyone'}, scope: null})).toMatch(/· everyone · all activity ·/);
  });

  it('lists rules with the example 15 columns', () => {
    const ended = {...RULE, status: 'expired', endedAt: RULE.expiresAt, events: 340};
    expect(ruleRow(ended)).toEqual({RULE: 'cap-17', AUTHOR: 'b.ivanov', SUBJECT: 'user alice.mendel', SCOPE: 'view Hit Triage',
      REASON: 'GROK-21044: campaign loses filters', ACTIVE: '2 d (expired)', EVENTS: 340});
    const server = {...RULE, number: 12, subject: {type: 'package', value: 'Snowflake'}, scope: null,
      expiresAt: new Date(Date.UTC(2026, 8, 28, 10, 50)).toISOString()};
    expect(ruleRow(server)).toMatchObject({SUBJECT: 'package Snowflake', SCOPE: 'server debug', ACTIVE: '30 min'});
    expect(ruleRow({...server, capture: {clicks: true}})).toMatchObject({SCOPE: 'all activity'});
  });

  it('renders timeline rows with milliseconds and short request ids', () => {
    const at = new Date();
    at.setHours(10, 14, 3, 112);
    expect(timelineRow({time: at.toISOString(), source: 'datlas', kind: 'request', summary: 'POST /api/projects/{id}/save',
      status: 403, ms: 28, requestId: 'mfz3k2a1b9x8y7kq.1'}))
      .toEqual({TIME: '09:14:03.112Z', SOURCE: 'datlas', KIND: 'request', SUMMARY: 'POST /api/projects/{id}/save', STATUS: 403, MS: 28, REQ: '…x8y7kq.1'});
  });

  it('takes exactly one timeline key', () => {
    expect(timelineQuery({report: 4820})).toEqual({report: '4820', limit: undefined});
    expect(timelineQuery({rule: 'cap-17', from: '-1h'}, NOW)).toEqual({rule: 'cap-17', limit: undefined, from: new Date(Date.UTC(2026, 8, 28, 9, 20)).toISOString()});
    expect(() => timelineQuery({})).toThrow(/exactly one of --action/);
    expect(() => timelineQuery({session: '', report: 1})).not.toThrow();
    expect(() => timelineQuery({action: 'a', session: 's'})).toThrow(/exactly one/);
  });
});

describe('handlers', () => {
  beforeEach(() => { vi.useFakeTimers({toFake: ['Date']}); vi.setSystemTime(NOW); });
  afterEach(() => vi.useRealTimers());

  it('adds a rule and prints its summary', async () => {
    const {connect, calls} = mockConnect(() => RULE);
    const {out} = await captureOutput(() => handleCapture(connect, 'add', [],
      {user: 'alice.mendel', view: 'Hit Triage', capture: 'clicks,errors', for: '2d', limit: 2000, reason: 'GROK-21044'}, 'table'));
    expect(calls[0]).toMatchObject({method: 'POST', path: '/logging/capture'});
    expect(out).toEqual(['rule cap-17  active until 2026-09-30 10:20Z · 1 user · 1 view · 0/2000 events']);
  });

  it('lists all rules since a time', async () => {
    const {connect, calls} = mockConnect(() => [RULE]);
    const {out} = await captureOutput(() => handleCapture(connect, 'list', [], {all: true, since: '90d'}, 'table'));
    expect(calls[0].path).toBe('/logging/capture?all=true&since=90d');
    expect(out[0].trimEnd().split(/\s{2,}/)).toEqual(['RULE', 'AUTHOR', 'SUBJECT', 'SCOPE', 'REASON', 'ACTIVE', 'EVENTS']);
  });

  it('shows a rule, its activations, and its timeline as csv', async () => {
    const {connect, calls} = mockConnect((_m, p) => p.startsWith('/log/timeline')
      ? [{time: utcIso(15, 2), source: 'client', kind: 'click', summary: 'Hit Triage / Filters / Reset', requestId: 'mfz3k2a1b9x8y7kq'}]
      : {...RULE, activations: [{sessionId: 'abcdef0123', user: 'alice.mendel', activatedAt: RULE.createdAt, until: RULE.expiresAt, triggerDetail: 'view Hit Triage'}]});
    const show = await captureOutput(() => handleCapture(connect, 'show', ['cap-17'], {}, 'table'));
    expect(calls[0].path).toBe('/logging/capture/cap-17');
    expect(show.out[0]).toMatch(/^rule\s+cap-17$/);
    expect(show.out.some((l) => /^capture\s+clicks,inputs,requests,calls,errors,server:debug=query,storage$/.test(l))).toBe(true);
    expect(show.out.some((l) => l.startsWith('abcdef01'))).toBe(true);
    const csv = await captureOutput(() => handleCapture(connect, 'show', ['cap-17'], {timeline: true}, 'csv'));
    expect(calls[1].path).toBe('/log/timeline?rule=cap-17');
    expect(csv.out[0]).toBe('TIME,SOURCE,KIND,SUMMARY,STATUS,MS,REQ');
    expect(csv.out[1]).toMatch(/,client,click,Hit Triage \/ Filters \/ Reset,,,…x8y7kq$/);
  });

  it('stops a rule with a reason', async () => {
    const {connect, calls} = mockConnect(() => ({...RULE, status: 'stopped'}));
    const {out} = await captureOutput(() => handleCapture(connect, 'stop', ['cap-17'], {reason: 'reproduced'}, 'table'));
    expect(calls[0]).toMatchObject({method: 'DELETE', path: '/logging/capture/cap-17', body: {reason: 'reproduced'}});
    expect(out).toEqual(['stopped cap-17 — reproduced']);
  });

  it('prints a report timeline and refuses a missing key', async () => {
    const {connect, calls} = mockConnect(() => []);
    await captureOutput(() => handleTimeline(connect, undefined, [], {report: 4820}, 'table'));
    expect(calls[0].path).toBe('/log/timeline?report=4820');
    const bad = await captureOutput(() => handleTimeline(connect, undefined, [], {}, 'table'));
    expect(bad.result).toBe(false);
    expect(bad.err.join('\n')).toMatch(/Usage: grok s timeline/);
  });
});
