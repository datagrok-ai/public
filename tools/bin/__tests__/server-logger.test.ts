import {describe, it, expect} from 'vitest';
import {handleLogger, parseScope, isOverride, propChanges, setArgs, revertBody, displayPath, diffMaps} from '../commands/server-logger';
import {applyListSpec, normalizeFlag, normalizeLevel} from '../utils/obs-format';
import {printError} from '../utils/server-output';
import {mockConnect, captureOutput, utcIso, apiError} from './obs-helpers';

const ALL = 'a4b45840-0000-0000-0000-00000000a11a';
const CHEM = 'c4e00000-0000-0000-0000-000000000c4e';

const POLICY = {
  version: 7,
  settings: {
    [`userGroupSettings.${ALL}.printLevels`]: ['error', 'warning', 'info'],
    [`userGroupSettings.${ALL}.postLevels`]: ['error', 'warning', 'audit', 'usage'],
    [`userGroupSettings.${ALL}.saveLevels`]: ['error', 'warning', 'audit', 'usage'],
    [`userGroupSettings.${ALL}.debugFlags`]: [],
    [`userGroupSettings.${CHEM}.saveLevels`]: ['error', 'warning', 'audit', 'usage', 'info'],
    'exportSettings.d1.endpoint': 'https://otel.acme.internal:4318',
  },
  defaults: {
    [`userGroupSettings.${ALL}.printLevels`]: ['error', 'warning', 'info'],
    [`userGroupSettings.${ALL}.postLevels`]: ['error', 'warning', 'audit', 'usage'],
    [`userGroupSettings.${ALL}.saveLevels`]: ['error', 'warning', 'audit', 'usage'],
    [`userGroupSettings.${ALL}.debugFlags`]: [],
  },
  locked: ['exportSettings', 'levels.audit'],
  overrides: [{id: 'o1', scope: 'package', scopeId: 'Snowflake', scopeName: 'Snowflake', changes: {debugFlags: ['query']},
    expiresAt: utcIso(11, 42), reason: 'prod timeouts', createdByLogin: 'k.lee'}],
  groups: {[ALL]: 'All users', [CHEM]: 'Chemists'},
  destinations: {d1: 'OpenTelemetry'},
};

describe('lists', () => {
  it('replaces with a plain list and edits with signed items', () => {
    expect(applyListSpec(['error'], 'error,warning', normalizeLevel, '--save-levels')).toEqual(['error', 'warning']);
    expect(applyListSpec(['error', 'audit'], '-audit,+info', normalizeLevel, '--save-levels')).toEqual(['error', 'info']);
    expect(applyListSpec([], '+queries,+files', normalizeFlag, '--debug-flags')).toEqual(['query', 'storage']);
    expect(applyListSpec(['query'], '+query', normalizeFlag, '--debug-flags')).toEqual(['query']);
  });

  it('refuses mixed forms and leaves unknown names to the server', () => {
    expect(() => applyListSpec([], 'error,+info', normalizeLevel, '--save-levels')).toThrow(/not both/);
    expect(applyListSpec([], '+queries,+connections', normalizeFlag, '--debug-flags')).toEqual(['query', 'connections']);
  });
});

describe('scope and target', () => {
  it('parses scopes', () => {
    expect(parseScope(undefined)).toEqual({type: 'all', value: undefined, label: 'all'});
    expect(parseScope('package:Snowflake')).toEqual({type: 'package', value: 'Snowflake', label: 'package:Snowflake'});
    expect(parseScope('group:Lab: West').value).toBe('Lab: West');
    expect(() => parseScope('team:x')).toThrow(/--scope is all/);
    expect(() => parseScope('user:')).toThrow(/--scope is all/);
    expect(() => parseScope('all:x')).toThrow(/--scope is all/);
  });

  it('decides base change or override', () => {
    expect(isOverride({}, parseScope('all'))).toBe(false);
    expect(isOverride({}, parseScope('group:Chemists'))).toBe(false);
    expect(isOverride({for: '30m'}, parseScope('all'))).toBe(true);
    expect(isOverride({until: '2026-10-01T10:00'}, parseScope('group:Chemists'))).toBe(true);
    expect(isOverride({}, parseScope('user:alice'))).toBe(true);
    expect(isOverride({}, parseScope('session:s1'))).toBe(true);
    expect(isOverride({}, parseScope('package:Chem'))).toBe(true);
  });

  it('reads property flags and --set pairs', () => {
    expect(propChanges({'print-format': 'json', 'print-details': 'false', 'post-levels': 'error'}, () => []))
      .toEqual({postLevels: ['error'], printFormat: 'json', printDetails: false});
    expect(() => propChanges({'print-format': 'xml'}, () => [])).toThrow(/text or json/);
    expect(setArgs(['exportFlushSeconds=5', 'exportSettings.d1.endpoint=https://x'])).toEqual({exportFlushSeconds: 5, 'exportSettings.d1.endpoint': 'https://x'});
    expect(() => setArgs('novalue')).toThrow(/<path>=<json>/);
  });

  it('builds revert bodies', () => {
    expect(revertBody([], {})).toEqual({reason: undefined});
    expect(revertBody(['5'], {reason: 'bad'})).toEqual({reason: 'bad', version: 5});
    expect(revertBody([], {override: 'o1'})).toEqual({reason: undefined, override: 'o1'});
    expect(revertBody([], {overrides: true})).toEqual({reason: undefined, allOverrides: true});
    expect(() => revertBody(['5'], {overrides: true})).toThrow(/one of/);
    expect(() => revertBody(['latest'], {})).toThrow(/version number/);
  });

  it('names paths the way operators read them', () => {
    expect(displayPath(`userGroupSettings.${ALL}.debugFlags`, POLICY)).toBe('server.debugFlags');
    expect(displayPath(`userGroupSettings.${CHEM}.saveLevels`, POLICY)).toBe('server.groups.Chemists.saveLevels');
    expect(displayPath('exportSettings.d1.endpoint', POLICY)).toBe('server.export.OpenTelemetry.endpoint');
    expect(displayPath('exportBatchSize', POLICY)).toBe('server.exportBatchSize');
  });

  it('diffs two flat maps', () => {
    expect(diffMaps({a: [1], b: 2, c: 'x'}, {a: [1], b: 3, d: []})).toEqual([
      {change: '~', path: 'b', value: '2 → 3'},
      {change: '-', path: 'c', value: 'x'},
      {change: '+', path: 'd', value: '(none)'},
    ]);
  });
});

describe('handleLogger', () => {

  it('refuses a target other than server', async () => {
    const {connect} = mockConnect(() => POLICY);
    await expect(handleLogger(connect, 'get', ['grok_connect'], {}, 'table')).rejects.toThrow(/Only 'server' is supported \(got 'grok_connect'\)/);
  });

  it('prints the policy block', async () => {
    const {connect, calls} = mockConnect(() => POLICY);
    const {out} = await captureOutput(() => handleLogger(connect, 'get', ['server'], {}, 'table'));
    expect(calls[0].path).toBe('/logging/policy');
    const label = (s: string) => s.padEnd('group Chemists'.length + 2);
    expect(out[0]).toBe(`${label('print')}error warning info          post  error warning audit usage`);
    expect(out[1]).toBe(`${label('save')}error warning audit usage   debug flags  (none)`);
    expect(out[2]).toBe(`${label('locked')}exportSettings, levels.audit   # from deployment configuration`);
    expect(out[3]).toBe(`${label('group Chemists')}save error warning audit usage info`);
    expect(out[4]).toBe(`${label('override')}package:Snowflake  debugFlags query  reverts 11:42Z  by k.lee  prod timeouts`);
  });

  it('shows effective settings for a scope with their sources', async () => {
    const {connect, calls} = mockConnect((_m, p) => p.startsWith('/logging/policy/effective')
      ? {settings: {debugFlags: ['query']}, sources: {debugFlags: 'override:o1'}} : POLICY);
    const {out} = await captureOutput(() => handleLogger(connect, 'get', [], {scope: 'package:Snowflake'}, 'table'));
    expect(calls[1].path).toBe('/logging/policy/effective?package=Snowflake');
    expect(out[2].trimEnd()).toBe('debugFlags  query  override:o1');
  });

  it('makes a package-scoped override from signed flags on top of the effective settings', async () => {
    const {connect, calls} = mockConnect((m, p) => {
      if (p === '/logging/policy') return POLICY;
      if (p.startsWith('/logging/policy/effective')) return {settings: {debugFlags: ['db']}};
      return {id: 'o2', expiresAt: utcIso(11, 42)};
    });
    const {out} = await captureOutput(() => handleLogger(connect, 'set', ['server'],
      {'debug-flags': '+queries', scope: 'package:Snowflake', for: '30m', reason: 'ELN'}, 'table'));
    expect(calls[2]).toMatchObject({method: 'POST', path: '/logging/policy/overrides',
      body: {scope: 'package', scopeId: 'Snowflake', set: {debugFlags: ['db', 'query']}, reason: 'ELN', forMinutes: 30}});
    expect(out).toEqual(['+ server.debugFlags  db, query  scope package:Snowflake  reverts 11:42Z']);
  });

  it('changes the base settings of All Users by path', async () => {
    const {connect, calls} = mockConnect((m) => m === 'PUT' ? {version: 8, changed: []} : POLICY);
    const {out} = await captureOutput(() => handleLogger(connect, 'set', [], {'save-levels': '+info', reason: 'more'}, 'table'));
    expect(calls[1]).toMatchObject({method: 'PUT', path: '/logging/policy',
      body: {set: {[`userGroupSettings.${ALL}.saveLevels`]: ['error', 'warning', 'audit', 'usage', 'info']}, reason: 'more'}});
    expect(out).toEqual(['~ server.saveLevels  error, warning, audit, usage, info', 'policy version 8']);
  });

  it('refuses --set with an override', async () => {
    const {connect} = mockConnect(() => POLICY);
    await expect(handleLogger(connect, 'set', [], {set: 'exportBatchSize=10', for: '1h'}, 'table')).rejects.toThrow(/--set changes the base settings/);
  });

  it('passes a lock refusal on as the server words it, without the HTTP status', async () => {
    const {connect} = mockConnect((m) => {
      if (m === 'PUT') throw apiError(409, 'levels.audit is locked by deployment configuration');
      return POLICY;
    });
    const refusal = await handleLogger(connect, 'set', ['server'], {'save-levels': '-audit'}, 'table').catch((e) => e);
    expect(refusal.apiError).toMatchObject({errorCode: 409, verbatim: true});
    const {err} = await captureOutput(async () => printError(refusal));
    expect(err).toEqual(['levels.audit is locked by deployment configuration']);
  });

  it('diffs against the deployment defaults, with active overrides', async () => {
    const {connect, calls} = mockConnect(() => POLICY);
    const {out} = await captureOutput(() => handleLogger(connect, 'diff', [], {}, 'table'));
    expect(calls[0].path).toBe('/logging/policy?defaults=true');
    expect(out.map((l) => l.replace(/\s+/g, ' '))).toEqual([
      '+ server.export.OpenTelemetry.endpoint https://otel.acme.internal:4318',
      '+ server.groups.Chemists.saveLevels error, warning, audit, usage, info',
      '+ server.debugFlags query scope package:Snowflake reverts 11:42Z',
    ]);
  });

  it('diffs a history version and two hosts', async () => {
    const old = {...POLICY, settings: {...POLICY.settings, [`userGroupSettings.${ALL}.debugFlags`]: ['db']}};
    const byVersion = mockConnect((_m, p) => p.includes('version=3') ? old : POLICY);
    const v = await captureOutput(() => handleLogger(byVersion.connect, 'diff', [], {version: '3'}, 'table'));
    expect(v.out.map((l) => l.replace(/\s+/g, ' ').trim())).toEqual(['~ server.debugFlags db → (none)']);

    const byHost = mockConnect((_m, _p, _b, host) => host === 'prod' ? POLICY : {...old, groups: {x1: 'All users', [CHEM]: 'Chemists'},
      settings: Object.fromEntries(Object.entries(old.settings).map(([k, val]) => [k.replace(ALL, 'x1'), val]))});
    const h = await captureOutput(() => handleLogger(byHost.connect, 'diff', [], {host: ['prod', 'val']}, 'table'));
    expect(h.out.map((l) => l.replace(/\s+/g, ' ').trim())).toEqual(['~ server.debugFlags (none) → db']);
  });

  it('lists overrides from several hosts and reverts', async () => {
    const {connect, calls} = mockConnect((m, p) => p === '/logging/policy/revert' ? {reverted: 'override o1'} : POLICY.overrides);
    const list = await captureOutput(() => handleLogger(connect, 'overrides', [], {host: ['prod', 'val']}, 'table'));
    expect(list.out[0].trimEnd().split(/\s{2,}/)).toEqual(['HOST', 'ID', 'SCOPE', 'CHANGES', 'REVERTS', 'BY', 'REASON']);
    expect(list.out[2]).toMatch(/^prod\s+o1\s+package:Snowflake\s+debugFlags=query\s+11:42Z\s+k\.lee\s+prod timeouts/);
    const rev = await captureOutput(() => handleLogger(connect, 'revert', [], {override: 'o1', reason: 'done'}, 'table'));
    expect(calls[calls.length - 1]).toMatchObject({method: 'POST', path: '/logging/policy/revert', body: {override: 'o1', reason: 'done'}});
    expect(rev.out).toEqual(['reverted override o1']);
  });
});
