import {describe, it, expect} from 'vitest';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import {handleRules, ruleBody, ruleRow, testRows} from '../commands/server-rules';
import {mockConnect, captureOutput, utcIso} from './obs-helpers';

const DEFINITION = {severity: 'warning', match: {source: 'audit', type: 'user-login-failed'}, groupBy: 'param:user',
  window: 15, when: {count: 5}, summary: '{count} failed logins for {group} in {window} min'};

const RULE = {name: 'failed-logins-per-user', kind: 'rule-failed-logins-per-user', source: 'settings', readOnly: false,
  enabled: true, shadowed: false, definition: DEFINITION, error: null, problems: {ongoing: 2, alerting: 1},
  createdBy: 'admin', createdAt: utcIso(9, 0), updatedBy: 'admin', updatedAt: utcIso(9, 30)};

const TEST = {rule: DEFINITION, kind: 'rule-failed-logins-per-user', at: utcIso(10, 0), ms: 42, matched: 17,
  now: [{key: 'alice', group: 'alice', severity: 'warning', summary: '6 failed logins for alice in 15 min', count: 6}],
  raised: [{key: 'bob', group: 'bob', firstAt: utcIso(4, 0), lastAt: utcIso(5, 0), summary: '5 failed logins for bob in 15 min'}],
  points: 25, truncated: false};

describe('ruleBody', () => {
  const file = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'rules-')), 'rule.json');

  it('reads a file or inline JSON; a name argument sets the name', () => {
    fs.writeFileSync(file, JSON.stringify({name: 'old', ...DEFINITION}));
    expect(ruleBody({json: file})).toEqual({name: 'old', ...DEFINITION});
    expect(ruleBody({json: file}, 'new-name').name).toBe('new-name');
    expect(ruleBody({json: ' {"window": 5}'}, 'x')).toEqual({window: 5, name: 'x'});
  });

  it('refuses no rule, invalid JSON and a non-object', () => {
    expect(() => ruleBody({})).toThrow('Pass the rule as --json');
    expect(() => ruleBody({json: true})).toThrow('Pass the rule as --json');
    expect(() => ruleBody({json: '{"window": }'})).toThrow(/^Invalid rule JSON in --json: /);
    fs.writeFileSync(file, '[1]');
    expect(() => ruleBody({json: file})).toThrow(`Invalid rule JSON in ${file}: a rule is one JSON object`);
    expect(() => ruleBody({json: path.join(path.dirname(file), 'missing.json')})).toThrow(/^Invalid rule JSON in .*missing\.json/);
  });
});

describe('ruleRow and testRows', () => {
  it('prints the rule columns, count and warning by default', () => {
    expect(ruleRow(RULE)).toEqual({NAME: 'failed-logins-per-user', SOURCE: 'settings', ON: 'yes', WHEN: 'count',
      GROUP: 'param:user', SEV: 'warning', ONGOING: 2, ALERTING: 1, ERROR: ''});
    const bare = ruleRow({name: 'x', enabled: false, definition: {groupBy: ['param:a', 'param:b'], when: {count: 3, users: 2}}});
    expect(bare).toMatchObject({ON: 'no', WHEN: 'count+users', GROUP: 'param:a,param:b', SEV: 'warning', ONGOING: 0});
    expect(ruleRow({definition: {}}).WHEN).toBe('count');
  });

  it('shapes what holds now and what would have raised', () => {
    const t = testRows(TEST);
    expect(t.now).toEqual([{GROUP: 'alice', SUMMARY: '6 failed logins for alice in 15 min', COUNT: 6}]);
    expect(t.raised).toEqual([{FIRST: '04:00Z', GROUP: 'bob', SUMMARY: '5 failed logins for bob in 15 min'}]);
    expect(testRows({})).toEqual({now: [], raised: []});
  });
});

describe('handleRules', () => {
  it('lists rules and gets one with its definition as JSON', async () => {
    const {connect, calls} = mockConnect((_m, p) => p === '/problems/rules' ? [RULE] : RULE);
    const list = await captureOutput(() => handleRules(connect, 'list', [], {}, 'table'));
    expect(list.out[0].trimEnd().split(/\s{2,}/)).toEqual(['NAME', 'SOURCE', 'ON', 'WHEN', 'GROUP', 'SEV', 'ONGOING', 'ALERTING', 'ERROR']);
    const get = await captureOutput(() => handleRules(connect, 'get', [RULE.name], {}, 'table'));
    expect(calls.map((c) => c.path)).toEqual(['/problems/rules', `/problems/rules/${RULE.name}`]);
    expect(get.out[0]).toMatch(/^rule\s+failed-logins-per-user  rule-failed-logins-per-user$/);
    expect(JSON.parse(get.out.slice(get.out.indexOf('') + 1).join('\n'))).toEqual(DEFINITION);
  });

  it('adds, edits, enables, disables and deletes', async () => {
    const {connect, calls} = mockConnect((m, p) => m === 'DELETE' ? {deleted: RULE.name} : {...RULE, enabled: !p.endsWith('/disable')});
    const json = JSON.stringify(DEFINITION);
    const out: string[] = [];
    for (const [verb, args] of [['add', [RULE.name]], ['edit', [RULE.name]], ['disable', [RULE.name]], ['enable', [RULE.name]],
      ['delete', [RULE.name]]] as [string, string[]][])
      out.push(...(await captureOutput(() => handleRules(connect, verb, args, {json}, 'table'))).out);
    expect(calls.map((c) => [c.method, c.path, c.body])).toEqual([
      ['POST', '/problems/rules', {...DEFINITION, name: RULE.name}],
      ['PUT', `/problems/rules/${RULE.name}`, DEFINITION],
      ['POST', `/problems/rules/${RULE.name}/disable`, {}],
      ['POST', `/problems/rules/${RULE.name}/enable`, {}],
      ['DELETE', `/problems/rules/${RULE.name}`, undefined],
    ]);
    expect(out).toEqual([`added rule ${RULE.name} (enabled)`, `replaced the definition of ${RULE.name}`,
      `disabled ${RULE.name}`, `enabled ${RULE.name}`, `deleted ${RULE.name}`]);
  });

  it('tests a stored rule by name or a definition, and prints both tables', async () => {
    const {connect, calls} = mockConnect(() => TEST);
    const {out} = await captureOutput(() => handleRules(connect, 'test', [RULE.name], {hours: 6}, 'table'));
    await captureOutput(() => handleRules(connect, 'test', [], {json: JSON.stringify({name: 'x', ...DEFINITION})}, 'json'));
    expect(calls.map((c) => [c.path, c.body])).toEqual([['/problems/rules/test?hours=6', {name: RULE.name}],
      ['/problems/rules/test', {name: 'x', ...DEFINITION}]]);
    expect(out[0]).toBe('Matched 17 events in 6 h (42 ms)');
    expect(out).toContain('Holds now');
    expect(out).toContain('Would have raised in the last 6 h');
    expect(out.join('\n')).toMatch(/alice\s+6 failed logins for alice in 15 min\s+6/);
    expect(out.join('\n')).toMatch(/04:00Z\s+bob/);
  });

  it('answers a missing name or rule, and an unknown verb, with the usage', async () => {
    const {connect, calls} = mockConnect(() => ({}));
    for (const [verb, line] of [['get', 'get <name>'], ['edit', 'edit <name>'], ['enable', 'enable <name>'],
      ['delete', 'delete <name>'], ['test', 'test (<name> | --json']]) {
      const {err, result} = await captureOutput(() => handleRules(connect, verb, [], {}, 'table'));
      expect(result).toBe(false);
      expect(err.join('\n')).toContain(`Usage: grok s observe rules ${line}`);
    }
    await expect(handleRules(connect, 'add', [], {}, 'table')).rejects.toThrow('Pass the rule as --json');
    const {err, result} = await captureOutput(() => handleRules(connect, 'rename', [], {}, 'table'));
    expect(result).toBe(false);
    expect(err.join('\n')).toMatch(/Usage: grok s observe rules <verb>/);
    expect(calls).toEqual([]);
  });
});
