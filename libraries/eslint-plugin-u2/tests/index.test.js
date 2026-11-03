'use strict';

const {test} = require('node:test');
const assert = require('node:assert/strict');
const {readdirSync} = require('node:fs');
const {join} = require('node:path');
const {Linter} = require('eslint');
const parser = require('@typescript-eslint/parser');
const plugin = require('..');

test('every rule file is registered, tested and configured', () => {
  const files = readdirSync(join(__dirname, '..', 'rules')).map((f) => f.replace(/\.js$/, ''));
  assert.deepEqual(Object.keys(plugin.rules).sort(), files.sort());
  const tests = readdirSync(__dirname);
  for (const rule of files) {
    assert.ok(tests.includes(`${rule}.test.js`), `${rule} has no RuleTester suite`);
    const modes = Object.values(plugin.configs).map((c) => c.rules[`u2/${rule}`]).filter((m) => m !== undefined);
    assert.ok(modes.length > 0, `${rule} is in no config`);
  }
  for (const config of Object.values(plugin.configs)) {
    for (const [rule, mode] of Object.entries(config.rules)) {
      assert.ok(rule.slice('u2/'.length) in plugin.rules, rule);
      assert.ok(mode === 'warn' || mode === 'error', `${rule}: ${mode}`);
    }
  }
});

test('the consumer config reports a HitTriage-shaped file through the flat Linter', () => {
  const code = `import {Control, Dialog, Scope} from '@datagrok-libraries/u2';
import {host} from '@datagrok-libraries/u2/src/dg/index.js';
export class HitHome extends Control {}
export class ProgressTracker { private readonly dialog = new Dialog('Progress Tracker'); }
export function init(view, o) {
  const home = Scope.runWith(undefined, () => new HitHome(o));
  view.root.replaceChildren(host(home));
  view.subs.push(grok.events.onCurrentViewChanged.subscribe(() => home.refresh()));
  return view.root.querySelector('.u2-catalog-table');
}`;
  const config = [{files: ['**/*.ts'], languageOptions: {parser}, ...plugin.configs.consumer}];
  const found = new Linter({configType: 'flat'}).verify(code, config, {filename: join(process.cwd(), 'home.ts')})
    .map((m) => m.ruleId).sort();
  assert.deepEqual(found, ['u2/import-entry', 'u2/no-internal-selectors', 'u2/no-raw-subscribe', 'u2/no-subclass',
    'u2/no-ui-in-class', 'u2/no-unowned-root']);
});
