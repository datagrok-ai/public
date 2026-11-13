import {afterAll, beforeAll, describe, expect, it} from 'vitest';
import {spawnSync} from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {check} from '../commands/check';

const u2 = path.resolve(__dirname, '../../../libraries/u2');
const OWN = '.zz-css-home {\n  --u2-page-gap: 8px;\n}\n';
const INTERNAL = '\n.zz-css-home > .u2-x {\n  color: red;\n}\n';
let dir: string;

const css = (text: string) => fs.writeFileSync(path.join(dir, 'css', 'zz-css.css'), text);
const grokCheck = () => check({_: ['check', dir], 'no-exit': true} as Parameters<typeof check>[0]);
const consumerCheck = () => {
  const run = spawnSync(process.execPath, [path.join(u2, 'scripts', 'check-consumer.mjs'), dir, '--json'],
    {encoding: 'utf-8', timeout: 120000});
  return JSON.parse(run.stdout) as {ok: boolean, trees: {problems: string[]}[]};
};

// a pre-1.0 package outside the repository that reaches u2 through node_modules, as a workspace package does
describe('grok check: u2 internal selectors in a package\'s CSS', {timeout: 60_000}, () => {
  beforeAll(() => {
    dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-check-u2-'));
    fs.writeFileSync(path.join(dir, 'package.json'), JSON.stringify({name: '@datagrok/zz-css', version: '0.0.1',
      description: 'u2 CSS check fixture', dependencies: {'@datagrok-libraries/u2': 'workspace:^'}}));
    fs.mkdirSync(path.join(dir, 'css'));
    fs.mkdirSync(path.join(dir, 'node_modules', '@datagrok-libraries'), {recursive: true});
    fs.symlinkSync(u2, path.join(dir, 'node_modules', '@datagrok-libraries', 'u2'), 'junction');
  });

  afterAll(() => fs.rmSync(dir, {recursive: true, force: true}));

  it('passes a package that styles its own classes and sets --u2-* properties', () => {
    css(OWN);
    expect(consumerCheck().ok).toBe(true);
    expect(grokCheck()).toBe(true);
  });

  it('fails the same package once its CSS targets a .u2-* class', () => {
    css(OWN + INTERNAL);
    const report = consumerCheck();
    expect(report.ok).toBe(false);
    expect(report.trees[0].problems).toHaveLength(1);
    expect(report.trees[0].problems[0]).toMatch(/zz-css\.css:5:1 {2}css-internal-selectors {2}`\.zz-css-home > \.u2-x \{`/);
    expect(grokCheck()).toBe(false);
  });

  it('passes it again once the rule is gone', () => {
    css(OWN);
    expect(grokCheck()).toBe(true);
  });
});
