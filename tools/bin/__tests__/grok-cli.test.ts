import {describe, it, expect, beforeAll, afterAll} from 'vitest';
import {spawnSync} from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';

const grokJs = path.join(__dirname, '..', 'grok.js');
const stub = path.join(__dirname, 'fixtures', 'stub-commands.cjs');

const commandNames = [...fs.readFileSync(grokJs, 'utf-8').matchAll(/^\s+'?([\w-]+)'?: lazy\(/gm)].map((m) => m[1]);

function parsedOptions(...args: string[]): Record<string, unknown> {
  const run = spawnSync(process.execPath, ['-r', stub, grokJs, ...args], {encoding: 'utf-8'});
  expect(run.status, run.stderr).toBe(0);
  return JSON.parse(run.stdout);
}

describe('grok.js argument parsing', () => {
  it('finds the command table', () => {
    expect(commandNames).toEqual(expect.arrayContaining(['api', 'check', 'publish', 'test', 'server', 's']));
  });

  // Commands such as `grok api` and `grok create` reject any option they do not expect, so a
  // default that grok.js adds for one command breaks every other one (6.7.7: `send-email: false`).
  const defaults: Record<string, Record<string, unknown>> = {
    test: {dartium: false},
    stresstest: {dartium: false},
    server: {'send-email': false},
    s: {'send-email': false},
  };

  it.each(commandNames)('`grok %s` receives no options it was not given', (command) => {
    expect(parsedOptions(command)).toEqual(defaults[command] ?? {});
  });

  it('passes --send-email and --prefix to grok s', () => {
    expect(parsedOptions('s', 'users', 'invite', 'a@b.c', '--send-email', '--prefix', '007'))
      .toEqual({'send-email': true, 'prefix': '007'});
  });

  it('keeps --send-email from swallowing the next argument', () => {
    const run = spawnSync(process.execPath, ['-r', stub, grokJs, 's', 'users', 'invite', '--send-email', 'a@b.c'],
      {encoding: 'utf-8'});
    expect(JSON.parse(run.stdout)).toEqual({'send-email': true});
  });

  it('passes --dartium to grok test', () => {
    expect(parsedOptions('test', '--dartium')).toEqual({dartium: true});
  });
});

describe('grok api in a package directory', () => {
  let dir: string;

  beforeAll(() => {
    dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-api-'));
    fs.writeFileSync(path.join(dir, 'package.json'), JSON.stringify({name: 'cli-probe', version: '0.0.1'}));
    fs.mkdirSync(path.join(dir, 'src'));
    fs.writeFileSync(path.join(dir, 'src', 'package.ts'), 'export const _package = null;\n');
  });

  afterAll(() => fs.rmSync(dir, {recursive: true, force: true}));

  it('exits with code 0', () => {
    const run = spawnSync(process.execPath, [grokJs, 'api'], {cwd: dir, encoding: 'utf-8'});
    expect(run.status, run.stdout + run.stderr).toBe(0);
  });
});
