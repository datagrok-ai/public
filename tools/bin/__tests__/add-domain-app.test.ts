import {describe, it, expect, vi, afterEach} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {add} from '../commands/add';

const testsDir = path.dirname(fileURLToPath(import.meta.url));
const fixtureDir = path.join(testsDir, 'fixtures', 'domain-package');
const manifestPath = path.join(fixtureDir, 'databases', 'testdb', 'schema.json');

const cwd = process.cwd();
afterEach(() => process.chdir(cwd));

/** A scratch package: [withManifest] brings the fixture's `databases/testdb/schema.json`
 * (a package declaring its own schema), [ts] its tsconfig (a TypeScript package). */
function makePackage(options?: {withManifest?: boolean, ts?: boolean}): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-add-domain-app-'));
  if (options?.withManifest)
    fs.cpSync(fixtureDir, dir, {recursive: true});
  else
    fs.writeFileSync(path.join(dir, 'package.json'),
      '{"name": "@datagrok/scratch", "version": "1.0.0"}');
  if (options?.ts !== false)
    fs.writeFileSync(path.join(dir, 'tsconfig.json'), '{"compilerOptions": {"strict": true}}');
  process.chdir(dir);
  return dir;
}

/** Runs the command with console output captured, so message tests stay quiet. */
function run(args: {_: string[], domain?: string | boolean}): {result: any, output: string} {
  const spy = vi.spyOn(console, 'log').mockImplementation(() => {});
  try {
    return {result: add(args), output: spy.mock.calls.map((a) => a.join(' ')).join('\n')};
  } finally {
    spy.mockRestore();
  }
}

const addApp = (domain: string | boolean, name?: string) =>
  run({_: name == null ? ['add', 'app'] : ['add', 'app', name], domain: domain});

const entry = (dir: string, ts: boolean = true) =>
  fs.readFileSync(path.join(dir, 'src', ts ? 'package.ts' : 'package.js'), 'utf8');

describe('grok add app --domain', () => {
  it('scaffolds a one-expression app over one table of a schema the package does not declare', () => {
    const dir = makePackage();
    expect(addApp('apitests.item').result).toBe(true);

    const code = entry(dir);
    expect(code).toContain(`import {domains} from '@datagrok-libraries/domain-ui';`);
    expect(code).toContain('//name: Item');
    expect(code).toContain('//meta.role: app');
    expect(code).toContain('//input: string path {meta.url: true; optional: true}');
    expect(code).toContain('//output: view result');
    expect(code).toContain('export async function Item(): Promise<DG.ViewBase> {');
    expect(code).toContain(`  return (await domains.table('apitests.item')).app();`);
    // the import goes with the template's own imports (webpack externals depend on them),
    // after the last of them and before anything else
    expect(code.indexOf(`export * from './package.g';`))
      .toBeLessThan(code.indexOf('@datagrok-libraries/domain-ui'));
    expect(code.indexOf('@datagrok-libraries/domain-ui'))
      .toBeLessThan(code.indexOf('export const _package'));

    const packageObj = JSON.parse(fs.readFileSync(path.join(dir, 'package.json'), 'utf8'));
    expect(packageObj.dependencies['@datagrok-libraries/domain-ui']).toBeTruthy();
    // nothing to generate: the schema is another package's
    expect(fs.existsSync(path.join(dir, 'src', 'generated'))).toBe(false);
  });

  it('names the app when told to', () => {
    const dir = makePackage();
    expect(addApp('apitests.item', 'Tracker').result).toBe(true);
    expect(entry(dir)).toContain('export async function Tracker(): Promise<DG.ViewBase> {');
    expect(entry(dir)).toContain('//name: Tracker');
  });

  it('a JavaScript package gets the untyped template', () => {
    const dir = makePackage({ts: false});
    expect(addApp('apitests.item').result).toBe(true);
    expect(entry(dir, false)).toContain('export async function Item() {');
  });

  it('a whole schema the package declares: one app per table, plus db.ts and db-ui.ts', () => {
    const dir = makePackage({withManifest: true});
    expect(addApp('testdb').result).toBe(true);

    const code = entry(dir);
    expect(code).toContain(`return (await domains.table('testdb.sample')).app();`);
    expect(code).toContain(`return (await domains.table('testdb.sample_event')).app();`);
    expect(code).toContain('export async function Sample(');
    expect(code).toContain('export async function SampleEvent(');
    // the typed clients AND the typed UI wrappers the app code can switch to
    expect(fs.existsSync(path.join(dir, 'src', 'generated', 'db.ts'))).toBe(true);
    expect(fs.readFileSync(path.join(dir, 'src', 'generated', 'db-ui.ts'), 'utf8'))
      .toContain('export const sampleUi = new SampleUi();');
  });

  it('copies a schema manifest given by path into databases/<schema>/', () => {
    const dir = makePackage();
    expect(addApp(manifestPath).result).toBe(true);
    expect(JSON.parse(fs.readFileSync(path.join(dir, 'databases', 'testdb', 'schema.json'), 'utf8')).name)
      .toBe('testdb');
    expect(entry(dir)).toContain(`domains.table('testdb.sample')`);
    expect(fs.existsSync(path.join(dir, 'src', 'generated', 'db-ui.ts'))).toBe(true);
  });

  it('is rerunnable: an app that is already declared is skipped, not duplicated', () => {
    const dir = makePackage();
    expect(addApp('apitests.item').result).toBe(true);
    expect(addApp('apitests.item').result).toBe(true);
    expect(entry(dir).match(/export async function Item\(/g)).toHaveLength(1);
    expect(entry(dir).match(/@datagrok-libraries\/domain-ui/g)).toHaveLength(1);
  });

  it('rejects --domain on anything but an app, and a malformed address', () => {
    makePackage();
    expect(run({_: ['add', 'function'], domain: 'testdb.sample'}).result).toBeFalsy();
    expect(addApp('testdb.sample.column').result).toBeFalsy();
    expect(addApp('not a schema').result).toBeFalsy();
    expect(addApp(true).result).toBeFalsy();
  });

  it('a whole schema the package does not declare names no tables — say so', () => {
    const dir = makePackage();
    expect(addApp('apitests').result).toBeFalsy();
    expect(fs.existsSync(path.join(dir, 'src'))).toBe(false);
  });

  it('refuses to name one app for a schema of several tables', () => {
    makePackage({withManifest: true});
    expect(addApp('testdb', 'Tracker').result).toBeFalsy();
  });
});
