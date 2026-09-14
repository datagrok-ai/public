import {describe, it, expect, vi, afterEach} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {add} from '../commands/add';

const testsDir = path.dirname(fileURLToPath(import.meta.url));
const toolsDir = path.dirname(path.dirname(testsDir));
const fixtureDir = path.join(testsDir, 'fixtures', 'domain-package');
const manifestPath = path.join(fixtureDir, 'databases', 'testdb', 'schema.json');
const u2Module = '@datagrok-libraries/u2/src/dg/index.js';

const cwd = process.cwd();
afterEach(() => process.chdir(cwd));

/** A scratch package: [withManifest] brings the fixture's `databases/testdb/schema.json`
 * (a package declaring its own schema), [ts] its tsconfig (a TypeScript package),
 * [template] the `grok create --ts` webpack config and tsconfig the u2 wiring edits. */
function makePackage(options?: {withManifest?: boolean, ts?: boolean, template?: boolean}): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-add-domain-app-'));
  if (options?.withManifest)
    fs.cpSync(fixtureDir, dir, {recursive: true});
  else
    fs.writeFileSync(path.join(dir, 'package.json'),
      '{"name": "@datagrok/scratch", "version": "1.0.0"}');
  if (options?.template) {
    fs.copyFileSync(path.join(toolsDir, 'package-template', 'ts.webpack.config.js'), path.join(dir, 'webpack.config.js'));
    fs.copyFileSync(path.join(toolsDir, 'package-template', 'tsconfig.json'), path.join(dir, 'tsconfig.json'));
  } else if (options?.ts !== false)
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
const read = (dir: string, ...parts: string[]) => fs.readFileSync(path.join(dir, ...parts), 'utf8');
const json = (dir: string, ...parts: string[]) => JSON.parse(read(dir, ...parts));

describe('grok add app --domain', () => {
  it('scaffolds the three-line app, the spec and a starter schema over one table of a new schema', () => {
    const dir = makePackage();
    expect(addApp('tracker.item').result).toBe(true);

    const code = entry(dir);
    expect(code).toContain(`import {domains} from '${u2Module}';`);
    expect(code).toContain(`import '@datagrok-libraries/u2/css/tokens.css';`);
    expect(code).toContain(`import '@datagrok-libraries/u2/css/domain.css';`);
    expect(code).toContain('//name: Item');
    expect(code).toContain('//tags: app');
    expect(code).toContain('//output: view result');
    expect(code).toContain('export async function Item(): Promise<DG.ViewBase> {');
    expect(code).toContain(`  return (await domains.table('tracker.item')).app();`);
    // the imports go with the template's own imports (webpack externals depend on them),
    // after the last of them and before anything else
    expect(code.indexOf(`export * from './package.g';`)).toBeLessThan(code.indexOf(u2Module));
    expect(code.indexOf(u2Module)).toBeLessThan(code.indexOf('export const _package'));

    const packageObj = json(dir, 'package.json');
    expect(packageObj.dependencies['@datagrok-libraries/u2']).toBe('../../libraries/u2');
    expect(packageObj.devDependencies['css-loader']).toBeTruthy();
    expect(packageObj.devDependencies['style-loader']).toBeTruthy();

    // the schema is the app's own: a one-table starter, and the typed clients + handles over it
    const manifest = json(dir, 'databases', 'tracker', 'schema.json');
    expect(manifest.name).toBe('tracker');
    expect(Object.keys(manifest.tables)).toEqual(['item']);
    expect(manifest.tables.item.columns.name.isName).toBe(true);
    expect(read(dir, 'src', 'generated', 'db.ts')).toContain('export interface ItemRow {');
    expect(read(dir, 'src', 'generated', 'db-ui.ts')).toContain('export function getTrackerDb(): Promise<TrackerDb>');

    // the same app as a designer-editable spec
    const spec = json(dir, 'src', 'app.spec.json');
    expect(spec.$schema).toBe('dg-ui/1');
    const source = spec.root.children[0];
    expect(source.tag).toBe('u2-domain-source');
    expect(source.name).toBe('items');
    expect(source.props.table).toBe('tracker.item');
    expect(read(dir, 'src', 'app.spec.json')).toContain('"$.items.source"');
    expect(read(dir, 'src', 'app.spec.json')).toContain('"cmd:items.save"');
  });

  it('wires webpack and tsconfig for u2 the way the Stockroom reference is wired', () => {
    const dir = makePackage({template: true});
    expect(addApp('tracker.item').result).toBe(true);

    const webpack = read(dir, 'webpack.config.js');
    expect(webpack).toContain(`alias: {'@datagrok-libraries/u2': path.resolve(__dirname, '../../libraries/u2')},`);
    expect(webpack).toContain(`{test: /\\.css$/i, use: ['style-loader', 'css-loader']},`);
    expect(webpack).toContain(`'datagrok-api/u2core': 'DG.U2',`);
    // inside the blocks they belong to
    expect(webpack.indexOf('resolve: {')).toBeLessThan(webpack.indexOf('alias:'));
    expect(webpack.indexOf('alias:')).toBeLessThan(webpack.indexOf('module: {'));
    expect(webpack.indexOf('rules: [')).toBeLessThan(webpack.indexOf('css-loader'));
    expect(webpack.indexOf('css-loader')).toBeLessThan(webpack.indexOf('plugins: ['));

    const tsconfig = read(dir, 'tsconfig.json');
    expect(tsconfig).toContain('"lib": ["ES2022", "ESNext.Disposable", "dom"],');
    expect(tsconfig).toContain('"@datagrok-libraries/u2": ["../../libraries/u2"],');
    expect(tsconfig).toContain('"@datagrok-libraries/u2/*": ["../../libraries/u2/*"],');
    expect(tsconfig).toContain('"datagrok-api/u2core": ["./node_modules/datagrok-api/src/u2core/index"]');
    expect(tsconfig.indexOf('"moduleResolution"')).toBeLessThan(tsconfig.indexOf('"paths"'));

    // rerun: every edit is made once
    expect(addApp('tracker.item').result).toBe(true);
    expect(read(dir, 'webpack.config.js').match(/alias:/g)).toHaveLength(1);
    expect(read(dir, 'webpack.config.js').match(/css-loader/g)).toHaveLength(1);
    expect(read(dir, 'tsconfig.json').match(/^\s*"paths"/gm)).toHaveLength(1);
    expect(read(dir, 'tsconfig.json').match(/ESNext\.Disposable/g)).toHaveLength(1);
  });

  it('names the app when told to', () => {
    const dir = makePackage();
    expect(addApp('tracker.item', 'Tracker').result).toBe(true);
    expect(entry(dir)).toContain('export async function Tracker(): Promise<DG.ViewBase> {');
    expect(entry(dir)).toContain('//name: Tracker');
  });

  it('a JavaScript package gets the untyped template', () => {
    const dir = makePackage({ts: false});
    expect(addApp('tracker.item').result).toBe(true);
    expect(entry(dir, false)).toContain('export async function Item() {');
    expect(entry(dir, false)).toContain(`  return (await domains.table('tracker.item')).app();`);
  });

  it('a whole schema the package declares: one app per table, plus db.ts and db-ui.ts', () => {
    const dir = makePackage({withManifest: true});
    expect(addApp('testdb').result).toBe(true);

    const code = entry(dir);
    expect(code).toContain(`return (await domains.table('testdb.sample')).app();`);
    expect(code).toContain(`return (await domains.table('testdb.sample_event')).app();`);
    expect(code).toContain('export async function Sample(');
    expect(code).toContain('export async function SampleEvent(');
    // the typed clients AND the typed u2 handles the app code can switch to
    expect(fs.existsSync(path.join(dir, 'src', 'generated', 'db.ts'))).toBe(true);
    expect(read(dir, 'src', 'generated', 'db-ui.ts')).toContain('export function getTestdbDb(): Promise<TestdbDb>');
    // the spec is over the first table
    expect(json(dir, 'src', 'app.spec.json').root.children[0].props.table).toBe('testdb.sample');
  });

  it('a table of a schema the package declares keeps the manifest as it is', () => {
    const dir = makePackage({withManifest: true});
    const before = read(dir, 'databases', 'testdb', 'schema.json');
    expect(addApp('testdb.sample').result).toBe(true);
    expect(read(dir, 'databases', 'testdb', 'schema.json')).toBe(before);
    expect(entry(dir)).not.toContain('sample_event');
  });

  it('copies a schema manifest given by path into databases/<schema>/', () => {
    const dir = makePackage();
    expect(addApp(manifestPath).result).toBe(true);
    expect(json(dir, 'databases', 'testdb', 'schema.json').name).toBe('testdb');
    expect(entry(dir)).toContain(`domains.table('testdb.sample')`);
    expect(fs.existsSync(path.join(dir, 'src', 'generated', 'db-ui.ts'))).toBe(true);
  });

  it('is rerunnable: an app that is already declared is skipped, not duplicated', () => {
    const dir = makePackage();
    expect(addApp('tracker.item').result).toBe(true);
    const spec = read(dir, 'src', 'app.spec.json');
    const manifest = read(dir, 'databases', 'tracker', 'schema.json');
    expect(addApp('tracker.item').result).toBe(true);
    expect(entry(dir).match(/export async function Item\(/g)).toHaveLength(1);
    expect(entry(dir).match(/@datagrok-libraries\/u2\/src\/dg/g)).toHaveLength(1);
    expect(entry(dir).match(/css\/tokens\.css/g)).toHaveLength(1);
    expect(read(dir, 'src', 'app.spec.json')).toBe(spec);
    expect(read(dir, 'databases', 'tracker', 'schema.json')).toBe(manifest);
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
    expect(fs.existsSync(path.join(dir, 'databases'))).toBe(false);
  });

  it('refuses to name one app for a schema of several tables', () => {
    makePackage({withManifest: true});
    expect(addApp('testdb', 'Tracker').result).toBeFalsy();
  });
});
