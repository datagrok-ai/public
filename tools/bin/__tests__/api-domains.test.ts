import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import ts from 'typescript';
import {createRequire} from 'module';
import {fileURLToPath} from 'url';
import {spawnSync} from 'child_process';
import {generateDomainClients} from '../commands/api';

const testsDir = path.dirname(fileURLToPath(import.meta.url));
const fixtureDir = path.join(testsDir, 'fixtures', 'domain-package');
const toolsDir = path.dirname(path.dirname(testsDir));
const jsApiDir = path.resolve(toolsDir, '..', 'js-api');
const domainUiDir = path.resolve(toolsDir, '..', 'libraries', 'domain-ui');
const dbPath = (dir: string) => path.join(dir, 'src', 'generated', 'db.ts');
const dbUiPath = (dir: string) => path.join(dir, 'src', 'generated', 'db-ui.ts');

function makePackage(): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-api-domains-'));
  fs.cpSync(fixtureDir, dir, {recursive: true});
  return dir;
}

function mutateManifest(dir: string, mutate: (manifest: any) => void): void {
  const manifestPath = path.join(dir, 'databases', 'testdb', 'schema.json');
  const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
  mutate(manifest);
  fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
}

/** Adds a `label` target table, a `sample_label` junction, and the `labels` many-to-many
 * relation on `sample` — the fixture every relation test starts from. */
function addRelation(manifest: any,
  relation: {[key: string]: string} = {via: 'sample_label', target: 'label'}): void {
  manifest.tables.label = {
    businessKey: ['name'], columns: {name: {type: 'string', required: true}}};
  manifest.tables.sample_label = {
    securityMode: 'master', delegate: 'sample_id', businessKey: ['sample_id', 'label_id'],
    columns: {
      sample_id: {type: 'ref', ref: 'sample', onDelete: 'cascade', required: true},
      label_id: {type: 'ref', ref: 'label', onDelete: 'cascade', required: true},
    }};
  manifest.tables.sample.relations = {labels: relation};
}

/** Runs the generator with console output captured, so error-message tests stay quiet. */
function runCapturingLog(dir: string): {result: boolean, output: string} {
  const spy = vi.spyOn(console, 'log').mockImplementation(() => {});
  try {
    const result = generateDomainClients(dir);
    return {result, output: spy.mock.calls.map((args) => args.join(' ')).join('\n')};
  } finally {
    spy.mockRestore();
  }
}

/** Copies a package's `.d.ts` tree (declarations only, so tsc never pulls the `.ts`
 * sources) from [src] to [dst]. */
function copyDts(src: string, dst: string): void {
  for (const entry of fs.readdirSync(src, {withFileTypes: true})) {
    if (entry.isDirectory()) {
      if (entry.name !== 'node_modules')
        copyDts(path.join(src, entry.name), path.join(dst, entry.name));
    } else if (entry.name.endsWith('.d.ts')) {
      fs.mkdirSync(dst, {recursive: true});
      fs.copyFileSync(path.join(src, entry.name), path.join(dst, entry.name));
    }
  }
}

/** Stages the local js-api types into `<dir>/node_modules/datagrok-api`. */
function stageDatagrokApiTypes(dir: string): void {
  const dstRoot = path.join(dir, 'node_modules', 'datagrok-api');
  copyDts(jsApiDir, dstRoot);
  fs.writeFileSync(path.join(dstRoot, 'package.json'),
    '{"name": "datagrok-api", "version": "0.0.0", "types": "datagrok.d.ts"}');
  // dayjs types (the generated db.ts imports `type {Dayjs}`).
  const dayjsSrc = path.join(jsApiDir, 'node_modules', 'dayjs');
  const dayjsDst = path.join(dir, 'node_modules', 'dayjs');
  fs.mkdirSync(path.join(dayjsDst, 'locale'), {recursive: true});
  fs.copyFileSync(path.join(dayjsSrc, 'package.json'), path.join(dayjsDst, 'package.json'));
  fs.copyFileSync(path.join(dayjsSrc, 'index.d.ts'), path.join(dayjsDst, 'index.d.ts'));
  fs.copyFileSync(path.join(dayjsSrc, 'locale', 'index.d.ts'),
    path.join(dayjsDst, 'locale', 'index.d.ts'));
  fs.copyFileSync(path.join(dayjsSrc, 'locale', 'types.d.ts'),
    path.join(dayjsDst, 'locale', 'types.d.ts'));
}

/** Stages the local domain-ui types (what the generated db-ui.ts imports) into
 * `<dir>/node_modules/@datagrok-libraries/domain-ui`. Its own `rxjs` imports stay
 * unresolved on purpose — `skipLibCheck` covers them, and staging 670 rxjs
 * declarations per test would not check anything this file cares about. */
function stageDomainUiTypes(dir: string): void {
  const dstRoot = path.join(dir, 'node_modules', '@datagrok-libraries', 'domain-ui');
  copyDts(path.join(domainUiDir, 'src'), path.join(dstRoot, 'src'));
  fs.writeFileSync(path.join(dstRoot, 'package.json'),
    '{"name": "@datagrok-libraries/domain-ui", "version": "0.0.0", "types": "./src/index.d.ts"}');
}

/** Typechecks the generated files plus one usage fixture against the local js-api (and,
 * for the `--ui` fixtures, domain-ui) types. */
function runTsc(dir: string, usageFile: string): {status: number | null, output: string} {
  stageDatagrokApiTypes(dir);
  const ui = fs.existsSync(dbUiPath(dir));
  if (ui)
    stageDomainUiTypes(dir);
  fs.writeFileSync(path.join(dir, 'tsconfig.json'), JSON.stringify({
    compilerOptions: {
      target: 'es2020',
      module: 'es2020',
      moduleResolution: 'bundler',
      strict: true,
      noEmit: true,
      skipLibCheck: true,
    },
    include: ['src/generated/db.ts', ...(ui ? ['src/generated/db-ui.ts'] : []), usageFile],
  }, null, 2));
  const tsc = path.join(toolsDir, 'node_modules', 'typescript', 'bin', 'tsc');
  const res = spawnSync(process.execPath, [tsc, '-p', dir], {encoding: 'utf8'});
  return {status: res.status, output: (res.stdout ?? '') + (res.stderr ?? '')};
}

describe('generateDomainClients', () => {
  it('generates typed row/insert interfaces, column unions, and schema clients', () => {
    const dir = makePackage();
    expect(generateDomainClients(dir)).toBe(true);
    const code = fs.readFileSync(dbPath(dir), 'utf8');

    // generated header and imports
    expect(code).toContain('auto-generated by the grok api command');
    expect(code).toContain(`import * as grok from 'datagrok-api/grok';`);
    expect(code).toContain(`import * as DG from 'datagrok-api/dg';`);
    expect(code).toContain(`import type {Dayjs} from 'dayjs';`);

    // row interface: system fields + relational columns + jsonb properties, optionality from `required`
    expect(code).toContain('export interface SampleRow {');
    expect(code).toContain('  id: string;');
    expect(code).toContain('  version: number;');
    expect(code).toContain('  created_on: Dayjs;');     // system datetime → dayjs
    expect(code).toContain('  name: string;');          // required
    expect(code).toContain('  count?: number;');        // int, optional
    expect(code).toContain('  ratio?: number;');        // float
    expect(code).toContain('  active?: boolean;');      // bool
    expect(code).toContain('  measured_on?: Dayjs;');   // datetime → dayjs on reads
    expect(code).toContain('  tags?: string[];');       // string_list
    expect(code).toContain('  owner?: string;');        // user
    expect(code).toContain('  note?: string;');         // jsonb property
    expect(code).toContain('  score: number;');         // required jsonb property
    expect(code).toContain('  sample_id: string;');     // ref
    expect(code).toContain('  status?: SampleStatus;'); // choices → named literal union
    expect(code).toContain(`export type SampleStatus = 'new' | 'done';`);

    // insert interface: system fields omitted; datetimes accept dayjs OR ISO strings;
    // idempotencyKey only where enabled
    expect(code).toContain('export interface SampleInsert {');
    expect(code).not.toMatch(/interface SampleInsert \{[^}]*\bid:/);
    expect(code).toContain('  measured_on?: Dayjs | string;');
    expect(code.match(/idempotencyKey\?: string;/g)).toHaveLength(1);
    expect(code).toMatch(/interface SampleInsert \{[^}]*idempotencyKey\?: string;/);
    expect(code).toContain('export interface SampleEventInsert {');
    expect(code).not.toMatch(/interface SampleEventInsert \{[^}]*idempotencyKey/);

    // column-name unions
    expect(code).toMatch(/export type SampleColumn = 'id' \| 'version' \| 'created_on' \| 'updated_on' \| 'author_id' \|/);
    expect(code).toMatch(/export type SampleEventColumn = [^;]*'sample_id' \| 'kind';/);
    expect(code).toContain(`'measured_on'`);

    // expand maps: master keys with '<fk>.<col>' fields, details keys with child-row arrays;
    // type ALIASES, not interfaces (the TExpand index-signature constraint)
    expect(code).toContain('export type SampleExpand = {');
    expect(code).toContain(`  'details:sample_event': {sample_event?: SampleEventRow[]};`);
    expect(code).toContain('export type SampleEventExpand = {');
    expect(code).toMatch(/'sample_id': \{'sample_id\.name'\?: string;/);
    expect(code).toContain(`'sample_id.measured_on'?: Dayjs;`);
    expect(code).toContain(`'sample_id.status'?: SampleStatus;`);
    expect(code).not.toContain('interface SampleExpand');

    // transaction union: one insert/update/delete arm per table
    expect(code).toContain('export type TestdbTransactionOp =');
    expect(code).toContain(
      `  {op: 'insert'; table: 'sample'; ref?: string; values: DG.DomainTxValues<SampleInsert>} |`);
    expect(code).toContain(`  {op: 'delete'; table: 'sample_event'; id: string};`);

    // per-schema typed clients: LAZY getters (no import-time side effects), four
    // generics, PLURAL property names, no datetime config (the client resolves
    // datetime columns from the registry)
    expect(code).toContain('export const testdbDb = {');
    expect(code).toContain('  get samples() {');
    expect(code).toContain('  get sampleEvents() {');
    expect(code).toContain(
      `    return grok.dapi.domains.table<SampleRow, SampleInsert, SampleColumn, SampleExpand>('testdb.sample');`);
    expect(code).not.toContain('datetimeColumns');
    expect(code).not.toMatch(/ as\r?\n?\s*DG\.DomainTableClient/);
    // mapped-tuple transaction: per-op result types for tuple ops literals
    expect(code).toContain(`  transaction<T extends TestdbTransactionOp[]>(ops: [...T]):`);
    expect(code).toContain(`      Promise<{[K in keyof T]: DG.DomainOpResultFor<T[K]>}> {`);

    // CRLF line endings per the repo code style
    expect(code).not.toMatch(/[^\r]\n/);
  });

  it('choices aliases dedupe identical value sets and disambiguate different ones', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => {
      m.tables.issue = {columns: {
        state_x: {type: 'string', choices: ['x', 'y']},
        state_y: {type: 'string', choices: ['a', 'b']},
      }};
      m.tables.issue_state = {columns: {
        x: {type: 'string', choices: ['x', 'y']},
        y: {type: 'string', choices: ['c', 'd']},
      }};
    });
    expect(generateDomainClients(dir)).toBe(true);
    const code = fs.readFileSync(dbPath(dir), 'utf8');
    // identical sets share ONE alias; different sets get the table-name suffix
    expect(code.match(/export type IssueStateX = /g)).toHaveLength(1);
    expect(code).toContain(`export type IssueStateX = 'x' | 'y';`);
    expect(code).toContain(`export type IssueStateY = 'a' | 'b';`);
    expect(code).toContain(`export type IssueStateYIssueState = 'c' | 'd';`);
    expect(code).toMatch(/interface IssueStateRow \{[^}]*x\?: IssueStateX;/);
    expect(code).toMatch(/interface IssueStateRow \{[^}]*y\?: IssueStateYIssueState;/);
  });

  it('expand maps: a multi-FK child emits one qualified details key per FK', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => m.tables.link = {columns: {
      a_id: {type: 'ref', ref: 'sample', onDelete: 'cascade'},
      b_id: {type: 'ref', ref: 'sample', onDelete: 'cascade'},
    }});
    expect(generateDomainClients(dir)).toBe(true);
    const code = fs.readFileSync(dbPath(dir), 'utf8');
    expect(code).toContain(`  'details:link.a_id': {link?: LinkRow[]};`);
    expect(code).toContain(`  'details:link.b_id': {link?: LinkRow[]};`);
    expect(code).not.toContain(`'details:link':`);
    // the child's own master expands stay per-column
    expect(code).toMatch(/'a_id': \{'a_id\.name'\?: string;/);
  });

  it('relations: expand entry on the read side, link sets on the write side', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => addRelation(m));
    expect(generateDomainClients(dir)).toBe(true);
    const code = fs.readFileSync(dbPath(dir), 'utf8');

    // read side: the link array lands under the relation's own name
    expect(code).toContain(`  'labels': {labels?: {id: string; name: string}[]};`);
    // write side: the insert payload and the update type the tx arm uses
    expect(code).toMatch(/interface SampleInsert \{[^}]*labels\?: string\[\];/);
    expect(code).toContain('export type SampleUpdate = Partial<SampleRow> & {');
    expect(code).toMatch(/export type SampleUpdate = [^;]*labels\?: string\[\];/);
    expect(code).toContain(
      `values: DG.DomainTxValues<SampleUpdate>; expectedVersion?: number} |`);
    // a relation is NOT a column: it never joins the column union
    expect(code).not.toMatch(/export type SampleColumn = [^;]*'labels'/);
    // the client takes the update type as its fifth generic, so update() accepts
    // the link sets; a relation-less table keeps the client's Partial<Row> default
    expect(code).toContain(
      `table<SampleRow, SampleInsert, SampleColumn, SampleExpand, SampleUpdate>('testdb.sample')`);
    expect(code).toContain(
      `table<LabelRow, LabelInsert, LabelColumn, LabelExpand>('testdb.label')`);
    // relation-less tables are untouched — the update arm keeps Partial<Row>
    expect(code).toContain(
      `values: DG.DomainTxValues<Partial<LabelRow>>; expectedVersion?: number} |`);
    expect(code).not.toContain('export type LabelUpdate');
    // the junction is an ordinary table and still expands as a details child
    expect(code).toContain(`  'details:sample_label': {sample_label?: SampleLabelRow[]};`);
  });

  it('relations: an ambiguous junction side must be named explicitly', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => {
      addRelation(m);
      m.tables.sample_label.columns.other_label_id = {type: 'ref', ref: 'label', onDelete: 'cascade'};
    });
    const {result, output} = runCapturingLog(dir);
    expect(result).toBe(false);
    expect(output).toContain(`has more than one ref column targeting 'label'`);
    expect(output).toContain(`declare 'viaTarget' explicitly`);
    expect(fs.existsSync(dbPath(dir))).toBe(false);

    // naming the side resolves it; a name that is not a ref to the target does not
    mutateManifest(dir, (m) => m.tables.sample.relations.labels.viaTarget = 'label_id');
    expect(generateDomainClients(dir)).toBe(true);
    mutateManifest(dir, (m) => m.tables.sample.relations.labels.viaTarget = 'sample_id');
    const bad = runCapturingLog(dir);
    expect(bad.result).toBe(false);
    expect(bad.output).toContain(
      `'sample_id' is not a ref column of junction table 'sample_label' targeting 'label'`);
  });

  it('relations: via and target must be declared in the same manifest', () => {
    const missingTarget = makePackage();
    mutateManifest(missingTarget, (m) => {
      addRelation(m);
      m.tables.sample.relations.labels.target = 'nope';
    });
    let res = runCapturingLog(missingTarget);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`target table 'nope' is not declared in this manifest`);
    expect(fs.existsSync(dbPath(missingTarget))).toBe(false);

    const missingVia = makePackage();
    mutateManifest(missingVia, (m) => {
      addRelation(m);
      m.tables.sample.relations.labels.via = 'nope';
    });
    res = runCapturingLog(missingVia);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`junction table 'nope' is not declared in this manifest`);
  });

  it('relations: the junction business key must cover both FKs', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => {
      addRelation(m);
      m.tables.sample_label.businessKey = ['sample_id'];
    });
    const {result, output} = runCapturingLog(dir);
    expect(result).toBe(false);
    expect(output).toContain(
      `must declare a 'businessKey' containing both 'sample_id' and 'label_id'`);
  });

  it('relations: a name may not collide with a column, and self-referential needs both sides', () => {
    const collision = makePackage();
    mutateManifest(collision, (m) => {
      addRelation(m);
      m.tables.sample.relations = {name: {via: 'sample_label', target: 'label'}};
    });
    let res = runCapturingLog(collision);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`relation 'sample.name' collides with a column of 'sample'`);

    // a relation back at its own table cannot auto-resolve the two sides
    const self = makePackage();
    mutateManifest(self, (m) => {
      addRelation(m);
      m.tables.sample_link = {businessKey: ['from_id', 'to_id'], columns: {
        from_id: {type: 'ref', ref: 'sample', onDelete: 'cascade', required: true},
        to_id: {type: 'ref', ref: 'sample', onDelete: 'cascade', required: true},
      }};
      m.tables.sample.relations = {blocks: {via: 'sample_link', target: 'sample'}};
    });
    res = runCapturingLog(self);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`a self-referential relation must name both 'viaSelf' and 'viaTarget'`);

    mutateManifest(self, (m) =>
      m.tables.sample.relations = {blocks:
        {via: 'sample_link', target: 'sample', viaSelf: 'from_id', viaTarget: 'to_id'}});
    expect(generateDomainClients(self)).toBe(true);
    expect(fs.readFileSync(dbPath(self), 'utf8'))
      .toContain(`  'blocks': {blocks?: {id: string; name: string}[]};`);
  });

  it('relations: the junction must differ from the owner and the target', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => {
      addRelation(m);
      m.tables.sample.relations.labels.target = 'sample_label';
    });
    const {result, output} = runCapturingLog(dir);
    expect(result).toBe(false);
    expect(output).toContain(`must differ from the owner and the target table`);
  });

  it('qualified refs: a Core target types the expand from the sealed Core declaration', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => {
      m.tables.sample.columns.source_query = {type: 'ref', ref: 'Core.queries'};
      m.tables.sample.columns.package = {type: 'ref', ref: 'Core.packages', required: true};
    });
    expect(generateDomainClients(dir)).toBe(true);
    const code = fs.readFileSync(dbPath(dir), 'utf8');
    expect(code).toContain('  source_query?: string;');
    expect(code).toContain('  package: string;');
    // the expand entry lists the target's declared columns with the same tsTypes
    expect(code).toContain(`  'source_query': {'source_query.friendlyName'?: string; 'source_query.createdOn'?: Dayjs;`);
    expect(code).toContain(`'source_query.updatedOn'?: Dayjs; 'source_query.author'?: string; ` +
      `'source_query.connection'?: string};`);
    // a Core column named like a generated system column is the target's own, not a collision
    expect(code).toMatch(/'package': \{[^}]*'package\.version'\?: string;/);
    // the external table exists only to type the expand: no accessor, types, or tx arms
    expect(code).not.toContain('QueriesRow');
    expect(code).not.toContain('PackagesRow');
    expect(code).not.toMatch(/table: '(queries|packages)'/);
    expect(code).not.toContain('Core.');
  });

  it('qualified refs: another plugin\'s table resolves from its installed manifest', () => {
    const dir = makePackage();
    const depDir = path.join(dir, 'node_modules', '@datagrok', 'other-plugin', 'databases', 'other');
    fs.mkdirSync(depDir, {recursive: true});
    fs.writeFileSync(path.join(depDir, 'schema.json'), JSON.stringify({
      name: 'other', version: '1.0.0',
      tables: {thing: {
        columns: {name: {type: 'string', required: true}, kind: {type: 'string', choices: ['a', 'b']}},
        schemas: ['extra', 'sealed'],
      }},
      // a manifest keys property-schema columns by name; a sealed declaration lists them
      propertySchemas: {extra: {weight: {type: 'float'}}, sealed: [{name: 'rank', type: 'int'}]},
    }));
    mutateManifest(dir, (m) => m.tables.sample.columns.thing_id = {type: 'ref', ref: 'other.thing'});
    expect(generateDomainClients(dir)).toBe(true);
    const code = fs.readFileSync(dbPath(dir), 'utf8');
    expect(code).toContain(`export type OtherThingKind = 'a' | 'b';`);
    expect(code).toContain(`  'thing_id': {'thing_id.name'?: string; 'thing_id.kind'?: OtherThingKind; ` +
      `'thing_id.weight'?: number;`);
    expect(code).toContain(`'thing_id.rank'?: number};`);
    expect(code).not.toContain('ThingRow');
  });

  it('qualified refs: type user/group is the same column as ref Core.users/Core.groups', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => {
      m.tables.sample.columns.owner = {type: 'user', onDelete: 'setnull'};
      m.tables.sample.columns.team = {type: 'group'};
    });
    expect(generateDomainClients(dir)).toBe(true);
    const aliased = fs.readFileSync(dbPath(dir));
    expect(aliased.toString()).toContain(`  'owner': {'owner.login'?: string; 'owner.firstName'?: string;`);
    expect(aliased.toString()).toMatch(/'team': \{'team\.friendlyName'\?: string;/);

    mutateManifest(dir, (m) => {
      m.tables.sample.columns.owner = {type: 'ref', ref: 'Core.users'};
      m.tables.sample.columns.team = {type: 'ref', ref: 'Core.groups'};
    });
    expect(generateDomainClients(dir)).toBe(true);
    expect(fs.readFileSync(dbPath(dir)).equals(aliased)).toBe(true);
  });

  it('qualified refs: an unresolvable target is an error', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => m.tables.sample.columns.x = {type: 'ref', ref: 'nope.t'});
    let res = runCapturingLog(dir);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`table 'sample' column 'x' references 'nope.t', but no installed dependency ` +
      `declares schema 'nope'`);
    expect(res.output).toContain(`install the package that declares schema 'nope' as a dependency`);
    expect(fs.existsSync(dbPath(dir))).toBe(false);

    mutateManifest(dir, (m) => m.tables.sample.columns.x = {type: 'ref', ref: 'Core.nope'});
    res = runCapturingLog(dir);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`schema 'Core' declares no table 'nope'`);

    mutateManifest(dir, (m) => m.tables.sample.columns.x = {type: 'ref'});
    res = runCapturingLog(dir);
    expect(res.result).toBe(false);
    expect(res.output).toContain(`table 'sample' column 'x' is a ref column without a 'ref' target`);
  });

  it('qualified refs: a missing sealed Core declaration is an error, not a crash', () => {
    // the fixture's `owner: {type: user}` is the first Core ref the generator meets
    const dir = makePackage();
    const exists = fs.existsSync;
    const spy = vi.spyOn(fs, 'existsSync').mockImplementation((p) => !String(p).endsWith('Core.json') && exists(p));
    try {
      const {result, output} = runCapturingLog(dir);
      expect(result).toBe(false);
      expect(output).toMatch(
        /column 'owner' references 'Core\.users', but the sealed Core declaration is missing: .*Core\.json/);
      expect(fs.existsSync(dbPath(dir))).toBe(false);
    } finally {
      spy.mockRestore();
    }
  });

  it('relation usage compiles under strict tsc', {timeout: 180_000}, () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => addRelation(m));
    expect(generateDomainClients(dir)).toBe(true);
    const res = runTsc(dir, 'usage-relations.ts');
    expect(res.output.trim(), res.output).toBe('');
    expect(res.status).toBe(0);
  });

  it('importing the generated module without a live grok is side-effect free', () => {
    const dir = makePackage();
    expect(generateDomainClients(dir)).toBe(true);
    const js = ts.transpileModule(fs.readFileSync(dbPath(dir), 'utf8'), {compilerOptions: {
      module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2020,
    }}).outputText;
    const jsPath = path.join(dir, 'db.cjs');
    fs.writeFileSync(jsPath, js);
    // Stub 'datagrok-api/grok' with an empty module: importing must succeed (lazy
    // getters), touching a table getter must throw (it dereferences grok.dapi).
    const nm = path.join(dir, 'node_modules', 'datagrok-api');
    fs.mkdirSync(nm, {recursive: true});
    fs.writeFileSync(path.join(nm, 'package.json'), '{"name": "datagrok-api", "version": "0.0.0"}');
    fs.writeFileSync(path.join(nm, 'grok.js'), 'module.exports = {};');
    fs.writeFileSync(path.join(nm, 'dg.js'), 'module.exports = {};');
    const requireFrom = createRequire(jsPath);
    const mod = requireFrom(jsPath);
    expect(Object.keys(mod)).toContain('testdbDb');
    expect(() => mod.testdbDb.samples).toThrow();
  });

  it('is idempotent: rerunning produces byte-identical output', () => {
    const dir = makePackage();
    expect(generateDomainClients(dir)).toBe(true);
    const first = fs.readFileSync(dbPath(dir));
    expect(generateDomainClients(dir)).toBe(true);
    const second = fs.readFileSync(dbPath(dir));
    expect(second.equals(first)).toBe(true);
  });

  it('leaves packages without databases/ untouched', () => {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-api-nodb-'));
    fs.writeFileSync(path.join(dir, 'package.json'), '{"name": "@datagrok/empty", "version": "1.0.0"}');
    expect(generateDomainClients(dir)).toBe(true);
    expect(fs.existsSync(path.join(dir, 'src'))).toBe(false);
  });

  it('ignores legacy databases/ dirs that have no schema.json', () => {
    const dir = makePackage();
    fs.rmSync(path.join(dir, 'databases', 'testdb', 'schema.json'));
    fs.writeFileSync(path.join(dir, 'databases', 'testdb', '0000_init.sql'), 'CREATE TABLE t (id int);');
    expect(generateDomainClients(dir)).toBe(true);
    expect(fs.existsSync(dbPath(dir))).toBe(false);
  });

  it('rejects a manifest that violates the JSON Schema', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => m.tables.sample.columns.count.type = 'text');
    expect(generateDomainClients(dir)).toBe(false);
    expect(fs.existsSync(dbPath(dir))).toBe(false);
  });

  it('rejects a table referencing an unknown property schema', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => m.tables.sample.schemas = ['missing']);
    expect(generateDomainClients(dir)).toBe(false);
    expect(fs.existsSync(dbPath(dir))).toBe(false);
  });

  it('rejects a column named like a generated system field', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => m.tables.sample.columns.id = {type: 'string'});
    const {result, output} = runCapturingLog(dir);
    expect(result).toBe(false);
    expect(output).toContain(`table 'sample' column 'id' collides with a generated system column`);
    expect(fs.existsSync(dbPath(dir))).toBe(false);
  });

  it('rejects a relational column duplicated by a property-schema property', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => m.tables.sample.columns.note = {type: 'string'});
    const {result, output} = runCapturingLog(dir);
    expect(result).toBe(false);
    expect(output).toContain(`table 'sample' declares duplicate column 'note'`);
    expect(fs.existsSync(dbPath(dir))).toBe(false);
  });

  it('rejects two manifests in one package sharing a table name', () => {
    const dir = makePackage();
    const dbDir = path.join(dir, 'databases');
    fs.cpSync(path.join(dbDir, 'testdb'), path.join(dbDir, 'zzdb'), {recursive: true});
    const manifestPath = path.join(dbDir, 'zzdb', 'schema.json');
    const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
    manifest.name = 'zzdb';
    fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
    const {result, output} = runCapturingLog(dir);
    expect(result).toBe(false);
    expect(output).toContain(`table 'sample' emits interface 'SampleRow' already generated`);
    expect(fs.existsSync(dbPath(dir))).toBe(false);
  });

  it('generated code and well-typed usage compile under strict tsc', {timeout: 180_000}, () => {
    const dir = makePackage();
    expect(generateDomainClients(dir)).toBe(true);
    const res = runTsc(dir, 'usage-good.ts');
    expect(res.output.trim(), res.output).toBe('');
    expect(res.status).toBe(0);
  });

  it('wrong-typed usage fails tsc: values, choices, expand, columns, groupBy, datetime, tx', {timeout: 180_000}, () => {
    const dir = makePackage();
    expect(generateDomainClients(dir)).toBe(true);
    const res = runTsc(dir, 'usage-bad.ts');
    expect(res.status).not.toBe(0);
    // exactly one error per intended negative — no accidental extra breakage
    expect(res.output.match(/error TS/g)).toHaveLength(15);
    // wrong-typed payload values: `count: string` and `active: string`
    expect(res.output).toContain('usage-bad.ts(5,');
    expect(res.output).toContain('usage-bad.ts(7,');
    // the TInsert generic gates required columns: an empty insert payload must be rejected
    expect(res.output).toContain('usage-bad.ts(8,');
    expect(res.output).toMatch(/error TS2345: Argument of type '\{\}' is not assignable/);
    // wrong choices literal
    expect(res.output).toMatch(/'"bogus"' is not assignable/);
    // wrong expand key / column / groupBy entry
    expect(res.output).toMatch(/'"details:nope"' is not assignable/);
    expect(res.output).toContain('usage-bad.ts(11,');
    expect(res.output).toContain('usage-bad.ts(12,');
    // a datetime Row field rejects a bare string (dayjs on reads)
    expect(res.output).toMatch(/Type 'string' is not assignable to type 'Dayjs/);
    // transaction values must pair with their table
    expect(res.output).toContain('usage-bad.ts(16,');
    // builder negatives: unknown where/orderBy/select/expand columns...
    expect(res.output).toContain('usage-bad.ts(20,');
    expect(res.output).toContain('usage-bad.ts(21,');
    expect(res.output).toContain('usage-bad.ts(22,');
    expect(res.output).toContain('usage-bad.ts(23,');
    // ...non-selected field access after select()...
    expect(res.output).toMatch(/Property 'count' does not exist on type 'Pick<SampleRow/);
    // ...and a field absent from the mapped-tuple insert result
    expect(res.output).toMatch(/Property 'deleted' does not exist on type 'DomainInsertResult'/);
  });
});

describe('generateDomainClients --ui', () => {
  it('emits db-ui.ts only when asked, and never touches db.ts', () => {
    const dir = makePackage();
    expect(generateDomainClients(dir)).toBe(true);
    expect(fs.existsSync(dbUiPath(dir))).toBe(false);
    const dbOnly = fs.readFileSync(dbPath(dir));

    expect(generateDomainClients(dir, {ui: true})).toBe(true);
    expect(fs.existsSync(dbUiPath(dir))).toBe(true);
    // data-only consumers are unaffected: same db.ts, and it imports nothing UI
    expect(fs.readFileSync(dbPath(dir)).equals(dbOnly)).toBe(true);
    expect(fs.readFileSync(dbPath(dir), 'utf8')).not.toContain('domain-ui');
  });

  it('emits per-table option types and a <table>Ui wrapper over domain-ui', () => {
    const dir = makePackage();
    expect(generateDomainClients(dir, {ui: true})).toBe(true);
    const code = fs.readFileSync(dbUiPath(dir), 'utf8');

    expect(code).toContain('auto-generated by the grok api command (--ui)');
    expect(code).toContain(`from '@datagrok-libraries/domain-ui';`);
    // the typed names it reuses come from db.ts — nothing is re-derived here
    expect(code).toMatch(/import \{testdbDb, SampleColumn, SampleExpand, SampleInsert, SampleRow,/);
    expect(code).toMatch(/from '\.\/db';/);

    // per-table option types narrowed to this table's columns / insert payload
    expect(code).toContain(
      'export type SampleQuerySpec = DG.DomainQuerySpec<SampleColumn, keyof SampleExpand & string>;');
    expect(code).toContain(`  columns?: SampleColumn[];`);
    expect(code).toContain(
      `export interface SampleGridOptions extends Omit<DomainGridOptions, 'query' | 'defaults'> {`);
    expect(code).toContain('  defaults?: Partial<SampleInsert>;');
    expect(code).toContain(`export interface SampleListOptions extends Omit<EntityListOptions, 'query'> {`);
    expect(code).toContain(`export interface SampleAppViewOptions extends Omit<DomainAppViewOptions, 'query'> {`);
    expect(code).toContain(`export interface SampleFormOptions extends Omit<DomainFormOptions, 'values'> {`);
    expect(code).toContain('  values?: Partial<SampleInsert>;');

    // the wrapper itself: typed client + the typed handle every factory delegates to
    expect(code).toContain('export class SampleUi {');
    expect(code).toContain(`  readonly address: string = 'testdb.sample';`);
    expect(code).toContain(
      '  get client(): DG.DomainTableClient<SampleRow, SampleInsert, SampleColumn, SampleExpand> {');
    expect(code).toContain('    return testdbDb.samples;');
    expect(code).toContain(
      '  table(): Promise<DomainTable<SampleRow, SampleInsert, SampleColumn, SampleExpand>> {');
    expect(code).toContain(
      '    return domains.table<SampleRow, SampleInsert, SampleColumn, SampleExpand>(this.client);');
    // every widget member goes through that handle — no class constructor is named here
    expect(code).toContain('  async form(options?: SampleFormOptions): Promise<DomainForm> {');
    expect(code).toContain('    return (await this.table()).form(options);');
    expect(code).toContain(
      '  async formDialog(options?: SampleFormOptions & DomainDialogOptions): Promise<boolean> {');
    expect(code).toContain('  async listView(options?: SampleAppViewOptions): Promise<DomainAppView> {');
    expect(code).toContain('  app(options?: SampleAppViewOptions): Promise<DomainAppView> {');
    expect(code).toContain('  async grid(options?: SampleGridOptions): Promise<DomainGrid> {');
    // WO-9R reverted EntityListWidget.create's null return — the nit is gone
    expect(code).toContain('  async list(options?: SampleListOptions): Promise<EntityListWidget> {');
    expect(code).not.toContain('EntityListWidget | null');
    expect(code).not.toMatch(/new Domain(AppView|EntityAppView)\(/);
    expect(code).toContain('  row(values: Partial<SampleRow> | null): DG.DomainRow {');
    expect(code).toContain(
      `    return new DG.DomainQuery({...params, schema: 'testdb', table: 'sample'});`);
    expect(code).toContain('export const sampleUi = new SampleUi();');
    // the snake_case table gets the camelCase client getter and its own wrapper
    expect(code).toContain('    return testdbDb.sampleEvents;');
    expect(code).toContain('export const sampleEventUi = new SampleEventUi();');

    expect(code).not.toMatch(/[^\r]\n/);   // CRLF, per the repo code style
  });

  it('the opt-in persists: a plain rerun keeps db-ui.ts up to date, byte-identically', () => {
    const dir = makePackage();
    expect(generateDomainClients(dir, {ui: true})).toBe(true);
    const first = fs.readFileSync(dbUiPath(dir));
    // no flag this time — the existing file is what opts the package in
    expect(generateDomainClients(dir)).toBe(true);
    expect(fs.readFileSync(dbUiPath(dir)).equals(first)).toBe(true);

    mutateManifest(dir, (m) => m.tables.sample.columns.extra_note = {type: 'string'});
    expect(generateDomainClients(dir)).toBe(true);
    expect(fs.readFileSync(dbUiPath(dir), 'utf8')).toContain('SampleUi');
    expect(fs.readFileSync(dbUiPath(dir), 'utf8')).toContain('export function testdbUiDb(): Promise<TestdbUiDb>');
    expect(fs.readFileSync(dbPath(dir), 'utf8')).toContain('extra_note?: string;');
  });

  it('typed UI usage compiles under strict tsc', {timeout: 180_000}, () => {
    const dir = makePackage();
    expect(generateDomainClients(dir, {ui: true})).toBe(true);
    const res = runTsc(dir, 'usage-ui-good.ts');
    expect(res.output.trim(), res.output).toBe('');
    expect(res.status).toBe(0);
  });

  it('relations: the handle takes the update payload as its fifth generic', () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => addRelation(m));
    expect(generateDomainClients(dir, {ui: true})).toBe(true);
    const code = fs.readFileSync(dbUiPath(dir), 'utf8');

    expect(code).toMatch(/SampleRow, SampleUpdate,/);
    expect(code).toContain('  get client(): DG.DomainTableClient<SampleRow, SampleInsert, ' +
      'SampleColumn, SampleExpand, SampleUpdate> {');
    expect(code).toContain('  table(): Promise<DomainTable<SampleRow, SampleInsert, ' +
      'SampleColumn, SampleExpand, SampleUpdate>> {');
    expect(code).toContain('    return domains.table<SampleRow, SampleInsert, ' +
      'SampleColumn, SampleExpand, SampleUpdate>(this.client);');
    // the schema handle's per-table property carries it too
    expect(code).toContain('  readonly samples: DomainTable<SampleRow, SampleInsert, ' +
      'SampleColumn, SampleExpand, SampleUpdate>;');
    // a relation-less table keeps four — its update payload is the client's Partial<Row>
    expect(code).toContain(
      '  get client(): DG.DomainTableClient<LabelRow, LabelInsert, LabelColumn, LabelExpand> {');
    expect(code).not.toContain('LabelUpdate');
  });

  it('relation UI usage compiles under strict tsc', {timeout: 180_000}, () => {
    const dir = makePackage();
    mutateManifest(dir, (m) => addRelation(m));
    expect(generateDomainClients(dir, {ui: true})).toBe(true);
    const res = runTsc(dir, 'usage-ui-relations.ts');
    expect(res.output.trim(), res.output).toBe('');
    expect(res.status).toBe(0);
  });

  it('wrong-typed UI usage fails tsc: columns, expand, filter, insert defaults, choices',
    {timeout: 180_000}, () => {
      const dir = makePackage();
      expect(generateDomainClients(dir, {ui: true})).toBe(true);
      const res = runTsc(dir, 'usage-ui-bad.ts');
      expect(res.status).not.toBe(0);
      // exactly one error per intended negative — no accidental extra breakage
      expect(res.output.match(/error TS/g)).toHaveLength(14);
      // a wrong column name in a query spec, an expand key, and a condition
      expect(res.output).toMatch(/usage-ui-bad\.ts\(5,.*'"nope"' is not assignable to type 'SampleColumn'/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(8,.*'"details:nope"' is not assignable/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(9,/);
      expect(res.output).toMatch(/Type '"nope"' is not assignable to type 'DomainColumnRef<SampleColumn>/);
      // insert defaults: unknown column, wrong type, a column of ANOTHER table
      expect(res.output).toMatch(/usage-ui-bad\.ts\(6,.*'nope' does not exist in type 'Partial<SampleInsert>'/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(7,/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(11,.*'kind' does not exist in type 'Partial<SampleInsert>'/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(12,/);
      // rowFrom values and DomainQuery columns are this table's too
      expect(res.output).toMatch(/usage-ui-bad\.ts\(13,.*'nope' does not exist in type 'Partial<SampleRow>'/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(14,.*'"nope"' is not assignable to type 'SampleColumn'/);
      // choices stay the generated literal union
      expect(res.output).toMatch(/usage-ui-bad\.ts\(16,.*'"bogus"' is not assignable to type 'SampleStatus/);
      // Â§1.7: the form's seed values are this table's insert payload
      expect(res.output).toMatch(/usage-ui-bad\.ts\(18,.*'nmae' does not exist in type 'Partial<SampleInsert>'/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(19,/);
      // the schema handle's per-table properties carry the same typing
      expect(res.output).toMatch(/usage-ui-bad\.ts\(22,.*'nmae' does not exist in type 'Partial<SampleInsert>'/);
      expect(res.output).toMatch(/usage-ui-bad\.ts\(23,/);
    });
});
