import {describe, it, expect, vi} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {fileURLToPath} from 'url';
import {Project} from 'ts-morph';
import {buildManifests, readEntityClasses, snakeCase, EntityClass} from '../utils/entity-classes';
import {buildSnapshot, schemaDdl} from '../utils/domain-snapshot';
import {schema} from '../commands/schema';

const testsDir = path.dirname(fileURLToPath(import.meta.url));
const fixtureDir = path.join(testsDir, 'fixtures', 'entity-package');

const expectedManifest = {
  $schema: 'https://datagrok.ai/schemas/domain-schema.schema.json',
  name: 'lab',
  version: '1.2.3',
  tables: {
    sample: {
      securityMode: 'row', idempotency: true, businessKey: ['name'],
      columns: {
        name: {type: 'string', required: true, unique: true, isName: true},
        count: {type: 'float', min: 0},
        active: {type: 'bool'},
        measured_on: {type: 'datetime'},
        tags: {type: 'string_list'},
        owner: {type: 'user'},
        status: {type: 'string', choices: ['new', 'done'], default: 'new'},
      },
    },
    plate: {
      friendlyName: 'Plates', singularName: 'Plate',
      columns: {
        barcode: {type: 'string', required: true, unique: true},
        rows: {type: 'int', min: 1, max: 32},
      },
    },
    plate_well: {
      securityMode: 'master', delegate: 'plate_id', audit: false, businessKey: ['plate_id', 'well_position'],
      columns: {
        plate_id: {type: 'ref', ref: 'plate', onDelete: 'cascade', required: true},
        sample_id: {type: 'ref', ref: 'sample', onDelete: 'setnull'},
        well_position: {type: 'string', required: true, friendlyName: 'Position'},
      },
    },
    reader: {
      description: 'Plate readers',
      columns: {
        name: {type: 'string', required: true},
        config: {type: 'file'},
      },
    },
  },
};

const assayClass = `
@grok.decorators.entity({schema: 'lab'})
export class Assay {
  @grok.decorators.column({ref: () => Sample, onDelete: 'cascade', required: true})
  sample_id!: string;

  @grok.decorators.column({ref: () => Reader})
  reader_id?: string;

  @grok.decorators.column()
  result!: number;
}
`;

function makePackage(): string {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-schema-generate-'));
  fs.cpSync(fixtureDir, dir, {recursive: true});
  return dir;
}

/** The entity classes of one in-memory source file. */
function classesOf(code: string): EntityClass[] {
  const project = new Project({useInMemoryFileSystem: true, skipFileDependencyResolution: true});
  project.createSourceFile('entities.ts', code);
  return readEntityClasses(project);
}

const manifestOf = (code: string) => buildManifests(classesOf(code), '0.0.1');

/** Runs the command with console output captured; returns its result, exit code, and output. */
async function run(dir: string, verb: string, options: {[key: string]: any} = {}):
Promise<{result: boolean, exitCode: number | undefined, output: string}> {
  const spy = vi.spyOn(console, 'log').mockImplementation(() => {});
  process.exitCode = undefined;
  try {
    const result = await schema({_: ['schema', verb], ...options}, dir);
    return {result, exitCode: process.exitCode as number | undefined,
      output: spy.mock.calls.map((args) => args.join(' ')).join('\n')};
  } finally {
    process.exitCode = undefined;
    spy.mockRestore();
  }
}

const manifestPath = (dir: string) => path.join(dir, 'databases', 'lab', 'schema.json');

describe('grok schema generate', () => {
  it('writes the manifest of the four fixture classes', async () => {
    const dir = makePackage();
    const res = await run(dir, 'generate');
    expect(res.result).toBe(true);
    expect(res.exitCode).toBeUndefined();
    expect(res.output).toContain(`Wrote ${path.join('databases', 'lab', 'schema.json')} (4 table(s) from 4 class(es))`);
    const text = fs.readFileSync(manifestPath(dir), 'utf8');
    expect(text.startsWith('{\n  "$schema"')).toBe(true);
    expect(text.endsWith('}\n')).toBe(true);
    const manifest = JSON.parse(text);
    expect(manifest).toEqual(expectedManifest);
    // declaration order is column order; class order is table order
    expect(Object.keys(manifest.tables)).toEqual(['sample', 'plate', 'plate_well', 'reader']);
    expect(Object.keys(manifest.tables.sample.columns))
      .toEqual(['name', 'count', 'active', 'measured_on', 'tags', 'owner', 'status']);
    expect(Object.keys(manifest.tables.plate_well.columns.plate_id)).toEqual(['type', 'ref', 'onDelete', 'required']);
    // a rerun is byte-identical
    await run(dir, 'generate');
    expect(fs.readFileSync(manifestPath(dir), 'utf8')).toBe(text);
  });

  it('keeps the sections classes cannot declare, and groups classes by schema', async () => {
    const dir = makePackage();
    await run(dir, 'generate');
    const manifest = JSON.parse(fs.readFileSync(manifestPath(dir), 'utf8'));
    manifest.extensible = {tables: true};
    manifest.propertySchemas = {extra: {note: {type: 'string'}}};
    fs.writeFileSync(manifestPath(dir), JSON.stringify(manifest));
    fs.writeFileSync(path.join(dir, 'src', 'qc.ts'), `import * as grok from 'datagrok-api/grok';
@grok.decorators.entity({schema: 'qc', schemas: ['checks']})
export class Run { @grok.decorators.column() started_on!: Date; }
`);
    fs.mkdirSync(path.join(dir, 'databases', 'qc'), {recursive: true});
    fs.writeFileSync(path.join(dir, 'databases', 'qc', 'schema.json'),
      JSON.stringify({name: 'qc', version: '0', tables: {x: {columns: {}}}, propertySchemas: {checks: {ok: {type: 'bool'}}}}));
    const res = await run(dir, 'generate');
    expect(res.exitCode).toBeUndefined();
    const lab = JSON.parse(fs.readFileSync(manifestPath(dir), 'utf8'));
    expect(Object.keys(lab)).toEqual(['$schema', 'name', 'version', 'extensible', 'tables', 'propertySchemas']);
    expect(lab.extensible).toEqual({tables: true});
    expect(lab.propertySchemas).toEqual({extra: {note: {type: 'string'}}});
    expect(lab.tables).toEqual(expectedManifest.tables);
    const qc = JSON.parse(fs.readFileSync(path.join(dir, 'databases', 'qc', 'schema.json'), 'utf8'));
    expect(qc.tables).toEqual({run: {schemas: ['checks'], columns: {started_on: {type: 'datetime'}}}});
    expect(qc.propertySchemas).toEqual({checks: {ok: {type: 'bool'}}});
    expect(qc.version).toBe('1.2.3');
  });

  it('infers column types from the property types, and names tables in snake_case', () => {
    const m = manifestOf(`
@grok.decorators.entity({schema: 's'})
class HTTPReaderConfig {
  @grok.decorators.column() a!: string;
  @grok.decorators.column() b!: number;
  @grok.decorators.column() c!: boolean;
  @grok.decorators.column() d!: Date;
  @grok.decorators.column() e!: string[];
  @grok.decorators.column() f!: Array<string>;
  @grok.decorators.column() g = 0;
  @grok.decorators.column() h?: string;
  @grok.decorators.column({type: 'int'}) i!: number;
  @grok.decorators.column({type: 'group'}) j!: string;
  @grok.decorators.column({ref: 'other', type: 'user'}) k!: string;
  @grok.decorators.column({min: -1.5}) l!: number;
}`);
    expect(Object.keys(m.s.tables)).toEqual(['http_reader_config']);
    const columns = m.s.tables.http_reader_config.columns;
    expect(Object.keys(columns).map((n) => columns[n].type))
      .toEqual(['string', 'float', 'bool', 'datetime', 'string_list', 'string_list', 'float', 'string', 'int', 'group', 'user', 'float']);
    // a user/group column drops a stray ref (the manifest allows ref on type ref only)
    expect(columns.k).toEqual({type: 'user'});
    expect(columns.l).toEqual({type: 'float', min: -1.5});
    expect(snakeCase('PlateWell')).toBe('plate_well');
    expect(snakeCase('Sample')).toBe('sample');
    expect(snakeCase('HTTPReader')).toBe('http_reader');
  });

  it('accepts decorators.entity when the namespace is imported by name, and skips other classes', () => {
    const classes = classesOf(`
import {decorators} from 'datagrok-api/grok';
@decorators.entity({schema: 's', table: 'things'})
class Thing { @decorators.column() name!: string; }
class Helper { name!: string; }
@somethingElse()
class Other { @decorators.column() name!: string; }`);
    expect(classes.map((c) => [c.className, c.table])).toEqual([['Thing', 'things']]);
  });

  it('rejects what the manifest cannot express, naming the class and option', () => {
    const entity = (body: string, options: string = `{schema: 's'}`) =>
      `@grok.decorators.entity(${options})\nclass Thing {\n${body}\n}`;
    expect(() => manifestOf(entity(`@grok.decorators.column() blob!: Uint8Array;`)))
      .toThrow(`Thing.blob: cannot infer a column type from 'Uint8Array'; declare {type: ...}`);
    expect(() => manifestOf(entity(`@grok.decorators.column({ref: () => Helper}) h!: string;`) + `\nclass Helper {}`))
      .toThrow(`Thing.h: ref target 'Helper' is not an @entity class`);
    expect(() => manifestOf(entity(`@grok.decorators.column({ref: () => Far}) h!: string;`) +
      `\n@grok.decorators.entity({schema: 'other'}) class Far {}`))
      .toThrow(`Thing.h: ref target 'Far' belongs to schema 'other', not 's'`);
    expect(() => manifestOf(entity(`@grok.decorators.column({ref: () => this.x}) h!: string;`)))
      .toThrow(`Thing.h: @column option 'ref' must be a table name or '() => EntityClass'`);
    expect(() => manifestOf(entity(`@grok.decorators.column({type: 'ref'}) h!: string;`)))
      .toThrow(`Thing.h: a 'ref' column needs its target ('ref')`);
    expect(() => manifestOf(entity(`@grok.decorators.column() n!: string;`) +
      `\n@grok.decorators.entity({schema: 's', table: 'thing'}) class Thing2 {}`))
      .toThrow(`Thing2 and Thing both map to table 'thing'`);
    expect(() => manifestOf(`const keys = ['a'];\n` + entity(`@grok.decorators.column() a!: string;`, `{schema: 's', businessKey: keys}`)))
      .toThrow(`Thing: @entity option 'businessKey' must be a string, number, boolean, array or object literal (found 'keys')`);
    expect(() => manifestOf(entity(`@grok.decorators.column({friendlyName: 'a' + 'b'}) a!: string;`)))
      .toThrow(`Thing.a: @column option 'friendlyName' must be a string, number, boolean, array or object literal (found ''a' + 'b'')`);
    expect(() => manifestOf(entity(`@grok.decorators.column() a!: string;`, `{}`))).toThrow(`Thing: @entity needs a 'schema'`);
    expect(() => manifestOf(entity(`@grok.decorators.column() a!: string;`, `{schema: 's', colour: 'red'}`)))
      .toThrow(`Thing: unknown @entity option 'colour'`);
    expect(() => manifestOf(entity(`@grok.decorators.column({colour: 'red'}) a!: string;`)))
      .toThrow(`Thing.a: unknown @column option 'colour'`);
    expect(() => manifestOf(entity(`@grok.decorators.column() a!: string;`, `{schema: 's', delegate: 'b'}`)))
      .toThrow(`Thing: delegate 'b' is not a @column of the class`);
    expect(() => manifestOf(entity(`@grok.decorators.column() a!: string;`, `{schema: 's', businessKey: ['a', 'b']}`)))
      .toThrow(`Thing: businessKey column 'b' is not a @column of the class`);
    expect(() => manifestOf(entity(`@grok.decorators.column({isName: true}) a!: number;`)))
      .toThrow(`Thing: isName column 'a' must be a string`);
    expect(() => manifestOf(entity(`@grok.decorators.column({name: 'x'}) a!: string;\n@grok.decorators.column() x!: string;`)))
      .toThrow(`Thing: two properties map to column 'x'`);
  });
});

describe('grok schema ddl', () => {
  it('prints CREATE TABLE for every table, each after the tables it refs', async () => {
    const dir = makePackage();
    await run(dir, 'generate');
    const res = await run(dir, 'ddl');
    expect(res.result).toBe(true);
    expect(res.exitCode).toBeUndefined();
    const manifest = JSON.parse(fs.readFileSync(manifestPath(dir), 'utf8'));
    const hash = buildSnapshot(manifest).hash!;
    expect(res.output.startsWith(`-- lab (${hash.substring(0, 12)})\n\nCREATE TABLE IF NOT EXISTS lab.plate (\n`)).toBe(true);
    const tables = [...res.output.matchAll(/CREATE TABLE IF NOT EXISTS lab\.(\w+) \(/g)].map((m) => m[1]);
    expect(tables).toEqual(['plate', 'sample', 'plate_well', 'reader']);
    expect(res.output).toContain('CREATE TABLE IF NOT EXISTS lab.plate_well (\n' +
      '  id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),\n' +
      '  is_deleted bool NOT NULL DEFAULT false,\n' +
      '  version int NOT NULL DEFAULT 1,\n' +
      '  created_on timestamp without time zone NOT NULL DEFAULT now(),\n' +
      '  updated_on timestamp without time zone NOT NULL DEFAULT now(),\n' +
      '  author_id uuid,\n' +
      "  data jsonb NOT NULL DEFAULT '{}',\n" +
      '  plate_id uuid NOT NULL,\n' +
      '  sample_id uuid,\n' +
      '  well_position text NOT NULL\n' +
      ');\n' +
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_plate_well_business_key ON lab.plate_well (plate_id, well_position) ' +
      'WHERE NOT is_deleted;\n');
    expect(res.output).toContain('REFERENCES lab.plate(id);');
    expect(res.output).toContain('REFERENCES lab.sample(id);');
    expect(res.output).toContain('REFERENCES public.users(id);');
    expect(res.output).toContain('  idempotency_key uuid,\n');
    expect(res.output).toContain('CREATE UNIQUE INDEX IF NOT EXISTS ux_sample_idempotency_key ON lab.sample (idempotency_key) ' +
      'WHERE idempotency_key IS NOT NULL;\n');

    // --out writes the same text instead of printing it
    const out = await run(dir, 'ddl', {out: 'lab.sql'});
    expect(out.output).toBe('Wrote lab.sql');
    expect(fs.readFileSync(path.join(dir, 'lab.sql'), 'utf8').trimEnd()).toBe(res.output);
  });

  it('refuses a reference cycle', () => {
    const s = buildSnapshot({name: 'x', version: '1', tables: {
      a: {columns: {b_id: {type: 'ref', ref: 'b'}, self_id: {type: 'ref', ref: 'a'}}},
      b: {columns: {a_id: {type: 'ref', ref: 'a'}}},
    }});
    expect(() => schemaDdl(s)).toThrow('reference cycle: a -> b -> a');
  });
});

describe('generate → seal → migrate', () => {
  it('scaffolds the migration for a fifth class', async () => {
    const dir = makePackage();
    await run(dir, 'generate');
    await run(dir, 'seal');
    expect((await run(dir, 'check')).exitCode).toBeUndefined();
    fs.appendFileSync(path.join(dir, 'src', 'entities.ts'), assayClass);
    let res = await run(dir, 'generate');
    expect(res.output).toContain('(5 table(s) from 5 class(es))');
    const manifest = JSON.parse(fs.readFileSync(manifestPath(dir), 'utf8'));
    expect(manifest.tables.assay).toEqual({columns: {
      sample_id: {type: 'ref', ref: 'sample', onDelete: 'cascade', required: true},
      reader_id: {type: 'ref', ref: 'reader'},
      result: {type: 'float'},
    }});
    expect((await run(dir, 'check')).exitCode).toBe(1);

    res = await run(dir, 'migrate', {name: 'add_assay'});
    expect(res.exitCode).toBeUndefined();
    expect(res.output).toContain('add-table assay: created [physical]');
    expect(res.output).toContain('No change needs a migration script');
    const up = fs.readFileSync(path.join(dir, 'databases', 'lab', '0001_add_assay.sql'), 'utf8');
    expect(up).toContain('-- lab: add_assay\n-- ems-migration: from=');
    expect(up).toContain('-- add-table assay: created [physical]\nCREATE TABLE IF NOT EXISTS lab.assay (\n');
    expect(up).toContain('  sample_id uuid NOT NULL,\n  reader_id uuid,\n  result float8\n);\n');
    expect(up).toContain('ALTER TABLE lab.assay ADD CONSTRAINT fk_assay_sample_id FOREIGN KEY (sample_id) REFERENCES lab.sample(id);');
    expect(up).toContain('ALTER TABLE lab.assay ADD CONSTRAINT fk_assay_reader_id FOREIGN KEY (reader_id) REFERENCES lab.reader(id);');
    expect(fs.readFileSync(path.join(dir, 'databases', 'lab', 'down', '0001_add_assay.sql'), 'utf8'))
      .toContain('-- data loss\nDROP TABLE IF EXISTS lab.assay;\n');
    expect((await run(dir, 'check')).exitCode).toBeUndefined();
  });
});

describe('grok schema generate — the server\'s quick-wins keys', () => {
  it('carries json and immutable columns, constraints and grants into the manifest', () => {
    const m = manifestOf(`
import * as grok from 'datagrok-api/grok';

@grok.decorators.entity({schema: 'lab', securityMode: 'row',
  constraints: {weight_positive: {check: 'weight IS NULL OR weight > 0'}},
  grants: {'#all-users': ['view']}})
export class Term {
  @grok.decorators.column({required: true, isName: true, immutable: true}) code!: string;
  @grok.decorators.column() weight?: number;
  @grok.decorators.column({type: 'json'}) meta?: object;
  @grok.decorators.column() props?: object;
}
`);
    const t = m.lab.tables.term;
    expect(t.constraints).toEqual({weight_positive: {check: 'weight IS NULL OR weight > 0'}});
    expect(t.grants).toEqual({'#all-users': ['view']});
    expect(t.columns.code).toEqual({type: 'string', required: true, isName: true, immutable: true});
    expect(t.columns.meta).toEqual({type: 'json'});
    expect(t.columns.props).toEqual({type: 'json'});
    expect(buildSnapshot(m.lab).tables.term.constraints).toEqual({weight_positive: 'weight is null or weight > 0'});
  });
});
