import {describe, it, expect, vi, afterEach} from 'vitest';
import fs from 'fs';
import os from 'os';
import path from 'path';
import {createHash} from 'crypto';
import {fileURLToPath} from 'url';
import {buildSnapshot, canonical, diff, emptySnapshot, hashOf, scaffold, serialize, snapshotBody,
  Change} from '../utils/domain-snapshot';
import {schema} from '../commands/schema';
import {createClient} from '../utils/server-client';

vi.mock('../utils/server-client', () => ({createClient: vi.fn()}));

const testsDir = path.dirname(fileURLToPath(import.meta.url));
const fixtureDir = path.join(testsDir, 'fixtures', 'domain-package');
const fixturePath = path.join(fixtureDir, 'databases', 'testdb', 'schema.json');

function fixture(mutate?: (m: any) => void): any {
  const m = JSON.parse(fs.readFileSync(fixturePath, 'utf8'));
  if (mutate)
    mutate(m);
  return m;
}

/** Adds a `label` table, a `sample_label` junction, and the `labels` relation on `sample`. */
function addRelation(m: any, relation: {[key: string]: any} = {via: 'sample_label', target: 'label'}): void {
  m.tables.label = {businessKey: ['name'], columns: {name: {type: 'string', required: true}}};
  m.tables.sample_label = {
    securityMode: 'master', delegate: 'sample_id', businessKey: ['sample_id', 'label_id'],
    columns: {
      sample_id: {type: 'ref', ref: 'sample', onDelete: 'cascade', required: true},
      label_id: {type: 'ref', ref: 'label', onDelete: 'cascade', required: true},
    }};
  m.tables.sample.relations = {labels: relation};
}

/** The changes between the fixture and a mutation of it. */
function changesOf(mutate: (m: any) => void): Change[] {
  return diff(buildSnapshot(fixture()), buildSnapshot(fixture(mutate)));
}

const flags = (c: Change) => [c.kind, c.physical, c.auto];

describe('buildSnapshot', () => {
  it('normalizes the fixture manifest: sorted tables, defaults omitted, arrays in key order', () => {
    const s = buildSnapshot(fixture());
    expect(s.format).toBe(1);
    expect(s.schema).toBe('testdb');
    expect(s.version).toBe('1.0.0');
    expect(s.extensibleTables).toBeUndefined();
    expect(Object.keys(s.tables)).toEqual(['sample', 'sample_event']);

    const sample = s.tables.sample;
    expect(sample.securityMode).toBe('row');
    expect(sample.promotion).toBeUndefined();
    expect(sample.businessKey).toEqual(['name']);
    expect(sample.idempotency).toBe(true);
    expect(sample.schemas).toEqual(['extra']);
    expect(sample.nameColumn).toBe('name');
    expect(sample.columns.map((c) => c.name))
      .toEqual(['name', 'count', 'ratio', 'active', 'measured_on', 'tags', 'owner', 'status']);
    expect(sample.columns[0]).toEqual({name: 'name', type: 'string', required: true, unique: true});
    expect(sample.columns[1]).toEqual({name: 'count', type: 'int', min: 0});
    expect(sample.columns[6]).toEqual({name: 'owner', type: 'user'});
    expect(sample.columns[7]).toEqual({name: 'status', type: 'string', choices: ['new', 'done']});

    const event = s.tables.sample_event;
    expect(event).toMatchObject({securityMode: 'master', delegate: 'sample_id', audit: false});
    expect(event.nameColumn).toBeUndefined();
    expect(event.columns[0]).toEqual(
      {name: 'sample_id', type: 'ref', ref: 'sample', onDelete: 'cascade', required: true});

    expect(s.propertySchemas).toEqual({extra: [
      {name: 'note', type: 'string'}, {name: 'score', type: 'float', required: true}]});
  });

  it('resolves nameColumn: isName wins, then a string column named name, else absent', () => {
    expect(buildSnapshot(fixture((m) => m.tables.sample.columns.kind = {type: 'string', isName: true}))
      .tables.sample.nameColumn).toBe('kind');
    expect(buildSnapshot(fixture((m) => m.propertySchemas.extra.note.isName = true))
      .tables.sample.nameColumn).toBe('note');
    expect(buildSnapshot(fixture((m) => m.tables.sample_event.columns.name = {type: 'int'}))
      .tables.sample_event.nameColumn).toBeUndefined();
    expect(buildSnapshot(fixture((m) => m.tables.sample_event.columns.name = {type: 'string'}))
      .tables.sample_event.nameColumn).toBe('name');
  });

  it('emits promotion only for row-mode tables, and only when not lazy', () => {
    expect(buildSnapshot(fixture((m) => m.tables.sample.promotion = 'eager')).tables.sample.promotion).toBe('eager');
    expect(buildSnapshot(fixture((m) => m.tables.sample.promotion = 'lazy')).tables.sample.promotion).toBeUndefined();
    expect(buildSnapshot(fixture((m) => m.tables.sample_event.promotion = 'eager'))
      .tables.sample_event.promotion).toBeUndefined();
  });

  it('resolves relation FK sides from the junction and emits relations as an array', () => {
    const s = buildSnapshot(fixture((m) => addRelation(m)));
    expect(s.tables.sample.relations).toEqual([
      {name: 'labels', via: 'sample_label', target: 'label', viaSelf: 'sample_id', viaTarget: 'label_id'}]);
    expect(buildSnapshot(fixture((m) => addRelation(m, {via: 'sample_label', target: 'label', allowCreate: false})))
      .tables.sample.relations![0].allowCreate).toBe(false);
    expect(buildSnapshot(fixture((m) => addRelation(m, {via: 'sample_label', target: 'label', allowCreate: true})))
      .tables.sample.relations![0].allowCreate).toBeUndefined();
    expect(() => buildSnapshot(fixture((m) => {
      addRelation(m);
      m.tables.sample_label.columns.other_label_id = {type: 'ref', ref: 'label'};
    }))).toThrow(/more than one ref column targeting 'label'/);
    expect(buildSnapshot(fixture((m) => {
      addRelation(m, {via: 'sample_label', target: 'label', viaTarget: 'other_label_id'});
      m.tables.sample_label.columns.other_label_id = {type: 'ref', ref: 'label'};
    })).tables.sample.relations![0].viaTarget).toBe('other_label_id');
  });

  it('normalizes filters and display names the way the server does', () => {
    const s = buildSnapshot(fixture((m) => {
      m.tables.sample.filters = [{column: 'status', label: '  Status '}, {column: 'count', type: 'histogram', bins: 10}];
      m.tables.sample.singularName = '  Sample ';
      m.tables.sample.pluralName = '   ';
      m.tables.sample.columns.count.default = 5;
      m.tables.sample.columns.active.default = true;
    }));
    expect(s.tables.sample.filters).toEqual([{column: 'status', label: 'Status'}, {column: 'count', type: 'histogram', bins: 10}]);
    expect(s.tables.sample.singularName).toBe('Sample');
    expect(s.tables.sample.pluralName).toBeUndefined();
    expect(s.tables.sample.columns[1].default).toBe('5');
    expect(s.tables.sample.columns[3].default).toBe('true');
  });
});

describe('hash', () => {
  const base = buildSnapshot(fixture()).hash!;

  it('is the SHA-256 of the canonical body and survives a serialize round trip', () => {
    const s = buildSnapshot(fixture());
    expect(base).toBe(createHash('sha256').update(canonical(snapshotBody(s))).digest('hex'));
    const parsed = JSON.parse(serialize(s));
    expect(Object.keys(parsed)).toEqual(['format', 'schema', 'version', 'hash', 'propertySchemas', 'tables']);
    expect(parsed.hash).toBe(base);
    expect(hashOf(parsed)).toBe(base);
    expect(serialize(s).endsWith('\n')).toBe(true);
  });

  it('ignores object key order, explicit defaults, version and schema name', () => {
    expect(buildSnapshot(fixture((m) => {
      m.tables.sample.columns.name = {unique: true, required: true, type: 'string'};
      m.tables.sample = {columns: m.tables.sample.columns, schemas: m.tables.sample.schemas,
        idempotency: true, businessKey: ['name'], securityMode: 'row'};
    })).hash).toBe(base);
    expect(buildSnapshot(fixture((m) => {
      m.version = '2.0.0';
      m.name = 'other';
      m.tables.sample.promotion = 'lazy';
      m.tables.sample.defaultRowVisibility = 'table';
      m.tables.sample.softDelete = true;
      m.tables.sample.audit = true;
      m.tables.sample.extensible = false;
      m.tables.sample.ginIndex = false;
      m.tables.sample.filters = [];
      m.tables.sample.columns.count.required = false;
      m.tables.sample.columns.count.unique = false;
      m.tables.sample.columns.count.isName = false;
      m.tables.sample.columns.count.choices = [];
      m.tables.sample_event.schemas = [];
      m.tables.sample_event.idempotency = false;
    })).hash).toBe(base);
  });

  it('changes with column order and with any declared field', () => {
    expect(buildSnapshot(fixture((m) => {
      const {name, ...rest} = m.tables.sample.columns;
      m.tables.sample.columns = {...rest, name};
    })).hash).not.toBe(base);
    expect(buildSnapshot(fixture((m) => m.tables.sample.columns.count.min = 1)).hash).not.toBe(base);
    expect(buildSnapshot(fixture((m) => m.extensible = {tables: true})).hash).not.toBe(base);
  });

  it('canonical: sorted keys by code unit, integral numbers as integers, no undefined', () => {
    expect(canonical({b: 1.0, a: [2.5, 'x', true, null], B: {}, _: [], d: undefined}))
      .toBe('{"B":{},"_":[],"a":[2.5,"x",true,null],"b":1}');
    expect(canonical('q"\\\n')).toBe('"q\\"\\\\\\n"');
  });
});

describe('diff', () => {
  it('is empty between equal snapshots and all add-table from an empty baseline', () => {
    const s = buildSnapshot(fixture());
    expect(diff(s, s)).toEqual([]);
    expect(diff(emptySnapshot('testdb'), s).map(flags))
      .toEqual([['add-table', true, true], ['add-table', true, true], ['property-schema', false, true]]);
  });

  it('classifies column additions and drops, with a rename hint for a same-typed pair', () => {
    expect(changesOf((m) => m.tables.sample.columns.note2 = {type: 'string'}).map(flags))
      .toEqual([['add-column', true, true]]);
    expect(changesOf((m) => m.tables.sample.columns.code = {type: 'string', required: true})[0])
      .toMatchObject({kind: 'add-column', table: 'sample', column: 'code', physical: true, auto: false,
        detail: 'added (string, required without a default)'});
    expect(changesOf((m) => m.tables.sample.columns.code = {type: 'string', required: true, default: 'x'})[0])
      .toMatchObject({kind: 'add-column', auto: true});
    expect(changesOf((m) => delete m.tables.sample.columns.count).map(flags)).toEqual([['drop-column', true, false]]);

    const renamed = changesOf((m) => {
      delete m.tables.sample.columns.count;
      m.tables.sample.columns.count2 = {type: 'int', min: 0};
    });
    expect(renamed.map((c) => c.kind)).toEqual(['drop-column', 'add-column']);
    for (const c of renamed)
      expect(c.detail).toContain('(rename? count -> count2)');
    const notRenamed = changesOf((m) => {
      delete m.tables.sample.columns.count;
      m.tables.sample.columns.label = {type: 'string'};
    });
    expect(notRenamed.every((c) => !c.detail.includes('rename?'))).toBe(true);
  });

  it('classifies type, required, unique, ref and onDelete changes', () => {
    expect(changesOf((m) => m.tables.sample.columns.count.type = 'float')[0])
      .toMatchObject({kind: 'change-type', column: 'count', detail: 'type int -> float', physical: true, auto: false});
    expect(changesOf((m) => m.tables.sample.columns.name.required = false).map(flags))
      .toEqual([['change-required', true, false]]);
    expect(changesOf((m) => m.tables.sample.columns.count.unique = true).map(flags))
      .toEqual([['change-unique', true, true]]);
    expect(changesOf((m) => m.tables.sample.columns.name.unique = false).map(flags))
      .toEqual([['change-unique', true, false]]);
    expect(changesOf((m) => m.tables.sample_event.columns.sample_id.ref = 'label')[0])
      .toMatchObject({kind: 'change-ref', detail: 'ref "sample" -> "label"', physical: true, auto: false});
    expect(changesOf((m) => m.tables.sample_event.columns.sample_id.onDelete = 'restrict')[0])
      .toMatchObject({kind: 'change-on-delete', physical: false, auto: false});
    expect(changesOf((m) => m.tables.sample.columns.count.min = 1)[0])
      .toMatchObject({kind: 'column-meta', detail: 'min 0 -> 1', physical: false, auto: true});
  });

  it('classifies table knobs: refused ones manual, index ones physical, the rest registry-only', () => {
    expect(changesOf((m) => m.tables.sample.securityMode = 'table')[0])
      .toMatchObject({kind: 'table-knob', detail: 'securityMode "row" -> (none)', physical: false, auto: false});
    expect(changesOf((m) => m.tables.sample_event.delegate = 'kind').map(flags)).toEqual([['table-knob', false, false]]);
    expect(changesOf((m) => m.tables.sample.businessKey = ['name', 'status']).map(flags))
      .toEqual([['table-knob', true, false]]);
    expect(changesOf((m) => m.tables.sample.idempotency = false).map(flags)).toEqual([['table-knob', true, false]]);
    expect(changesOf((m) => m.tables.sample.ginIndex = true).map(flags)).toEqual([['table-knob', true, false]]);
    expect(changesOf((m) => m.tables.sample.friendlyName = 'Samples').map(flags)).toEqual([['table-knob', false, true]]);
    expect(changesOf((m) => m.tables.sample_event.audit = true).map(flags)).toEqual([['table-knob', false, true]]);
  });

  it('reports tables, relations, filters, reorders, property schemas and extensibility', () => {
    expect(changesOf((m) => addRelation(m)).map((c) => [c.kind, c.table]))
      .toEqual([['add-table', 'label'], ['relations', 'sample'], ['add-table', 'sample_label']]);
    expect(changesOf((m) => delete m.tables.sample_event).map(flags)).toEqual([['drop-table', true, false]]);
    expect(changesOf((m) => m.tables.sample.filters = [{column: 'status'}]).map(flags)).toEqual([['filters', false, true]]);
    expect(changesOf((m) => {
      const {name, ...rest} = m.tables.sample.columns;
      m.tables.sample.columns = {...rest, name};
    }).map(flags)).toEqual([['reorder-columns', false, true]]);
    expect(changesOf((m) => m.propertySchemas.extra.note.type = 'int').map((c) => [c.kind, c.table, c.detail]))
      .toEqual([['property-schema', 'extra', 'changed']]);
    expect(changesOf((m) => m.extensible = {tables: true}).map(flags)).toEqual([['schema-extensible', false, true]]);
  });
});

describe('scaffold', () => {
  const from = buildSnapshot(fixture());

  it('writes up statements for manual changes only, and reverses every physical change in the down', () => {
    const to = buildSnapshot(fixture((m) => {
      delete m.tables.sample.columns.count;
      m.tables.sample.columns.note2 = {type: 'string'};
      m.tables.sample.columns.code = {type: 'string', required: true};
      m.tables.sample.columns.status.unique = true;
      m.tables.sample.columns.name.required = false;
      m.tables.sample.businessKey = ['name', 'status'];
      m.tables.sample_event.audit = true;
    }));
    const {up, down} = scaffold(diff(from, to), 'testdb', {from: from.hash, to: to.hash, title: 'testdb: t1',
      generatedBy: 'grok schema migrate'});

    // the machine-readable transition the deployer matches against the recorded snapshot
    expect(up.startsWith(`-- testdb: t1\n-- ems-migration: from=${from.hash} to=${to.hash}\n-- Generated by`)).toBe(true);
    expect(down.startsWith(`-- testdb: t1 (down)\n-- ems-migration: from=${to.hash} to=${from.hash}\n-- Generated by`))
      .toBe(true);
    const header = /^--\s*ems-migration:\s*from=(\S+)\s+to=(\S+)/m;
    expect(up.match(header)!.slice(1)).toEqual([from.hash, to.hash]);
    expect(down.match(header)!.slice(1)).toEqual([to.hash, from.hash]);
    const initial = scaffold(diff(emptySnapshot('testdb'), to), 'testdb', {to: to.hash});
    expect(initial.up).toContain(`\n-- ems-migration: from=none to=${to.hash}\n`);
    expect(initial.down).toContain(`\n-- ems-migration: from=${to.hash} to=none\n`);
    expect(up).toContain('-- testdb: t1\n');
    expect(up).toContain(`-- Generated by grok schema migrate from snapshot ${from.hash!.substring(0, 12)} ` +
      `to ${to.hash!.substring(0, 12)}.\n`);
    expect(up).toContain('-- Registry-only changes (no SQL; converged by the next deploy):\n' +
      '--   * table-knob sample_event: audit false -> (none)\n');
    expect(up).toContain('-- Applied by the deployer itself (no SQL here; reversed in the down script):\n' +
      '--   * change-unique sample.status: unique false -> true [physical]\n' +
      '--   * add-column sample.note2: added (string) [physical]\n');
    expect(up).toContain('-- data loss\nALTER TABLE testdb.sample DROP COLUMN IF EXISTS count;\n');
    expect(up).not.toContain('note2 text');
    expect(up).not.toContain('ux_sample_status');
    expect(up).toContain('-- TODO: a required column with no default cannot be added to a populated table\n' +
      'ALTER TABLE testdb.sample ADD COLUMN IF NOT EXISTS code text NOT NULL;\n');
    expect(up).toContain('ALTER TABLE testdb.sample ALTER COLUMN name DROP NOT NULL;\n');
    expect(up).toContain('DROP INDEX IF EXISTS testdb.ux_sample_business_key;\n' +
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_sample_business_key ON testdb.sample (name, status) WHERE NOT is_deleted;\n');

    expect(down).toContain('-- testdb: t1 (down)\n');
    expect(down).toContain(`from snapshot ${to.hash!.substring(0, 12)} to ${from.hash!.substring(0, 12)}.`);
    expect(down).not.toContain('Applied by the deployer');
    const order = [
      'ALTER TABLE testdb.sample DROP COLUMN IF EXISTS code;',
      'ALTER TABLE testdb.sample DROP COLUMN IF EXISTS note2;',
      'DROP INDEX IF EXISTS testdb.ux_sample_status;',
      '-- TODO: backfill NULLs in testdb.sample.name before SET NOT NULL\n' +
        'ALTER TABLE testdb.sample ALTER COLUMN name SET NOT NULL;',
      '-- data is not restored\nALTER TABLE testdb.sample ADD COLUMN IF NOT EXISTS count int;',
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_sample_business_key ON testdb.sample (name) WHERE NOT is_deleted;',
    ].map((s) => down.indexOf(s));
    expect(order.every((i) => i >= 0)).toBe(true);
    expect([...order].sort((a, b) => a - b)).toEqual(order);
  });

  it('renders CREATE TABLE in the engine vocabulary (from a drop-table down) and the FK DO block', () => {
    const withLabel = buildSnapshot(fixture((m) => m.tables.label = {
      idempotency: true, ginIndex: true, businessKey: ['name'],
      columns: {
        name: {type: 'string', required: true, unique: true},
        sample_id: {type: 'ref', ref: 'sample', onDelete: 'cascade'},
        owner: {type: 'user'},
        weight: {type: 'float', default: 1.5},
      }}));
    const {up, down} = scaffold(diff(withLabel, from), 'testdb');
    expect(up).toContain('-- drop-table label: dropped [physical] [manual]\n-- data loss\nDROP TABLE IF EXISTS testdb.label;\n');
    expect(down).toContain('-- data is not restored\nCREATE TABLE IF NOT EXISTS testdb.label (\n' +
      '  id uuid PRIMARY KEY DEFAULT uuid_generate_v4(),\n' +
      '  is_deleted bool NOT NULL DEFAULT false,\n' +
      '  version int NOT NULL DEFAULT 1,\n' +
      '  created_on timestamp without time zone NOT NULL DEFAULT now(),\n' +
      '  updated_on timestamp without time zone NOT NULL DEFAULT now(),\n' +
      '  author_id uuid,\n' +
      '  idempotency_key uuid,\n' +
      "  data jsonb NOT NULL DEFAULT '{}',\n" +
      '  name text NOT NULL,\n' +
      '  sample_id uuid,\n' +
      '  owner uuid,\n' +
      '  weight float8\n' +
      ');\n' +
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_label_name ON testdb.label (name) WHERE NOT is_deleted;\n' +
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_label_business_key ON testdb.label (name) WHERE NOT is_deleted;\n' +
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_label_idempotency_key ON testdb.label (idempotency_key) ' +
      'WHERE idempotency_key IS NOT NULL;\n' +
      'DO $$\nBEGIN\n' +
      "  IF NOT EXISTS (SELECT 1 FROM pg_constraint WHERE conname = 'fk_label_sample_id') THEN\n" +
      '    ALTER TABLE testdb.label ADD CONSTRAINT fk_label_sample_id FOREIGN KEY (sample_id) REFERENCES testdb.sample(id);\n' +
      '  END IF;\nEND$$;\n' +
      'CREATE INDEX IF NOT EXISTS ix_label_sample_id ON testdb.label (sample_id);\n' +
      'DO $$\nBEGIN\n' +
      "  IF NOT EXISTS (SELECT 1 FROM pg_constraint WHERE conname = 'fk_label_owner') THEN\n" +
      '    ALTER TABLE testdb.label ADD CONSTRAINT fk_label_owner FOREIGN KEY (owner) REFERENCES public.users(id);\n' +
      '  END IF;\nEND$$;\n' +
      'CREATE INDEX IF NOT EXISTS ix_label_owner ON testdb.label (owner);\n' +
      'CREATE INDEX IF NOT EXISTS ix_label_data ON testdb.label USING gin (data jsonb_path_ops);\n');
    // the additive direction is the deployer's: header only in the up, DROP TABLE in the down
    const added = scaffold(diff(from, withLabel), 'testdb');
    expect(added.up).toContain('--   * add-table label: created [physical]\n');
    expect(added.up).toContain('-- No manual statement.\n');
    expect(added.up).not.toContain('CREATE TABLE');
    expect(added.down).toContain('-- data loss\nDROP TABLE IF EXISTS testdb.label;\n');
  });

  it('renders type, ref and index-knob changes with their inverses', () => {
    const typed = scaffold(changesOf((m) => m.tables.sample.columns.count.type = 'float'), 'testdb');
    expect(typed.up).toContain('-- review: values that do not survive the cast fail the statement\n' +
      'ALTER TABLE testdb.sample ALTER COLUMN count TYPE float8 USING count::float8;\n');
    expect(typed.down).toContain('ALTER TABLE testdb.sample ALTER COLUMN count TYPE int USING count::int;\n');

    const ref = scaffold(changesOf((m) => {
      m.tables.label = {columns: {name: {type: 'string'}}};
      m.tables.sample_event.columns.sample_id.ref = 'label';
    }), 'testdb');
    expect(ref.up).toContain('ALTER TABLE testdb.sample_event DROP CONSTRAINT IF EXISTS fk_sample_event_sample_id;\n' +
      'DO $$\nBEGIN\n' +
      "  IF NOT EXISTS (SELECT 1 FROM pg_constraint WHERE conname = 'fk_sample_event_sample_id') THEN\n" +
      '    ALTER TABLE testdb.sample_event ADD CONSTRAINT fk_sample_event_sample_id FOREIGN KEY (sample_id) ' +
      'REFERENCES testdb.label(id);\n');
    expect(ref.down).toContain('REFERENCES testdb.sample(id);\n');

    const knobs = scaffold(changesOf((m) => {
      m.tables.sample.idempotency = false;
      m.tables.sample.ginIndex = true;
    }), 'testdb');
    expect(knobs.up).toContain('DROP INDEX IF EXISTS testdb.ux_sample_idempotency_key;\n' +
      'ALTER TABLE testdb.sample DROP COLUMN IF EXISTS idempotency_key;\n');
    expect(knobs.up).toContain('CREATE INDEX IF NOT EXISTS ix_sample_data ON testdb.sample USING gin (data jsonb_path_ops);\n');
    expect(knobs.down).toContain('DROP INDEX IF EXISTS testdb.ix_sample_data;\n');
    expect(knobs.down).toContain('ALTER TABLE testdb.sample ADD COLUMN IF NOT EXISTS idempotency_key uuid;\n' +
      'CREATE UNIQUE INDEX IF NOT EXISTS ux_sample_idempotency_key ON testdb.sample (idempotency_key) ' +
      'WHERE idempotency_key IS NOT NULL;\n');
  });
});

describe('grok schema', () => {
  function makePackage(): string {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-schema-'));
    fs.cpSync(fixtureDir, dir, {recursive: true});
    return dir;
  }

  function mutateManifest(dir: string, mutate: (manifest: any) => void): void {
    const manifestPath = path.join(dir, 'databases', 'testdb', 'schema.json');
    const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
    mutate(manifest);
    fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2));
  }

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

  const sealedPath = (dir: string) => path.join(dir, 'databases', 'testdb', 'snapshot.json');
  const schemaDir = (dir: string) => path.join(dir, 'databases', 'testdb');
  const sqlFiles = (dir: string) => fs.readdirSync(schemaDir(dir)).filter((f) => f.endsWith('.sql'));

  afterEach(() => vi.mocked(createClient).mockReset());

  it('check fails without a seal, seal writes the snapshot, check passes after it', async () => {
    const dir = makePackage();
    let res = await run(dir, 'check');
    expect(res.result).toBe(true);
    expect(res.exitCode).toBe(1);
    expect(res.output).toContain('no snapshot sealed in');

    res = await run(dir, 'seal');
    expect(res.exitCode).toBeUndefined();
    const sealed = fs.readFileSync(sealedPath(dir), 'utf8');
    const expected = buildSnapshot(fixture());
    expect(sealed).toBe(serialize(expected));
    expect(res.output).toContain(`Sealed ${path.join('databases', 'testdb', 'snapshot.json')} ` +
      `(${expected.hash!.substring(0, 12)})`);

    res = await run(dir, 'check');
    expect(res.exitCode).toBeUndefined();
    expect(res.output).toContain('testdb: up to date');

    mutateManifest(dir, (m) => m.tables.sample.columns.note2 = {type: 'string'});
    res = await run(dir, 'check');
    expect(res.exitCode).toBe(1);
    expect(res.output).toContain('add-column sample.note2: added (string) [physical]');
    res = await run(dir, 'diff');
    expect(res.exitCode).toBeUndefined();
    expect(res.output).toContain('testdb: 1 change(s)');
    expect(res.output).toContain('add-column sample.note2: added (string) [physical]');
    // --dir narrows to one schema directory; an unknown one is a usage error
    expect((await run(dir, 'check', {dir: 'databases/testdb'})).exitCode).toBe(1);
    expect((await run(dir, 'check', {dir: 'databases/nope'})).result).toBe(false);
  });

  it('migrate numbers the scripts after the existing ones, reseals, and skips registry-only changes', async () => {
    const dir = makePackage();
    await run(dir, 'seal');
    fs.writeFileSync(path.join(schemaDir(dir), '0003_legacy.sql'), 'select 1;');
    expect((await run(dir, 'migrate')).result).toBe(false);
    let res = await run(dir, 'migrate', {name: 'nothing'});
    expect(res.output).toContain('nothing to migrate');
    expect(sqlFiles(dir)).toEqual(['0003_legacy.sql']);

    mutateManifest(dir, (m) => m.tables.sample.friendlyName = 'Samples');
    res = await run(dir, 'migrate', {name: 'meta'});
    expect(res.output).toContain('nothing physical to migrate');
    expect(sqlFiles(dir)).toEqual(['0003_legacy.sql']);
    expect((await run(dir, 'check')).exitCode).toBeUndefined();

    mutateManifest(dir, (m) => m.tables.sample.columns.note2 = {type: 'string'});
    res = await run(dir, 'migrate', {name: 'add_note'});
    expect(res.result).toBe(true);
    expect(res.exitCode).toBeUndefined();
    // the up script sits next to schema.json (where the deployer runs NNNN_*.sql), the down one in down/
    expect(sqlFiles(dir)).toEqual(['0003_legacy.sql', '0004_add_note.sql']);
    const up1 = fs.readFileSync(path.join(schemaDir(dir), '0004_add_note.sql'), 'utf8');
    expect(up1).toContain('-- testdb: add_note\n');
    const sealedBefore = buildSnapshot(fixture((m) => m.tables.sample.friendlyName = 'Samples')).hash;
    const sealedAfter = buildSnapshot(fixture((m) => {
      m.tables.sample.friendlyName = 'Samples';
      m.tables.sample.columns.note2 = {type: 'string'};
    })).hash;
    expect(up1).toContain(`-- ems-migration: from=${sealedBefore} to=${sealedAfter}\n`);
    expect(up1).toContain('--   * add-column sample.note2: added (string) [physical]\n');
    expect(up1).toContain('-- No manual statement.\n');
    expect(fs.readFileSync(path.join(schemaDir(dir), 'down', '0004_add_note.sql'), 'utf8'))
      .toContain('ALTER TABLE testdb.sample DROP COLUMN IF EXISTS note2;\n');
    expect((await run(dir, 'check')).exitCode).toBeUndefined();

    mutateManifest(dir, (m) => delete m.tables.sample.columns.count);
    res = await run(dir, 'migrate', {name: 'drop-count'});
    expect(res.output).toContain(path.join('databases', 'testdb', '0005_drop-count.sql'));
    expect(fs.existsSync(path.join(schemaDir(dir), 'down', '0005_drop-count.sql'))).toBe(true);
    expect(fs.readFileSync(path.join(schemaDir(dir), '0005_drop-count.sql'), 'utf8'))
      .toContain('-- data loss\nALTER TABLE testdb.sample DROP COLUMN IF EXISTS count;\n');
    expect(JSON.parse(fs.readFileSync(sealedPath(dir), 'utf8')).hash)
      .toBe(buildSnapshot(JSON.parse(fs.readFileSync(path.join(dir, 'databases', 'testdb', 'schema.json'), 'utf8'))).hash);
  });

  it('--server diffs against the recorded snapshot; 404 means none recorded', async () => {
    const dir = makePackage();
    const get = vi.fn();
    vi.mocked(createClient).mockResolvedValue({get} as any);

    get.mockRejectedValue(Object.assign(new Error('Not found'), {apiError: {errorCode: 404}}));
    let res = await run(dir, 'check', {server: true});
    expect(res.exitCode).toBe(1);
    expect(res.output).toContain('no snapshot recorded on the server');
    expect(get).toHaveBeenCalledWith('/domains/schemas/testdb/snapshot');
    expect(vi.mocked(createClient)).toHaveBeenCalledWith(undefined);

    get.mockResolvedValue(JSON.parse(serialize(buildSnapshot(fixture()))));
    res = await run(dir, 'check', {server: 'dev'});
    expect(res.exitCode).toBeUndefined();
    expect(res.output).toContain('up to date');
    expect(vi.mocked(createClient)).toHaveBeenLastCalledWith('dev');

    mutateManifest(dir, (m) => m.tables.sample.columns.count.type = 'float');
    res = await run(dir, 'diff', {server: true});
    expect(res.output).toContain('change-type sample.count: type int -> float [physical] [manual]');

    get.mockRejectedValue(Object.assign(new Error('boom'), {apiError: {errorCode: 500}}));
    res = await run(dir, 'diff', {server: true});
    expect(res.result).toBe(true);
    expect(res.exitCode).toBe(1);
    expect(res.output).toContain('boom');
  });
});
