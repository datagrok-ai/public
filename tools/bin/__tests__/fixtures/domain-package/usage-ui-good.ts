// Positive tsc fixture (typed u2 handles): compiled by the --ui test.
import type {DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';
import {getTestdbDb, TestdbDb} from './src/generated/db-ui';
import type {SampleRow, SampleEventRow} from './src/generated/db';

export async function goodUi(): Promise<void> {
  // ONE await, then every handle is typed by its row type from db.ts
  const db: TestdbDb = await getTestdbDb();
  const schema: 'testdb' = db.schema;
  const samples: DomainTable<SampleRow> = db.tables.samples;
  const events: DomainTable<SampleEventRow> = db.tables.sampleEvents;
  void schema;

  // the data clients are the ones db.ts exports — rows keep their row type
  const rows: SampleRow[] = await db.data.samples.query({columns: ['name', 'count']});
  const name: string = rows[0].name;
  void name;

  // the registries check the row type: actions, validators, the renderer, sources and drafts
  samples.actions.add({name: 'Bump', requires: 'edit', when: (r) => r.count !== undefined, run: (r) => {
    const n: number | undefined = r.count;
    r.name = `${r.name} (${n})`;
  }});
  samples.validators.add('name', (v, r) => r.status === 'new' && String(v).length < 3 ? 'Too short' : null);
  samples.renderer = {...samples.renderer, caption: (r) => `${r.name}: ${r.count ?? 0}`};
  const src = samples.source({query: 'status = "new"', defaults: {status: 'new'}});
  const draft = samples.draft({name: 'a', count: 1});
  const created: SampleRow = draft.newRow({name: 'b'});
  void src;
  events.draft({sample_id: created.id, kind: 'created'});

  // the same handle again is the same object — the cache is per page
  const again: TestdbDb = await getTestdbDb();
  void again;
}
