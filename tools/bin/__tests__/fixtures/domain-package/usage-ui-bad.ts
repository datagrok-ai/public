// Negative tsc fixture: every wrong-typed usage below must FAIL compilation.
import {getTestdbDb} from './src/generated/db-ui';

export async function badUi(): Promise<void> {
  const db = await getTestdbDb();
  const schema: 'other' = db.schema;
  const missing = db.tables.nope;
  db.tables.samples.actions.add({name: 'x', run: (r) => { const n: number = r.name; void n; }});
  db.tables.samples.validators.add('name', (v, r) => { const k: number = r.status; void k; return null; });
  db.tables.samples.draft({count: 'not-a-number'});
  db.tables.sampleEvents.draft({sample_id: 42});
  await db.data.samples.query({columns: ['nope']});
  void schema; void missing;
}
