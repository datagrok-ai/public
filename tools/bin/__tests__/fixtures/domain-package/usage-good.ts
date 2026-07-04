// Positive tsc fixture: well-typed usage of the generated domain clients must compile.
import {testdbDb, SampleInsert, SampleColumn} from './src/generated/db';

export async function good(): Promise<void> {
  const ins: SampleInsert = {name: 'a', count: 3, tags: ['x'], score: 0.5, idempotencyKey: 'k'};
  await testdbDb.sample.insert(ins);
  const rows = await testdbDb.sample.query({filter: 'count > 1'});
  const count: number | undefined = rows[0].count;
  const col: SampleColumn = 'measured_on';
  void count; void col;
  await testdbDb.sampleEvent.insert({sample_id: rows[0].id, kind: 'created'});
}
