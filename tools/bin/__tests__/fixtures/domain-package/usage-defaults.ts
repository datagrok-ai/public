// Positive tsc fixture: an insert may omit the required columns the server fills from their default.
import {testdbDb, SampleRow} from './src/generated/db';

export async function defaults(row: SampleRow): Promise<void> {
  await testdbDb.samples.insert({name: 'a'});
  await testdbDb.transaction([{op: 'insert', table: 'sample', values: {name: 'b'}}]);
  const prio: number = row.prio;
  const score: number = row.score;
  void prio; void score;
}
