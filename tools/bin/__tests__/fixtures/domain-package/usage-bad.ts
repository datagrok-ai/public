// Negative tsc fixture: wrong-typed insert payloads must FAIL compilation.
import {testdbDb, SampleInsert} from './src/generated/db';

export async function bad(): Promise<void> {
  const ins: SampleInsert = {name: 'a', score: 1, count: 'not-a-number'};
  await testdbDb.sample.insert(ins);
  await testdbDb.sample.insert({active: 'yes'});
}
