// Negative tsc fixture: every wrong-typed usage below must FAIL compilation.
import {testdbDb, SampleInsert, TestdbTransactionOp} from './src/generated/db';

export async function bad(): Promise<void> {
  const ins: SampleInsert = {name: 'a', score: 1, count: 'not-a-number'};
  await testdbDb.sample.insert(ins);
  await testdbDb.sample.insert({name: 'a', score: 1, active: 'yes'});
  await testdbDb.sample.insert({});
  await testdbDb.sample.insert({name: 'a', score: 1, status: 'bogus'});
  await testdbDb.sample.query({expand: ['details:nope']});
  await testdbDb.sample.query({columns: ['nope']});
  await testdbDb.sample.aggregate({groupBy: ['nope'], measures: [{fn: 'count'}]});
  const rows = await testdbDb.sample.query({});
  rows[0].measured_on = '2026-01-01';
  const ops: TestdbTransactionOp[] = [
    {op: 'insert', table: 'sample_event', values: {name: 'wrong-table-values', score: 1}},
  ];
  await testdbDb.transaction(ops);
}
