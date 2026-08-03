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
  // builder negatives: every chain below must fail on the typed client
  testdbDb.sample.query().where('nope', '=', 1);
  testdbDb.sample.query().orderBy('nope');
  testdbDb.sample.query().select('nope');
  testdbDb.sample.query().expand('details:nope');
  const slim = await testdbDb.sample.query().select('name').top(1);
  const notSelected: number | undefined = slim[0].count;
  const [insR] = await testdbDb.transaction([
    {op: 'insert', table: 'sample', values: {name: 'x', score: 1}}]);
  const badField: boolean = insR.deleted;
  void notSelected; void badField;
}
