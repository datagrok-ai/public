// Positive tsc fixture: well-typed usage of the generated domain clients must compile.
import * as grok from 'datagrok-api/grok';
import type {Dayjs} from 'dayjs';
import {testdbDb, SampleRow, SampleInsert, SampleColumn, TestdbTransactionOp} from './src/generated/db';

export async function good(): Promise<void> {
  const ins: SampleInsert = {name: 'a', count: 3, tags: ['x'], score: 0.5, idempotencyKey: 'k',
    status: 'new', measured_on: '2026-01-01T00:00:00Z'};
  await testdbDb.sample.insert(ins);
  const rows = await testdbDb.sample.query({filter: 'count > 1', columns: ['name', 'count']});
  const count: number | undefined = rows[0].count;
  const measured: Dayjs | undefined = rows[0].measured_on;
  const created: Dayjs = rows[0].created_on;
  const col: SampleColumn = 'measured_on';
  void count; void col; void measured; void created;
  // expand keys are compile-checked per table
  await testdbDb.sample.query({expand: ['details:sample_event']});
  await testdbDb.sampleEvent.query({expand: ['sample_id']});
  await testdbDb.sampleEvent.insert({sample_id: rows[0].id, kind: 'created'});
  // typed transaction ops pair values with their table; $refs ride as strings
  const ops: TestdbTransactionOp[] = [
    {op: 'insert', table: 'sample', ref: 's', values: {name: 'tx', score: 1}},
    {op: 'insert', table: 'sample_event', values: {sample_id: '$s', kind: 'created'}},
    {op: 'delete', table: 'sample', id: '$s'},
  ];
  await testdbDb.transaction(ops);
  // single-generic (row-only) clients stay permissive on insert — backward compat
  await grok.dapi.domains.table<SampleRow>('testdb.sample').insert({name: 'z', idempotencyKey: 'k'});
  await grok.dapi.domains.table('testdb.sample').insert({});
}
