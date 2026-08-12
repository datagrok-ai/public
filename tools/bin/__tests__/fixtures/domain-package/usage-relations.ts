// Positive tsc fixture (many-to-many relations): compiled only by the relation tests,
// which add a `label` target and a `sample_label` junction to the manifest first.
import {testdbDb, SampleUpdate, TestdbTransactionOp} from './src/generated/db';

export async function relations(): Promise<void> {
  // read side: the relation is an expand key, and its fields intersect into the row
  const rows = await testdbDb.samples.query().expand('labels');
  const link = rows[0].labels?.[0];
  const linkId: string | undefined = link?.id;
  const linkName: string | undefined = link?.name;
  void linkId; void linkName;
  await testdbDb.samples.query({expand: ['labels']});
  // filtering travels the relation as a path, in the tree and in the string form alike
  await testdbDb.samples.query({filter: [{property: 'labels.name', operator: '=', value: 'bug'}]});
  await testdbDb.samples.query({filter: 'labels.name = "bug"'});
  await testdbDb.samples.facets({facets: [{id: 'l', kind: 'categories', column: 'labels.id'}]});
  // write side: an insert carries the link set...
  await testdbDb.samples.insert({name: 'a', score: 1, labels: ['id1', 'id2']});
  // ...a patch replaces the links the caller can see ([] clears them). The client
  // takes <Table>Update as its fifth generic, so both the variable and the object
  // literal (excess-property checked) go straight into update()
  const patch: SampleUpdate = {count: 2, labels: []};
  await testdbDb.samples.update('x', patch);
  await testdbDb.samples.update('x', {labels: ['id1'], name: 'a'}, {version: 3});
  // ...and a transaction links a target created in the same transaction, element-wise
  const ops: TestdbTransactionOp[] = [
    {op: 'insert', table: 'label', ref: 'l', values: {name: 'bug'}},
    {op: 'update', table: 'sample', id: 'x', values: {labels: ['$l'], count: 1}},
  ];
  await testdbDb.transaction(ops);
}
