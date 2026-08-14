// Positive tsc fixture (many-to-many relations + typed UI): compiled only by the
// relation --ui test, which adds a `label` target and a `sample_label` junction first.
import {DomainTable} from '@datagrok-libraries/domain-ui';
import {sampleUi, labelUi, testdbUiDb} from './src/generated/db-ui';
import type {SampleColumn, SampleExpand, SampleInsert, SampleRow,
  SampleUpdate} from './src/generated/db';

export async function relationsUi(): Promise<void> {
  // the handle carries the update payload, so its client takes the link sets —
  // as a variable AND as an object literal (excess-property checked)
  const samples: DomainTable<SampleRow, SampleInsert, SampleColumn, SampleExpand, SampleUpdate> =
    await sampleUi.table();
  const patch: SampleUpdate = {count: 2, labels: []};
  await samples.client.update('x', patch);
  await sampleUi.client.update('x', {labels: ['id1'], name: 'a'});
  await sampleUi.client.insert({name: 'a', score: 1, labels: ['id1']});
  // the relation is an expand key of every spec the widgets take
  await sampleUi.grid({query: {expand: ['labels']}});
  await sampleUi.list({query: {expand: ['labels'], limit: 10}});
  // the schema handle's per-table property carries the same typing
  const db = await testdbUiDb();
  await db.samples.client.update('x', {labels: ['id1']});
  // a relation-less table keeps the Partial<Row> default
  await labelUi.client.update('y', {name: 'bug'});
}
