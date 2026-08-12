// Positive tsc fixture: well-typed usage of the generated typed UI wrappers must compile.
import * as DG from 'datagrok-api/dg';
import {DomainAppView, DomainForm, DomainGrid, DomainTable,
  EntityListWidget} from '@datagrok-libraries/domain-ui';
import {sampleUi, sampleEventUi, testdbUiDb, TestdbUiDb, SampleQuerySpec} from './src/generated/db-ui';
import type {SampleColumn, SampleExpand, SampleInsert, SampleRow} from './src/generated/db';

export async function goodUi(): Promise<void> {
  // UI-FACADE §1.7: the typed handle is the ONE await, and every factory on it is
  // synchronous — the way a page that builds several widgets pays for one prefetch
  const samples: DomainTable<SampleRow, SampleInsert, SampleColumn, SampleExpand> =
    await sampleUi.table();
  const onHandle: DomainForm = samples.form({values: {name: 'a'}});
  const appOnHandle: DomainAppView = samples.app();
  const handleRows: SampleRow[] = await samples.query({columns: ['name', 'status']});
  void onHandle; void appOnHandle; void handleRows;

  // the shortcuts acquire a handle of their own — with this table's insert payload,
  // columns and expand keys compile-checked
  const form: DomainForm = await sampleUi.form({values: {name: 'a', status: 'new', count: 2}});
  const saved: boolean = await sampleUi.formDialog({values: {name: 'a'}, title: 'New sample'});
  void form; void saved;

  const view: DomainAppView = await sampleUi.app({query: {filter: 'count > 1', columns: ['name', 'status']}});
  const listPage: DomainAppView = await sampleUi.listView({query: {columns: ['name']}});
  const entity = await sampleUi.entityView('some-id');
  void view; void listPage; void entity;

  const grid: DomainGrid = await sampleUi.grid({query: {columns: ['name', 'measured_on']}, editable: false});
  const list: EntityListWidget = await sampleUi.list({mode: 'grid',
    query: {expand: ['details:sample_event'], limit: 10}});
  void grid; void list;
  // a detail grid prefilling its parent FK — the column is checked against the child's insert
  await sampleEventUi.grid({defaults: {sample_id: 'parent-id', kind: 'created'}});
  const spec: SampleQuerySpec = {filter: {property: 'status', operator: '=', value: 'new'},
    columns: ['measured_on'], expand: ['details:sample_event']};
  void spec;

  // the typed client is the one from db.ts — rows keep their row type
  const rows = await sampleUi.client.query({columns: ['name', 'count']});
  const name: string = rows[0].name;
  const row: DG.DomainRow = sampleUi.row(rows[0]);
  void name; void row;

  // dialogs and pickers go through the table's registered handler
  const handler: DG.DomainObjectHandler = sampleUi.handler();
  await sampleUi.edit(row);
  const picked: DG.DomainRow | null = await sampleUi.pick();
  void handler; void picked;

  // a serializable query with the schema and table implied
  const query: DG.DomainQuery = sampleUi.query({filters: ['status = "new"'], columns: ['name'], limit: 10});
  const address: string = sampleUi.address;
  void query; void address;

  // the schema-level handle: EVERY table prefetched under one await, each typed
  const db: TestdbUiDb = await testdbUiDb();
  const dbForm: DomainForm = db.samples.form({values: {name: 'a', status: 'new'}});
  const dbGrid: DomainGrid = db.sampleEvents.grid({defaults: {sample_id: 'parent-id', kind: 'created'}});
  const dbApp: DomainAppView = db.samples.app({query: {columns: ['name', 'status']}});
  const dbRows: SampleRow[] = await db.samples.query({columns: ['name']});
  const byName: DomainTable = db.table('sample_event');
  void dbForm; void dbGrid; void dbApp; void dbRows; void byName;
}
