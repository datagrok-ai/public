// Positive tsc fixture: well-typed usage of the generated typed UI wrappers must compile.
import * as DG from 'datagrok-api/dg';
import {DomainAppView, DomainGrid, EntityListWidget} from '@datagrok-libraries/domain-ui';
import {sampleUi, sampleEventUi, SampleQuerySpec} from './src/generated/db-ui';

export async function goodUi(): Promise<void> {
  // a whole browse/CRUD page, one expression — with a compile-checked query
  const view: DomainAppView = sampleUi.appView({query: {filter: 'count > 1', columns: ['name', 'status']}});
  const entity = sampleUi.entityView('some-id');
  void view; void entity;

  // grids and lists: columns, expand keys and insert defaults are this table's
  const grid: DomainGrid = await sampleUi.grid({query: {columns: ['name', 'measured_on']}, editable: false});
  const list: EntityListWidget | null = await sampleUi.list({mode: 'grid',
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
  const table: string = sampleUi.table;
  void query; void table;
}
