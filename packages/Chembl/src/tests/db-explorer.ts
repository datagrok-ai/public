import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {DBExplorerRenderer} from '@datagrok-libraries/db-explorer/src/renderer';
import {SchemaAndConnection} from '@datagrok-libraries/db-explorer/src/types';

category('DB Explorer', () => {
  test('Links pane for columns without DbColumn tags', async () => {
    const schema = {
      schema: {
        references: {},
        referencedBy: {public: {molecule_dictionary: {molregno: [
          {refSchema: 'public', refTable: 'compound_records', refColumn: 'molregno'}]}}},
      },
      connection: null,
    } as unknown as SchemaAndConnection;
    const schemaPromise = () => Promise.resolve(schema);
    const renderer = new DBExplorerRenderer(schemaPromise, 'public', {valueConverter: (v) => v, joinOptions: []});
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('molregno', new Int32Array([22942]))]);
    df.col('molregno')!.setTag(DG.Tags.DbTable, 'molecule_dictionary');
    const acc = ui.accordion('CHEMBL123');
    await renderer.renderAssociations(acc, schemaPromise, 'public', 'molecule_dictionary', df);
    const links = acc.getPane('Links');
    links.expanded = true;
    expect(links.root.textContent?.includes('public.molecule_dictionary.molregno'), true);
  });
});
