import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/test/src/test';

/* The enrichments this package ships (enrichments/*.json) are written by the server into PowerPack's folder
   when the package is published. Each is applied to the first rows of its key table the way the Enrich pane
   applies it, so a column the ChEMBL schema no longer has fails here: PowerPack reports a failed query in a
   balloon and adds no column. */
const ROOT = 'System:AppData/PowerPack/enrichments/Chembl_Chembl/';

category('Enrichments', () => {
  test('Every shipped enrichment applies', async () => {
    const run = DG.Func.find({package: 'PowerPack', name: 'runEnrichment'})[0];
    expect(run != null, true, 'PowerPack:runEnrichment is not registered');
    const conn = (await grok.dapi.connections.filter('namespace = "Chembl:" and shortName = "Chembl"').list())[0];
    expect(conn != null, true, 'no Chembl:Chembl connection');
    const files = (await grok.dapi.files.list(ROOT, true)).filter((f) => f.isFile && f.name.endsWith('.json'));
    expect(files.length > 0, true, `no enrichment is deployed under ${ROOT}`);
    const failures: string[] = [];
    for (const file of files) {
      // <db>/<schema>/<table>/<column>/<name>.json
      const [db, schema, table, column, fileName] = file.fullPath.slice(file.fullPath.indexOf(ROOT) + ROOT.length).split('/');
      const name = fileName.slice(0, -'.json'.length);
      const query = DG.TableQuery.create(conn);
      query.table = `${schema}.${table}`;
      query.fields = [`${schema}.${table}.${column}`];
      query.limit = 20;
      const df = await query.executeTable();
      const before = df.columns.length;
      await run.prepare({conn, schema, table, column, name, df, db}).call();
      if (df.columns.length === before)
        failures.push(`"${name}" on ${table}.${column} added no column`);
    }
    expect(failures.join('; '), '');
  }, {timeout: 300000});
});
