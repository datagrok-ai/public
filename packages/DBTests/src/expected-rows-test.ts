import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from './package-test';

category('Queries: expected rows', () => {
  test('Every query returns its meta.testExpectedRows', async () => {
    const queries = (await grok.dapi.queries.filter(`package.id = "${_package.id}"`)
      .include('params,connection').list()).filter((q) => q.options['testExpectedRows'] != null);
    expect(queries.length > 0, true, 'no query declares meta.testExpectedRows');
    const reachable: {[connectionId: string]: string} = {};
    const unreachable = new Set<string>();
    const wrong: string[] = [];
    for (const q of queries) {
      const connection = q.connection;
      reachable[connection.id] ??= await (await grok.dapi.connections.find(connection.id)).test();
      if (reachable[connection.id] !== 'ok') {
        unreachable.add(`${connection.name}: ${reachable[connection.id]}`);
        continue;
      }
      try {
        // the counts are for the inputs' defaults, which a call does not evaluate by itself
        const call = q.prepare();
        for (const p of q.inputs.filter((p) => p.initialValue != null))
          call.setParamValue(p.name, await call.evalParamDefault(p.name));
        const rows = (await call.call()).getOutputParamValue().rowCount;
        if (rows !== Number(q.options['testExpectedRows']))
          wrong.push(`${q.name} returned ${rows} rows, expected ${q.options['testExpectedRows']}`);
      } catch (e: any) {
        wrong.push(`${q.nqName} failed: ${String(e?.message ?? e).replace(/\s+/g, ' ')}`);
      }
    }
    expect(wrong.join(' | '), '');
    return unreachable.size ? `Skipped the queries of unreachable connections: ${[...unreachable].join('; ')}` : 'OK';
  }, {timeout: 120000});

  test('PostgresByStringChoices filters by the chosen country', async () => {
    for (const [country, rows] of [['France', 77], ['USA', 122]] as [string, number][]) {
      const df: DG.DataFrame = await grok.data.query('DbTests:PostgresByStringChoices', {shipCountry: country});
      expect(df.rowCount, rows, country);
      expect(df.col('shipcountry')!.toList().every((v) => v === country), true, `every shipcountry is ${country}`);
    }
  });
});
