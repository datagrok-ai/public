import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

category('ENA', () => {
  category('detector', () => {
    test('detects EnaID semType', async () => {
      const df = DG.DataFrame.fromCsv(
        `id
AA046425
AA046426
BX927145`
      );
      await grok.data.detectSemanticTypes(df);
      expect(df.col('id')!.semType, 'EnaID');
    });

    test('does not detect non-ENA strings', async () => {
      const df = DG.DataFrame.fromCsv(
        `id
hello
world`
      );
      await grok.data.detectSemanticTypes(df);
      expect(df.col('id')!.semType !== 'EnaID', true);
    });

    test('does not detect partial matches', async () => {
      const df = DG.DataFrame.fromCsv(
        `id
AA04642
AA0464255`
      );
      await grok.data.detectSemanticTypes(df);
      expect(df.col('id')!.semType !== 'EnaID', true);
    });
  });

  category('info panel', () => {
    test('enaSequence panel returns widget', async () => {
      const fn = DG.Func.find({name: 'enaSequence', package: 'valerij_developer_exercise_1'})[0];
      expect(fn != null, true);

      const result = await fn.apply({cellText: 'AA046425'});
      const widget = result.getOutputParamValue();
      expect(widget instanceof DG.Widget, true);
    });
  });

  category('fetchENASequence', () => {
    test('formENADataTable function is registered', async () => {
      const fn = DG.Func.find({name: 'formENADataTable', package: 'valerij_developer_exercise_1'})[0];
      expect(fn != null, true);
    });

    test('fetches sequences via text search', async () => {
      const url = 'https://www.ebi.ac.uk/ena/browser/api/embl/textsearch?result=sequence&query=coronavirus&limit=3&offset=0';
      const resp = await grok.dapi.fetchProxy(url);
      const text = await resp.text();
      expect(text.length > 0, true);
      expect(text.includes('ID'), true);
      expect(text.includes('SQ'), true);
    });

    test('ENA fasta endpoint returns data', async () => {
      const url = 'https://www.ebi.ac.uk/ena/browser/api/fasta/AA046425';
      const resp = await grok.dapi.fetchProxy(url);
      const text = await resp.text();
      expect(text.startsWith('>'), true);
      expect(text.length > 50, true);
    });
  });
});
