import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseUrlInputs, applyUrlInputs, missingMandatoryInputs, buildInputsUrl} from '../url-inputs';

// Fixture: Compute2:TestUrlInputsFixture (int a, double b, bool flag, string s, datetime when,
// dataframe df, optional double opt, nullable int nul) — see package.ts.
category('URL inputs', () => {
  const prepare = () => DG.Func.byName('Compute2:TestUrlInputsFixture').prepare();

  test('parses scalars, skips unmatched params, warns on unparsable values', async () => {
    const call = prepare();
    const {patch, warnings} = await parseUrlInputs(call,
      new URLSearchParams('a=5&b=2.5&flag=1&s=hello&bogus=1&nul=abc'));
    expect(patch.get('a'), 5);
    expect(patch.get('b'), 2.5);
    expect(patch.get('flag'), true);
    expect(patch.get('s'), 'hello');
    // unmatched params (platform params like q/layout/browse) are skipped without a warning
    expect(patch.has('bogus'), false);
    expect(patch.has('nul'), false);
    expect(warnings.length, 1, `expected a warning only for nul, got: ${warnings.join('; ')}`);
  });

  test('rejects non-integer and empty values for numeric inputs', async () => {
    const call = prepare();
    const {patch, warnings} = await parseUrlInputs(call, new URLSearchParams('a=5.5&b='));
    expect(patch.size, 0);
    expect(warnings.length, 2);
  });

  test('parses datetime into dayjs', async () => {
    const call = prepare();
    const {patch, warnings} = await parseUrlInputs(call, new URLSearchParams('when=2026-01-02T03:04:05Z'));
    expect(warnings.length, 0);
    expect(dayjs.isDayjs(patch.get('when')), true);
    expect(patch.get('when').toISOString(), '2026-01-02T03:04:05.000Z');
  });

  test('skips the id param silently', async () => {
    const call = prepare();
    const {patch, warnings} = await parseUrlInputs(call, new URLSearchParams('id=123'));
    expect(patch.size, 0);
    expect(warnings.length, 0);
  });

  test('loads dataframe inputs by entity id', async () => {
    const df = grok.data.demo.demog(10);
    const id = await grok.dapi.tables.uploadDataFrame(df);
    try {
      const call = prepare();
      const {patch, entityIds, warnings} = await parseUrlInputs(call, new URLSearchParams(`df=${id}`));
      expect(warnings.length, 0, warnings.join('; '));
      expect(patch.get('df') instanceof DG.DataFrame, true);
      expect(patch.get('df').rowCount, 10);
      expect(entityIds.get('df'), id);
    } finally {
      const tableInfo = await grok.dapi.tables.find(id);
      if (tableInfo)
        await grok.dapi.tables.delete(tableInfo);
    }
  });

  test('warns when a table id cannot be loaded', async () => {
    const call = prepare();
    const {patch, warnings} = await parseUrlInputs(call,
      new URLSearchParams('df=00000000-0000-0000-0000-000000000000'));
    expect(patch.size, 0);
    expect(warnings.length, 1);
  });

  test('missingMandatoryInputs exempts optional but not nullable', async () => {
    const call = prepare();
    const missing = missingMandatoryInputs(call);
    expect(missing.includes('a'), true);
    expect(missing.includes('df'), true);
    expect(missing.includes('opt'), false);
    // nullable defaults to true for non-strings, so it does not exempt an input
    expect(missing.includes('nul'), true);
    applyUrlInputs(call, new Map<string, any>(Object.entries(
      {a: 1, b: 1.5, flag: true, s: 'x', when: dayjs(), df: grok.data.demo.demog(1), nul: 7})));
    expect(missingMandatoryInputs(call).length, 0);
  });

  test('buildInputsUrl serializes scalars and known entity ids only', async () => {
    const call = prepare();
    applyUrlInputs(call, new Map<string, any>(Object.entries({a: 5, b: 2.5, flag: false, s: 'txt'})));
    const {url, skipped} = buildInputsUrl(call, new Map([['df', 'table-id-1']]));
    const params = new URL(url).searchParams;
    expect(params.get('a'), '5');
    expect(params.get('b'), '2.5');
    expect(params.get('flag'), 'false');
    expect(params.get('s'), 'txt');
    expect(params.get('df'), null);
    expect(skipped.length, 0);

    call.inputs['df'] = grok.data.demo.demog(1);
    const withValueNoId = buildInputsUrl(call, new Map());
    expect(withValueNoId.skipped.includes('df'), true);
    const withId = buildInputsUrl(call, new Map([['df', 'table-id-1']]));
    expect(new URL(withId.url).searchParams.get('df'), 'table-id-1');
  });
});
