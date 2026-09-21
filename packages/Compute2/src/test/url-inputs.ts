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
      const {patch, warnings} = await parseUrlInputs(call, new URLSearchParams(`df=${id}`));
      expect(warnings.length, 0, warnings.join('; '));
      expect(patch.get('df') instanceof DG.DataFrame, true);
      expect(patch.get('df').rowCount, 10);
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

  test('buildInputsUrl serializes scalars and current entity ids only', async () => {
    const call = prepare();
    applyUrlInputs(call, new Map<string, any>(Object.entries({a: 5, b: 2.5, flag: false, s: 'txt'})));
    const {url, skipped} = buildInputsUrl(call);
    const params = new URL(url).searchParams;
    expect(params.get('a'), '5');
    expect(params.get('b'), '2.5');
    expect(params.get('flag'), 'false');
    expect(params.get('s'), 'txt');
    expect(params.get('df'), null);
    expect(skipped.length, 0);

    // a local dataframe has no entity id and is skipped
    call.inputs['df'] = grok.data.demo.demog(1);
    expect(buildInputsUrl(call).skipped.includes('df'), true);
  });

  test('buildInputsUrl uses the current entity id of a loaded table', async () => {
    const id = await grok.dapi.tables.uploadDataFrame(grok.data.demo.demog(5));
    try {
      const call = prepare();
      applyUrlInputs(call, new Map<string, any>(Object.entries(
        {a: 1, b: 1.5, flag: true, s: 'x', df: await grok.dapi.tables.getTable(id)})));
      expect(new URL(buildInputsUrl(call).url).searchParams.get('df'), id);
      // replacing the table drops it from the link instead of reusing a stale id
      call.inputs['df'] = grok.data.demo.demog(1);
      const rebuilt = buildInputsUrl(call);
      expect(new URL(rebuilt.url).searchParams.get('df'), null);
      expect(rebuilt.skipped.includes('df'), true);
    } finally {
      const tableInfo = await grok.dapi.tables.find(id);
      if (tableInfo)
        await grok.dapi.tables.delete(tableInfo);
    }
  });
});
