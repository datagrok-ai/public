import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {DfExportEntry, prepareExportEntries} from '../project-export-data';

// Fixtures: Compute2:TestProjectExportSingleOut (int x -> res),
// Compute2:TestProjectExportMultiOut (dataframe df, int x -> res1, res2),
// Compute2:TestProjectExportFileIn (optional file, int x -> res) — see package.ts.
// The calls are only prepare()d, never run.
category('Project export data', () => {
  const mkDf = (name?: string | null) => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'a', [1, 2])]);
    if (name != null)
      df.name = name;
    return df;
  };

  const entry = (name: string, isInput: boolean, df: DG.DataFrame): DfExportEntry =>
    ({name, isInput, df, viewers: []});

  const script = (e: DfExportEntry) => e.df.getTag(DG.Tags.CreationScript) ?? null;
  const expr = (e: DfExportEntry) => script(e)!.split(' //')[0];

  test('single output: stamps variable name and creation script without accessor', async () => {
    const call = DG.Func.byName('Compute2:TestProjectExportSingleOut').prepare({x: 5});
    const [out] = prepareExportEntries(call, [entry('res', false, mkDf('Result'))]);
    expect(out.df.name, 'Result');
    expect(out.df.getTag(DG.Tags.VariableName), 'res');
    const s = script(out)!;
    expect(/^res = .*TestProjectExportSingleOut\(5\) \/\/\{"timestamp": \d+\}$/.test(s), true, s);
  });

  test('multi output with df input: accessors on outputs, input referenced by name, not stamped', async () => {
    const base = mkDf('Base Table');
    const call = DG.Func.byName('Compute2:TestProjectExportMultiOut').prepare({df: base, x: 2});
    const [inp, out1, out2] = prepareExportEntries(call, [
      entry('df', true, base), entry('res1', false, mkDf('R1')), entry('res2', false, mkDf('R2')),
    ]);
    expect(script(inp), null);
    expect(inp.df.getTag(DG.Tags.VariableName) ?? null, null);
    expect(expr(out1).endsWith('.res1'), true, expr(out1));
    expect(expr(out2).endsWith('.res2'), true, expr(out2));
    expect(expr(out1).includes('"Base Table"'), true, expr(out1));
    const callPart = (s: string) => s.slice(s.indexOf(' = ') + 3, s.lastIndexOf('.'));
    expect(callPart(expr(out1)), callPart(expr(out2)));
  });

  test('df input not exported: outputs publish as static snapshots', async () => {
    const call = DG.Func.byName('Compute2:TestProjectExportMultiOut').prepare({df: mkDf('Base'), x: 2});
    const [out1, out2] = prepareExportEntries(call, [
      entry('res1', false, mkDf('R1')), entry('res2', false, mkDf('R2')),
    ]);
    expect(script(out1), null);
    expect(script(out2), null);
  });

  test('unnamed and unsaved df input: clone falls back to the parameter name', async () => {
    const base = mkDf(); // fresh in-memory df, name is null, never uploaded
    expect(base.name == null || base.name === '', true, `unexpected name: ${base.name}`);
    const call = DG.Func.byName('Compute2:TestProjectExportMultiOut').prepare({df: base, x: 2});
    const [inp, out] = prepareExportEntries(call, [
      entry('df', true, base), entry('res1', false, mkDf('R1')),
    ]);
    expect(inp.df.name, 'df');
    expect(expr(out).includes('"df"'), true, expr(out));
    expect(expr(out).includes('null'), false, expr(out));
  });

  test('whitespace-only df name: clone falls back to the parameter name', async () => {
    const base = mkDf('   ');
    const call = DG.Func.byName('Compute2:TestProjectExportMultiOut').prepare({df: base, x: 2});
    const [inp] = prepareExportEntries(call, [entry('df', true, base)]);
    expect(inp.df.name, 'df');
  });

  test('name collision: dedup suffix, creation script references the deduped input name', async () => {
    const base = mkDf('Same');
    const call = DG.Func.byName('Compute2:TestProjectExportMultiOut').prepare({df: base, x: 2});
    const [inp, out] = prepareExportEntries(call, [
      entry('df', true, base), entry('res1', false, mkDf('Same')),
    ]);
    expect(inp.df.name, 'Same');
    expect(out.df.name, 'Same (2)');
    expect(expr(out).includes('"Same"'), true, expr(out));
  });

  test('quotes in df name: platform escaping is applied in the creation script', async () => {
    const base = mkDf('My "Q" Table');
    const call = DG.Func.byName('Compute2:TestProjectExportMultiOut').prepare({df: base, x: 2});
    const [inp, out] = prepareExportEntries(call, [
      entry('df', true, base), entry('res1', false, mkDf('R1')),
    ]);
    expect(inp.df.name, 'My "Q" Table');
    expect(expr(out).includes('My ^^Q^^ Table'), true, expr(out));
  });

  test('unset optional file input keeps sync on, a set one disables it', async () => {
    const unset = DG.Func.byName('Compute2:TestProjectExportFileIn').prepare({x: 1});
    const [out1] = prepareExportEntries(unset, [entry('res', false, mkDf('R'))]);
    expect(script(out1) != null, true);

    const set = DG.Func.byName('Compute2:TestProjectExportFileIn').prepare({
      file: DG.FileInfo.fromString('f.txt', 'data'), x: 1,
    });
    const [out2] = prepareExportEntries(set, [entry('res', false, mkDf('R'))]);
    expect(script(out2), null);
  });
});
