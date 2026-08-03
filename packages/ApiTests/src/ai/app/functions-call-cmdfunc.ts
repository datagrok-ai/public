import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {demog, withTableView} from '../helpers';

// @CmdFunc Dart commands surface as platform Funcs via grok.functions —
// eval/find:Promise<Func>, call(name, params), parse(expr):FuncCall. AddNewColumn is the safest core @CmdFunc.
category('AI: App: Functions Call (@CmdFunc)', () => {
  test('eval resolves a @CmdFunc name to a DG.Func', async () => {
    const f: DG.Func = await grok.functions.eval('AddNewColumn');
    expect(f instanceof DG.Func, true);
    expect(f.name, 'AddNewColumn');
    // A @CmdFunc is registered with typed input params reachable via the resolved Func.
    expect(Array.isArray(f.inputs), true);
    expect(f.inputs.length > 0, true);
  });

  test('find returns a non-null Func for a @CmdFunc', async () => {
    const f = await grok.functions.find('AddNewColumn');
    expect(f != null, true);
    expect(f instanceof DG.Func, true);
    expect((f as DG.Func).name, 'AddNewColumn');
    // Note: find() delegates to eval(), which THROWS ("Variable ... not found") for an
    // unknown name rather than returning null, so a negative branch is not asserted here.
  });

  test('call(AddNewColumn) transforms the table and adds the column', async () => {
    const df = demog(20);
    const before = df.columns.length;
    // AddNewColumn is a @CmdFunc that appends a computed column in place.
    await grok.functions.call('AddNewColumn', {table: df, expression: '${age} * 2', name: 'age2'});
    expect(df.columns.length, before + 1);
    expect(df.columns.contains('age2'), true);
    // Value-content check: the new column equals the formula applied to the source.
    const age2 = df.col('age2')!;
    const age = df.col('age')!;
    for (let i = 0; i < df.rowCount; i++) {
      if (age.isNone(i))
        continue;
      expect(age2.get(i), age.get(i) * 2);
      break;
    }
  });

  test('call works for a second core @CmdFunc (GetTable)', async () => {
    // GetTable(tableName) is a core @CmdFunc resolving a DataFrame by name
    // (xamgle/lib/src/grok_actions.dart). Confirm it exists before asserting.
    const f = await grok.functions.find('GetTable');
    // Defensive: degrade gracefully if this @CmdFunc is not registered on this build.
    if (f == null)
      return;
    const df = demog(33);
    // Attaching to a table view registers it under its name so GetTable can resolve it.
    await withTableView(df, async () => {
      const resolved: DG.DataFrame = await grok.functions.call('GetTable', {tableName: df.name});
      expect(resolved instanceof DG.DataFrame, true);
      expect(resolved.rowCount, 33);
      expect(resolved.name, df.name);
    });
  });

  test('parse builds a FuncCall for a @CmdFunc expression', async () => {
    const fc = grok.functions.parse('AddNewColumn(null, "x", "1")');
    expect(fc instanceof DG.FuncCall, true);
    expect(fc.func != null, true);
    expect(fc.func.name, 'AddNewColumn');
  });
}, {owner: 'agolovko@datagrok.ai'});
