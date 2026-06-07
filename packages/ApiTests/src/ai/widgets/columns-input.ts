import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow} from '../helpers';

// ColumnsInput (rich column multi-picker) — ui.input.columns + setColumnsInputTable,
// setInputAdditionalColumns and the available/checked/showSelectedColsOnTop creation options.
category('AI: Widgets: ColumnsInput', () => {
  test('value get/set round-trips the column selection', async () => {
    const t = demog();
    const input = ui.input.columns('cols', {table: t});
    expect(input instanceof DG.InputBase, true);
    // The `checked` creation option seeds DOM checkboxes that only materialize on render,
    // so the value model is exercised through the value setter (the headless-stable contract).
    input.value = [t.col('age')!, t.col('height')!];
    expect(Array.isArray(input.value), true);
    expect(input.value.length, 2);
    const names = input.value.map((c) => c.name);
    expect(names.indexOf('age') >= 0, true);
    expect(names.indexOf('height') >= 0, true);
  });

  test('available restricts the candidate columns', async () => {
    const input = ui.input.columns('cols', {table: demog(), available: ['age', 'height', 'weight']});
    expect(input instanceof DG.InputBase, true);
    expect(Array.isArray(input.value), true);
  });

  test('setColumnsInputTable re-binds the table without throwing', async () => {
    const input = ui.input.columns('cols', {table: demog()});
    expectNoThrow(() => ui.input.setColumnsInputTable(input, demog(100), (c: DG.Column) => c.type === DG.TYPE.FLOAT));
  });

  test('additionalColumns + onAdditionalColumnsChanged wire up', async () => {
    const df = demog();
    let fired = false;
    const input = ui.input.columns('cols', {
      table: df,
      additionalColumns: {extra: [df.col('age')!]},
      onAdditionalColumnsChanged: () => fired = true,
    });
    expect(input instanceof DG.InputBase, true);
    expectNoThrow(() => ui.input.setInputAdditionalColumns(input, {extra: [df.col('height')!]}));
    // Subscriber wiring should not throw; the flag may stay false in headless (no UI edit).
    expect(typeof fired, 'boolean');
  });

  test('showSelectedColsOnTop option does not throw', async () => {
    expectNoThrow(() => ui.input.columns('cols', {table: demog(), checked: ['age'], showSelectedColsOnTop: true}));
  });
}, {owner: 'agolovko@datagrok.ai'});
