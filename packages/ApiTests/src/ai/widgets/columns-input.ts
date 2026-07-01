import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog} from '../helpers';

// ColumnsInput (rich column multi-picker) — ui.input.columns + setColumnsInputTable,
// setInputAdditionalColumns and the available/checked/showSelectedColsOnTop creation options.
function names(input: DG.InputBase<DG.Column[]>): string[] {
  return input.value.map((c) => c.name).sort();
}

category('AI: Widgets: ColumnsInput', () => {
  test('value get/set round-trips the column selection', async () => {
    const t = demog();
    const input = ui.input.columns('cols', {table: t});
    expect(input instanceof DG.InputBase, true);
    input.value = [t.col('age')!, t.col('height')!];
    expect(input.value.length, 2);
    expect(names(input).join(','), 'age,height');
  });

  test('value setter round-trips columns from the available subset', async () => {
    const t = demog();
    const input = ui.input.columns('cols', {table: t, available: ['age', 'height', 'weight']});
    input.value = [t.col('age')!, t.col('weight')!];
    expect(input.value.length, 2);
    expect(names(input).join(','), 'age,weight');
  });

  test('setColumnsInputTable re-binds the table; value round-trips columns from the new table', async () => {
    const input = ui.input.columns('cols', {table: demog(50)});
    const t2 = demog(100);
    ui.input.setColumnsInputTable(input, t2, (c: DG.Column) => c.type === DG.TYPE.FLOAT);
    input.value = [t2.col('height')!, t2.col('weight')!];
    expect(input.value.length, 2);
    expect(names(input).join(','), 'height,weight');
  });

  test('additionalColumns wire up; value still round-trips after re-binding extras', async () => {
    const df = demog();
    const input = ui.input.columns('cols', {
      table: df,
      additionalColumns: {extra: [df.col('age')!]},
      onAdditionalColumnsChanged: () => {},
    });
    ui.input.setInputAdditionalColumns(input, {extra: [df.col('height')!]});
    input.value = [df.col('age')!, df.col('sex')!];
    expect(input.value.length, 2);
    expect(names(input).join(','), 'age,sex');
  });

  test('showSelectedColsOnTop input round-trips the value', async () => {
    const t = demog();
    const input = ui.input.columns('cols', {table: t, checked: ['age'], showSelectedColsOnTop: true});
    input.value = [t.col('age')!, t.col('height')!];
    expect(input.value.length, 2);
    expect(names(input).join(','), 'age,height');
  });
}, {owner: 'agolovko@datagrok.ai'});
