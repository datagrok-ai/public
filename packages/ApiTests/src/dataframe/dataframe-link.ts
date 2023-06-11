// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';

category('DataFrame: Link', () => {
  let DF1: DG.DataFrame;
  let DF2: DG.DataFrame;

  before(async () => {
    DF1 = grok.data.demo.demog(10);
    DF2 = grok.data.demo.demog(20);
  });

  test('row to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [7, 4, 7, 1, 4, 6, 6, 6, 2, 7];
    grok.data.linkTables(df1, df2, ['sex', 'race'], ['sex', 'race'], [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER]);
    for (let i = 0; i < 10; i++) {
      df1.currentRowIdx = i;
      expect(df2.filter.trueCount, arr[i], 'row ' + i);
    }
  });

  test('row to row', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    grok.data.linkTables(df1, df2, ['sex', 'race'], ['sex', 'race'], [DG.SYNC_TYPE.CURRENT_ROW_TO_ROW]);
    for (let i = 0; i < 10; i++) {
      df1.currentRowIdx = i;
      expect(df2.currentRowIdx, i, 'row ' + i);
    }
  }, {skipReason: 'GROK-11670'});

  test('row to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [7, 4, 7, 1, 4, 4, 4, 2, 2, 7];
    grok.data.linkTables(df1, df2, ['sex', 'race', 'disease'],
      ['sex', 'race', 'disease'], [DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    for (let i = 0; i < 10; i++) {
      df1.currentRowIdx = i;
      expect(df2.selection.trueCount, arr[i], 'row ' + i);
    }
  });

  test('filter to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const race = df1.getCol('race');
    const arr = [];
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    for (const r of ['Caucasian', 'Asian', 'Black', 'Other']) {
      df1.filter.init((i) => race.get(i) === r);
      arr.push(df2.filter.trueCount);
    }
    console.log(arr);
  }, {skipReason: 'GROK-11670'});

  test('filter to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const race = df1.getCol('race');
    const arr = [];
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.FILTER_TO_SELECTION]);
    for (const r of ['Caucasian', 'Asian', 'Black', 'Other']) {
      df1.filter.init((i) => race.get(i) === r);
      arr.push(df2.selection.trueCount);
    }
    console.log(arr);
  }, {skipReason: 'GROK-11670'});

  test('mouse-over to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [4, 4, 7, 5, 4, 5, 5, 7, 4, 7];
    grok.data.linkTables(df1, df2, ['sex', 'site'],
      ['sex', 'site'], [DG.SYNC_TYPE.MOUSE_OVER_ROW_TO_FILTER]);
    for (let i = 0; i < 10; i++) {
      df1.mouseOverRowIdx = i;
      expect(df2.filter.trueCount, arr[i], 'row ' + i);
    }
  });

  test('mouse-over to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [11, 4, 11, 5, 4, 5, 5, 11, 11, 11];
    grok.data.linkTables(df1, df2, ['sex', 'disease'],
      ['sex', 'disease'], [DG.SYNC_TYPE.MOUSE_OVER_ROW_TO_SELECTION]);
    for (let i = 0; i < 10; i++) {
      df1.mouseOverRowIdx = i;
      expect(df2.selection.trueCount, arr[i], 'row ' + i);
    }
  });

  test('selection to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const race = df1.getCol('race');
    const arr = [];
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);
    for (const r of ['Caucasian', 'Asian', 'Black', 'Other']) {
      df1.selection.init((i) => race.get(i) === r);
      arr.push(df2.filter.trueCount);
    }
    console.log(arr);
  }, {skipReason: 'GROK-11670'});

  test('selection to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const race = df1.getCol('race');
    const arr = [];
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.SELECTION_TO_SELECTION]);
    for (const r of ['Caucasian', 'Asian', 'Black', 'Other']) {
      df1.selection.init((i) => race.get(i) === r);
      arr.push(df2.selection.trueCount);
    }
    console.log(arr);
  }, {skipReason: 'GROK-11670'});
});
