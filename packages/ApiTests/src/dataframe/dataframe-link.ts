// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';

category('DataFrame: Methods: Link', () => {
  let DF1: DG.DataFrame;
  let DF2: DG.DataFrame;

  before(async () => {
    DF1 = grok.data.demo.demog(10);
    DF2 = grok.data.demo.demog(20);
  });

  test('row to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [14, 4, 14, 1, 4, 14, 14, 14, 1, 14];
    grok.data.linkTables(df1, df2, ['sex', 'race'], ['sex', 'race'], [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER]);
    for (let i = 0; i < 10; i++) {
      df1.currentRowIdx = i;
      expect(df2.filter.trueCount, arr[i], 'row ' + (i + 1));
    }
  });

  test('row to row', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    grok.data.linkTables(df1, df2, ['subj'], ['subj'], [DG.SYNC_TYPE.CURRENT_ROW_TO_ROW]);
    for (let i = 0; i < 10; i++) {
      df1.currentRowIdx = i;
      expect(df2.currentRowIdx, i, 'row ' + (i + 1));
    }
  });

  test('row to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [4, 4, 2, 1, 4, 8, 8, 8, 1, 2];
    grok.data.linkTables(df1, df2, ['sex', 'race', 'disease'],
      ['sex', 'race', 'disease'], [DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION]);
    for (let i = 0; i < 10; i++) {
      df1.currentRowIdx = i;
      expect(df2.selection.trueCount, arr[i], 'row ' + (i + 1));
    }
  });

  test('filter to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [14, 4, 1, 1];
    const races = ['Caucasian', 'Asian', 'Black', 'Other'];
    const race = df1.getCol('race');
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.FILTER_TO_FILTER]);
    for (let i = 0; i < 4; i++) {
      const r = races[i];
      df1.filter.init((j) => race.get(j) === r);
      expect(df2.filter.trueCount, arr[i], r);
    }
  });

  test('filter to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [14, 4, 1, 1];
    const races = ['Caucasian', 'Asian', 'Black', 'Other'];
    const race = df1.getCol('race');
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.FILTER_TO_SELECTION]);
    for (let i = 0; i < 4; i++) {
      const r = races[i];
      df1.filter.init((j) => race.get(j) === r);
      expect(df2.selection.trueCount, arr[i], r);
    }
  });

  test('mouse-over to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [8, 3, 8, 7, 3, 7, 7, 7, 1, 8];
    grok.data.linkTables(df1, df2, ['sex', 'site'],
      ['sex', 'site'], [DG.SYNC_TYPE.MOUSE_OVER_ROW_TO_FILTER]);
    for (let i = 0; i < 10; i++) {
      df1.mouseOverRowIdx = i;
      expect(df2.filter.trueCount, arr[i], 'row ' + (i + 1));
    }
  });

  test('mouse-over to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [5, 4, 2, 9, 4, 9, 9, 9, 5, 2];
    grok.data.linkTables(df1, df2, ['sex', 'disease'],
      ['sex', 'disease'], [DG.SYNC_TYPE.MOUSE_OVER_ROW_TO_SELECTION]);
    for (let i = 0; i < 10; i++) {
      df1.mouseOverRowIdx = i;
      expect(df2.selection.trueCount, arr[i], 'row ' + (i + 1));
    }
  });

  test('selection to filter', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [14, 4, 1, 1];
    const races = ['Caucasian', 'Asian', 'Black', 'Other'];
    const race = df1.getCol('race');
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);
    for (let i = 0; i < 4; i++) {
      const r = races[i];
      df1.selection.init((j) => race.get(j) === r);
      expect(df2.filter.trueCount, arr[i], r);
    }
  });

  test('selection to selection', async () => {
    const df1 = DF1.clone();
    const df2 = DF2.clone();
    const arr = [14, 4, 1, 1];
    const races = ['Caucasian', 'Asian', 'Black', 'Other'];
    const race = df1.getCol('race');
    grok.data.linkTables(df1, df2, ['race'], ['race'], [DG.SYNC_TYPE.SELECTION_TO_SELECTION]);
    for (let i = 0; i < 4; i++) {
      const r = races[i];
      df1.selection.init((j) => race.get(j) === r);
      expect(df2.selection.trueCount, arr[i], r);
    }
  });
}, {owner: 'drizhinashvili@datagrok.ai'});
