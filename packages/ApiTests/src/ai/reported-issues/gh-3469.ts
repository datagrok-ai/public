import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, look, until, wait, withTableView} from '../helpers';

// Regression coverage for #3469: Box plot color column auto-sync vs. manual override.
category('AI: gh-3469: Box plot color column auto-sync vs. manual override', () => {
  test('auto-sync: changing Category 1 propagates to markerColorColumnName when untouched', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await until(() => look(v)['category1ColumnName'] === 'race');
      expectLook(v, {category1ColumnName: 'race'});
      v.setOptions({category1ColumnName: 'sex'});
      await until(() => look(v)['markerColorColumnName'] === 'sex');
      expectLook(v, {category1ColumnName: 'sex', markerColorColumnName: 'sex'});
    });
  });

  test('manual override sticks: user-picked color column survives a Category 1 change', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await until(() => look(v)['markerColorColumnName'] != null);
      v.setOptions({markerColorColumnName: 'sex'});
      await until(() => look(v)['markerColorColumnName'] === 'sex');
      expectLook(v, {markerColorColumnName: 'sex'});
      v.setOptions({category1ColumnName: 'started'});
      await wait();
      expectLook(v, {category1ColumnName: 'started', markerColorColumnName: 'sex'});
    });
  });

  test('manual override sticks across Category 2 changes too', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race', category2: 'sex'});
      tv.addViewer(v);
      await until(() => look(v)['markerColorColumnName'] != null);
      v.setOptions({markerColorColumnName: 'started'});
      await until(() => look(v)['markerColorColumnName'] === 'started');
      expectLook(v, {markerColorColumnName: 'started'});
      v.setOptions({category2ColumnName: 'disease'});
      await wait();
      expectLook(v, {category2ColumnName: 'disease', markerColorColumnName: 'started'});
    });
  });

  test('global escape hatch: allowColorSynchronization=false freezes the color column on category change', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.boxPlot(tv.dataFrame, {value: 'age', category1: 'race'});
      tv.addViewer(v);
      await until(() => look(v)['markerColorColumnName'] === 'race');
      expectLook(v, {markerColorColumnName: 'race'});
      v.setOptions({allowColorSynchronization: false});
      await wait(100);
      v.setOptions({category1ColumnName: 'sex'});
      await wait();
      expectLook(v, {category1ColumnName: 'sex', markerColorColumnName: 'race'});
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
