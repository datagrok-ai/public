import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, wait, withTableView} from '../helpers';

// Layout backward-compatibility coverage for viewer Look classes that override
// `fromMap` to migrate legacy property names when loading an old layout JSON.
// Source (Dart):
//   core/client/d4/lib/src/viewers/box_plot/box_plot_look.dart    (categoryColumnName -> categoryColumnNames)
//   core/client/d4/lib/src/viewers/pie_chart/pie_chart_look.dart  (showInnerPercent -> showPercentage, showInnerLabel -> showLabel)
//   core/client/d4/lib/src/viewers/scatterplot/scatterplot_look.dart (labelFormColumnNames -> labelColumnNames)
//
// The migration path is only triggered when the layout JSON is fed through
// LookAndFeel.fromMap, which is invoked by ViewerDescriptor.init(json: ...)
// during TableView.loadState — i.e. via the layout-load code path, NOT
// setOptions. setOptions iterates options.keys and calls setProperty per
// key, so legacy keys silently no-op when fed through it. Tests therefore
// construct a ViewLayout from a JSON document with legacy keys, apply it
// via tv.loadLayout, then read v.getOptions(true).look on the recreated
// viewer to assert the migrated state.

function buildLayoutJson(viewerType: string, look: {[k: string]: any}, tableName: string): string {
  return JSON.stringify({
    '#type': 'ViewLayout',
    'viewStateMap': {
      'containerType': 'horizontal',
      'state': {'width': 1000, 'height': 800},
      'children': [{
        'containerType': 'fill',
        'state': {'width': 1000, 'height': 800, 'documentManager': true},
        'children': [{
          'containerType': 'panel',
          'state': {'width': 1000, 'height': 800,
            'element': {'id': '00000000-0000-0000-0000-000000000001', 'type': viewerType, 'look': look, 'title': viewerType}},
          'children': [],
        }],
      }],
      'floating': [],
      'tableName': tableName,
    },
    'type': 'TableView',
  });
}

function findViewer(tv: DG.TableView, type: string): DG.Viewer | null {
  for (var v of tv.viewers)
    if (v.type === type) return v;
  return null;
}

async function loadLayoutGetLook(tv: DG.TableView, viewerType: string, look: {[k: string]: any}, dfName: string,
  expectedType: string, noThrow: boolean = false): Promise<{[k: string]: any}> {
  const json = buildLayoutJson(viewerType, look, dfName);
  if (noThrow)
    expectNoThrow(() => tv.loadLayout(DG.ViewLayout.fromJson(json)));
  else
    tv.loadLayout(DG.ViewLayout.fromJson(json));
  await wait();
  const v = findViewer(tv, expectedType);
  expect(v != null, true);
  return v!.getOptions(true).look;
}

category('AI: Viewer layout backward compatibility', () => {
  // BoxPlotLook.fromMap migrates the legacy single-string `categoryColumnName`
  // (pre-trellis-categories) into the list-valued `categoryColumnNames`.
  //
  // The final categoryColumnNames value on the recreated viewer is NOT directly
  // observable from JS: after fromMap migrates the legacy key, the layout-load
  // flow runs `onFrameAttached → look.auto()` and then `setLook → copyFrom →
  // notify()`, and box_plot_core.onLookChanged (lines ~428-436) rebuilds
  // categoryColumnNames from category1ColumnName on a generic notify. So we
  // pin what IS observable: the layout loads without throwing, the viewer is
  // created with the correct type, fromMap-deserialized props that auto()
  // does not overwrite (valueColumnName) survive, and the legacy key does not
  // surface in the look bag.
  test('Box plot: legacy categoryColumnName layout loads without throwing', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Box plot',
        {'#type': 'BoxPlotLook', 'valueColumnName': 'age', 'categoryColumnName': 'race'},
        df.name, DG.VIEWER.BOX_PLOT, true);
      expect(look['valueColumnName'], 'age');
      expect(look['categoryColumnName'] === undefined, true);
      expect(Array.isArray(look['categoryColumnNames']), true);
    });
  });

  // When BOTH legacy and new keys are present, the new key wins inside fromMap
  // (the legacy branch is guarded by `categoryColumnNames == null`).
  test('Box plot: layout with both legacy and new keys loads without throwing', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Box plot',
        {'#type': 'BoxPlotLook', 'valueColumnName': 'age',
          'categoryColumnName': 'sex', 'categoryColumnNames': ['race']},
        df.name, DG.VIEWER.BOX_PLOT, true);
      expect(look['categoryColumnName'] === undefined, true);
      expect(Array.isArray(look['categoryColumnNames']), true);
    });
  });

  // PieChartLook.fromMap migrates:
  //   showInnerPercent -> showPercentage
  //   showInnerLabel   -> showLabel
  test('Pie chart: legacy showInnerPercent migrates to showPercentage', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Pie chart',
        {'#type': 'PieChartLook', 'categoryColumnName': 'race', 'showInnerPercent': true},
        df.name, DG.VIEWER.PIE_CHART);
      expect(look['showPercentage'], true);
      expect(look['showInnerPercent'] === undefined, true);
    });
  });

  test('Pie chart: legacy showInnerLabel migrates to showLabel', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Pie chart',
        {'#type': 'PieChartLook', 'categoryColumnName': 'race', 'showInnerLabel': false},
        df.name, DG.VIEWER.PIE_CHART);
      expect(look['showLabel'], false);
      expect(look['showInnerLabel'] === undefined, true);
    });
  });

  test('Pie chart: both legacy keys migrate together', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Pie chart',
        {'#type': 'PieChartLook', 'categoryColumnName': 'race',
          'showInnerPercent': true, 'showInnerLabel': true},
        df.name, DG.VIEWER.PIE_CHART);
      expect(look['showPercentage'], true);
      expect(look['showLabel'], true);
      expect(look['showInnerPercent'] === undefined, true);
      expect(look['showInnerLabel'] === undefined, true);
    });
  });

  // When the new key is already present, the legacy key is ignored.
  test('Pie chart: new showPercentage wins when both are present', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Pie chart',
        {'#type': 'PieChartLook', 'categoryColumnName': 'race',
          'showInnerPercent': true, 'showPercentage': false},
        df.name, DG.VIEWER.PIE_CHART);
      expect(look['showPercentage'], false);
    });
  });

  // ScatterPlotLook.fromMap migrates labelFormColumnNames -> labelColumnNames.
  test('Scatter plot: legacy labelFormColumnNames migrates to labelColumnNames', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Scatter plot',
        {'#type': 'ScatterPlotLook', 'xColumnName': 'age', 'yColumnName': 'height',
          'labelFormColumnNames': ['race', 'sex']},
        df.name, DG.VIEWER.SCATTER_PLOT);
      expect(Array.isArray(look['labelColumnNames']), true);
      expect(look['labelColumnNames'].length, 2);
      expect(look['labelColumnNames'][0], 'race');
      expect(look['labelColumnNames'][1], 'sex');
      expect(look['labelFormColumnNames'] === undefined, true);
    });
  });

  test('Scatter plot: new labelColumnNames wins when both are present', async () => {
    const df = demog();
    await withTableView(df, async (tv) => {
      const look = await loadLayoutGetLook(tv, 'Scatter plot',
        {'#type': 'ScatterPlotLook', 'xColumnName': 'age', 'yColumnName': 'height',
          'labelFormColumnNames': ['sex'], 'labelColumnNames': ['race']},
        df.name, DG.VIEWER.SCATTER_PLOT);
      expect(Array.isArray(look['labelColumnNames']), true);
      expect(look['labelColumnNames'].length, 1);
      expect(look['labelColumnNames'][0], 'race');
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
