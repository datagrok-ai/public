import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

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
          'state': {
            'width': 1000, 'height': 800,
            'element': {
              'id': '00000000-0000-0000-0000-000000000001',
              'type': viewerType,
              'look': look,
              'title': viewerType,
            },
          },
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

category('AI: Viewer layout backward compatibility', () => {
  // BoxPlotLook.fromMap migrates the legacy single-string `categoryColumnName`
  // (pre-trellis-categories) into the list-valued `categoryColumnNames`.
  // Condition (box_plot_look.dart:387):
  //   propMap['categoryColumnName'] != null && propMap['categoryColumnNames'] == null
  //
  // The final categoryColumnNames value on the recreated viewer is NOT directly
  // observable from JS: after fromMap migrates the legacy key, the layout-load
  // flow runs `onFrameAttached → look.auto()` and then `setLook → copyFrom →
  // notify()`, and box_plot_core.onLookChanged (lines ~428-436) rebuilds
  // categoryColumnNames from category1ColumnName on a generic notify. So we
  // pin what IS observable: the layout loads without throwing, the viewer is
  // created with the correct type, fromMap-deserialized props that auto()
  // does not overwrite (valueColumnName) survive, and the legacy key does not
  // surface in the look bag (it isn't a real property and gets consumed).
  test('Box plot: legacy categoryColumnName layout loads without throwing', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Box plot', {
        '#type': 'BoxPlotLook',
        'valueColumnName': 'age',
        'categoryColumnName': 'race',
      }, df.name);
      var threw = false;
      try {
        const layout = DG.ViewLayout.fromJson(json);
        tv.loadLayout(layout);
        await DG.delay(300);
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      const v = findViewer(tv, DG.VIEWER.BOX_PLOT);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(look['valueColumnName'], 'age');
      // Legacy key is consumed by the migration (propMap.remove) — it must
      // not surface in the look bag of the recreated viewer.
      expect(look['categoryColumnName'] === undefined, true);
      // The new key must exist as a List<String> regardless of which name(s)
      // ended up inside (the auto()/onLookChanged interaction described above).
      expect(Array.isArray(look['categoryColumnNames']), true);
    }
    finally {
      tv.close();
    }
  });

  // When BOTH legacy and new keys are present, the new key wins inside
  // fromMap (the legacy branch is guarded by `categoryColumnNames == null`).
  // Same observability caveat as above — we pin "no throw" + "no legacy key
  // surfaces".
  test('Box plot: layout with both legacy and new keys loads without throwing', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Box plot', {
        '#type': 'BoxPlotLook',
        'valueColumnName': 'age',
        'categoryColumnName': 'sex',          // legacy
        'categoryColumnNames': ['race'],      // new — wins inside fromMap
      }, df.name);
      var threw = false;
      try {
        const layout = DG.ViewLayout.fromJson(json);
        tv.loadLayout(layout);
        await DG.delay(300);
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      const v = findViewer(tv, DG.VIEWER.BOX_PLOT);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(look['categoryColumnName'] === undefined, true);
      expect(Array.isArray(look['categoryColumnNames']), true);
    }
    finally {
      tv.close();
    }
  });

  // PieChartLook.fromMap migrates two legacy keys:
  //   showInnerPercent -> showPercentage
  //   showInnerLabel   -> showLabel
  // (pie_chart_look.dart:155-163). Each is guarded by `legacy != null && new == null`.
  test('Pie chart: legacy showInnerPercent migrates to showPercentage', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Pie chart', {
        '#type': 'PieChartLook',
        'categoryColumnName': 'race',
        'showInnerPercent': true,
      }, df.name);
      const layout = DG.ViewLayout.fromJson(json);
      tv.loadLayout(layout);
      await DG.delay(300);
      const v = findViewer(tv, DG.VIEWER.PIE_CHART);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(look['showPercentage'], true);
      expect(look['showInnerPercent'] === undefined, true);
    }
    finally {
      tv.close();
    }
  });

  test('Pie chart: legacy showInnerLabel migrates to showLabel', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Pie chart', {
        '#type': 'PieChartLook',
        'categoryColumnName': 'race',
        'showInnerLabel': false,
      }, df.name);
      const layout = DG.ViewLayout.fromJson(json);
      tv.loadLayout(layout);
      await DG.delay(300);
      const v = findViewer(tv, DG.VIEWER.PIE_CHART);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(look['showLabel'], false);
      expect(look['showInnerLabel'] === undefined, true);
    }
    finally {
      tv.close();
    }
  });

  // Combined migration: both legacy keys migrate independently in one pass.
  test('Pie chart: both legacy keys migrate together', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Pie chart', {
        '#type': 'PieChartLook',
        'categoryColumnName': 'race',
        'showInnerPercent': true,
        'showInnerLabel': true,
      }, df.name);
      const layout = DG.ViewLayout.fromJson(json);
      tv.loadLayout(layout);
      await DG.delay(300);
      const v = findViewer(tv, DG.VIEWER.PIE_CHART);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(look['showPercentage'], true);
      expect(look['showLabel'], true);
      expect(look['showInnerPercent'] === undefined, true);
      expect(look['showInnerLabel'] === undefined, true);
    }
    finally {
      tv.close();
    }
  });

  // When the new key is already present, the legacy key is ignored (the
  // `propMap['showPercentage'] == null` guard prevents overwrite).
  test('Pie chart: new showPercentage wins when both are present', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Pie chart', {
        '#type': 'PieChartLook',
        'categoryColumnName': 'race',
        'showInnerPercent': true,    // legacy
        'showPercentage': false,     // new — should win
      }, df.name);
      const layout = DG.ViewLayout.fromJson(json);
      tv.loadLayout(layout);
      await DG.delay(300);
      const v = findViewer(tv, DG.VIEWER.PIE_CHART);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(look['showPercentage'], false);
    }
    finally {
      tv.close();
    }
  });

  // ScatterPlotLook.fromMap migrates labelFormColumnNames -> labelColumnNames
  // (scatterplot_look.dart:577-584).
  test('Scatter plot: legacy labelFormColumnNames migrates to labelColumnNames', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Scatter plot', {
        '#type': 'ScatterPlotLook',
        'xColumnName': 'age',
        'yColumnName': 'height',
        'labelFormColumnNames': ['race', 'sex'],
      }, df.name);
      const layout = DG.ViewLayout.fromJson(json);
      tv.loadLayout(layout);
      await DG.delay(300);
      const v = findViewer(tv, DG.VIEWER.SCATTER_PLOT);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(Array.isArray(look['labelColumnNames']), true);
      expect(look['labelColumnNames'].length, 2);
      expect(look['labelColumnNames'][0], 'race');
      expect(look['labelColumnNames'][1], 'sex');
      expect(look['labelFormColumnNames'] === undefined, true);
    }
    finally {
      tv.close();
    }
  });

  // Mirror of the box-plot/pie-chart "new wins" cases.
  test('Scatter plot: new labelColumnNames wins when both are present', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const json = buildLayoutJson('Scatter plot', {
        '#type': 'ScatterPlotLook',
        'xColumnName': 'age',
        'yColumnName': 'height',
        'labelFormColumnNames': ['sex'],      // legacy
        'labelColumnNames': ['race'],         // new — should win
      }, df.name);
      const layout = DG.ViewLayout.fromJson(json);
      tv.loadLayout(layout);
      await DG.delay(300);
      const v = findViewer(tv, DG.VIEWER.SCATTER_PLOT);
      expect(v != null, true);
      const look = v!.getOptions(true).look;
      expect(Array.isArray(look['labelColumnNames']), true);
      expect(look['labelColumnNames'].length, 1);
      expect(look['labelColumnNames'][0], 'race');
    }
    finally {
      tv.close();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
