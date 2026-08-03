import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, until, withTableView} from '../helpers';

// Legacy-property migration in viewer Look.fromMap, exercised via layout load
// (not setOptions, which no-ops on legacy keys).

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
              'type': viewerType, 'look': look, 'title': viewerType,
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
  for (const v of tv.viewers) {
    if (v.type === type)
      return v;
  }
  return null;
}

async function loadLayoutGetLook(tv: DG.TableView, viewerType: string, look: {[k: string]: any}, dfName: string,
  expectedType: string, noThrow: boolean = false): Promise<{[k: string]: any}> {
  const json = buildLayoutJson(viewerType, look, dfName);
  if (noThrow)
    expectNoThrow(() => tv.loadLayout(DG.ViewLayout.fromJson(json)));
  else
    tv.loadLayout(DG.ViewLayout.fromJson(json));
  await until(() => findViewer(tv, expectedType) != null);
  const v = findViewer(tv, expectedType);
  expect(v != null, true);
  return v!.getOptions(true).look;
}

category('AI: Viewer layout backward compatibility', () => {
  // Migrated categoryColumnNames isn't directly observable (auto/notify rebuilds it), so assert
  // what is: loads without throwing, survives valueColumnName, drops the legacy key.
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

  // When both legacy and new keys are present, the new key wins.
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
