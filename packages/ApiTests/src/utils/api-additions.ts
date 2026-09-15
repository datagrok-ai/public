import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Additive surface from phase 3 of the 2026-09 JS API audit (core/docs/reviews/js-api-audit/plan.md).
category('JS API: additions', () => {
  const table = () => DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.TYPE.INT, 'x', [1, 2, 3]),
    DG.Column.fromList(DG.TYPE.STRING, 'name', ['a', 'b', 'c']),
  ]);

  test('typed bags: tags hold strings, temp holds anything, map methods work', async () => {
    const df = table();
    df.tags['k'] = 'v';
    df.temp['obj'] = {n: 1};
    expect(df.tags['k'], 'v');
    expect(df.tags.has('k'), true);
    expect(df.temp['obj'].n, 1);
    expect([...df.tags.keys()].includes('k'), true);
    delete df.tags['k'];
    expect(df.tags.has('k'), false);
  });

  test('DataFrameMetaHelper.setGroups(null) clears the tag instead of throwing', async () => {
    const df = table();
    df.meta.setGroups({groups: [{name: 'g', columns: ['x']}]} as any);
    expect(df.tags.has('.columnGroups'), true);
    df.meta.setGroups(null);
    expect(df.tags.has('.columnGroups'), false);
  });

  test('Grid.cell accepts (row, column) like DataFrame.cell', async () => {
    const grid = DG.Viewer.grid(table());
    expect(grid.cell(2, 'x').cell.value, grid.cell('x', 2).cell.value);
    expect(grid.cell(2, 'x').cell.value, 3);
  });

  test('toCsvAsync equals toCsvEx', async () => {
    const df = table();
    expect(await df.toCsvAsync(), await df.toCsvEx());
  });

  test('Column factories return typed subclasses', async () => {
    expect(DG.Column.dateTime('d', 2) instanceof DG.DateTimeColumn, true);
    const s: DG.Column<string> = DG.Column.fromStrings('s', ['a']);
    expect(s.get(0), 'a');
    const n: DG.Column<number> = DG.Column.fromList(DG.TYPE.INT, 'n', [7]);
    expect(n.get(0), 7);
  });

  test('Viewer factories return the viewer classes, with *Viewer aliases', async () => {
    const df = grok.data.demo.demog(50);
    expect(DG.Viewer.barChart(df) instanceof DG.BarChartViewer, true);
    expect(DG.Viewer.pieChart(df) instanceof DG.PieChartViewer, true);
    expect(DG.Viewer.boxPlot(df) instanceof DG.BoxPlotViewer, true);
    expect(DG.Viewer.pcPlot(df) instanceof DG.PcPlotViewer, true);
    expect(DG.Viewer.trellisPlot(df) instanceof DG.TrellisPlotViewer, true);
    expect(DG.Viewer.pivotTable(df) instanceof DG.PivotViewer, true);
  });

  test('shell.currentTable / currentView / currentTableView / currentObject', async () => {
    const df = table();
    df.name = 'additions';
    const view = grok.shell.addTableView(df);
    try {
      expect(grok.shell.currentTable?.name, 'additions');
      expect(grok.shell.currentTableView?.dataFrame.name, 'additions');
      expect(grok.shell.currentView?.name, view.name);
      grok.shell.currentObject = df;
      expect(grok.shell.currentObject?.name, 'additions');
    }
    finally {
      view.close();
      grok.shell.closeTable(df);
    }
  });

  test('options-object overloads: clone, join, linkTables, addPane, rangeSlider', async () => {
    const df = table();
    const sub = df.clone({columns: ['x'], rows: DG.BitSet.create(3, (i) => i > 0)});
    expect(sub.columns.length, 1);
    expect(sub.rowCount, 2);
    const other = DG.DataFrame.fromColumns([DG.Column.fromList(DG.TYPE.INT, 'x', [1, 2]), DG.Column.fromList(DG.TYPE.STRING, 'v', ['p', 'q'])]);
    const joined = df.join(other, {keys: ['x'], columns: ['name'], columns2: ['v'], type: DG.JOIN_TYPE.LEFT});
    expect(joined.rowCount, 3);
    expect(joined.columns.names().includes('v'), true);
    grok.data.linkTables(df, other, ['x'], ['x'], [DG.SYNC_TYPE.SELECTION_TO_SELECTION], {initialSync: true});
    const acc = ui.accordion();
    const pane = acc.addPane('p', () => ui.div(), {expanded: true});
    expect(pane.expanded, true);
    const slider = ui.rangeSlider({minRange: 0, maxRange: 100, min: 10, max: 20});
    expect(slider.min, 10);
    expect(slider.max, 20);
  });

  test('ui.form layout and ui.icon dispatch', async () => {
    expect(ui.form([ui.input.string('a')], {layout: 'narrow'}).classList.contains('ui-form-condensed'), true);
    expect(ui.form([ui.input.string('a')], {layout: 'wide'}).classList.contains('ui-form-wide'), true);
    expect(ui.icon('plus').classList.contains('fa-plus'), true);
    expect(ui.icon('ai.svg').classList.contains('svg-icon'), true);
    expect(ui.icon('/images/x.png').classList.contains('image-icon'), true);
    let clicked = false;
    const i = ui.icon('plus', {onClick: () => clicked = true, tooltip: 'add'});
    i.click();
    expect(clicked, true);
  });

  test('Viewer.fromType returns the class its type maps to', async () => {
    const df = grok.data.demo.demog(20);
    const sp: DG.ScatterPlotViewer = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, df);
    expect(sp instanceof DG.ScatterPlotViewer, true);
    const byLiteral: DG.BarChartViewer = DG.Viewer.fromType('Bar chart', df);
    expect(byLiteral instanceof DG.BarChartViewer, true);
    const heat: DG.Grid = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, df);
    expect(heat instanceof DG.Grid, true);
    expect(DG.Viewer.fromType(DG.VIEWER.FILTERS, df) instanceof DG.FilterGroup, true);
    expect(DG.Viewer.fromType(DG.VIEWER.PIVOT_TABLE, df) instanceof DG.PivotViewer, true);
    const runtimeType: string = DG.VIEWER.PIE_CHART;
    const generic: DG.Viewer = DG.Viewer.fromType(runtimeType, df);
    expect(generic instanceof DG.PieChartViewer, true);
    for (const v of [sp, byLiteral, heat, generic])
      v.detach();
  });

  test('HttpDataSource verbs are immutable', async () => {
    const users = grok.dapi.users;
    const admins = users.filter('login = "admin"');
    expect(admins === users, false);
    const all = await users.list();
    expect(all.length > 1, true, 'the stand has more than one user');
    expect((await admins.list()).length, 1);
    expect((await admins.count()), 1);
    await users.first();
    expect((await users.list()).length, all.length, 'first() must not shrink later lists');
    expect((await users.list({pageSize: 1})).length, 1);
    expect((await users.list()).length, all.length, 'list options must not stick');
    const paged = users.by(1);
    expect((await paged.list()).length, 1);
    expect((await paged.nextPage().list()).length, 1);
    expect((await paged.list()).length, 1, 'nextPage() must not advance its source');
    expect((await users.order('login', true).first()).login >= (await users.order('login').first()).login, true);
    expect((await users.list()).length, all.length);
  });
});
