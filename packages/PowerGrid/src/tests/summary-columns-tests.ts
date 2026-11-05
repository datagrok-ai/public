import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {waitUntil, category, expect, expectArray, expectObject, test} from '@datagrok-libraries/test/src/test';
import {TagsCellRenderer} from '../cell-types/tags-cell-renderer';
import {BarCellRenderer} from '../cell-types/bar-cell-renderer';
import {layoutRoundTrip, projectRoundTrip} from './utils';

/** The renderer items of Add > Summary Columns, each with settings its settings panel can make. */
const SUMMARY_COLUMNS: {name: string, cellType: string, settings: {[key: string]: any}}[] = [
  {name: 'Sparklines', cellType: 'sparkline',
    settings: {columnNames: ['age', 'weight'], normalization: 'Row', zeroScale: true}},
  {name: 'Bar Chart', cellType: 'barchart', settings: {columnNames: ['height', 'weight'], colorCode: 'Off'}},
  {name: 'Pie Chart', cellType: 'piechart', settings: {columnNames: ['age', 'height', 'weight'], style: 'Angle'}},
  {name: 'VlaaiVis', cellType: 'vlaaivis', settings: {columnNames: ['age', 'height'], normalization: 'Global'}},
  {name: 'Radar', cellType: 'radar',
    settings: {columnNames: ['age', 'height', 'weight'], normalization: 'Global'}},
  {name: 'Smart Form', cellType: 'smartform',
    settings: {columnNames: ['sex', 'race', 'age'], showColumnNames: 'Always'}},
  {name: 'Tags', cellType: 'Tags', settings: {columnNames: ['elderly']}},
  {name: 'Confidence Interval', cellType: 'ConfidenceInterval',
    settings: {type: 'Value + Margin', estimateColumn: 'height', lowerColumn: 'age', centerMark: 'Diamond',
      showSerifs: false}},
];

const CHIP_COLOR = DG.Color.getCategoricalColor(0);

function demog(): DG.DataFrame {
  const df = grok.data.demo.demog(100);
  df.columns.addNewBool('elderly').init((i) => df.get('age', i) > 60);
  return df;
}

function summaryColumnsView(): DG.TableView {
  const view = grok.shell.addTableView(demog());
  for (const c of SUMMARY_COLUMNS) {
    const gc = view.grid.columns.add({gridColumnName: c.name, cellType: c.cellType});
    gc.settings = {[c.cellType]: structuredClone(c.settings)};
  }
  return view;
}

async function expectSummaryColumns(view: DG.TableView): Promise<void> {
  const grid = view.grid;
  const names = () => Array.from({length: grid.columns.length - 1}, (_, i) => grid.columns.byIndex(i + 1)!.name);
  await waitUntil(() => grid.columns.byName(SUMMARY_COLUMNS[SUMMARY_COLUMNS.length - 1].name) != null,
    'the summary columns did not come back', 5000);
  expectArray(names(), [...view.dataFrame.columns.names(), ...SUMMARY_COLUMNS.map((c) => c.name)]);
  for (const c of SUMMARY_COLUMNS) {
    const gc = grid.columns.byName(c.name)!;
    expect(gc.cellType, c.cellType.toLowerCase(), c.name);
    await waitUntil(() => gc.renderer?.cellType === c.cellType,
      `${c.name} is not drawn by the ${c.cellType} renderer`, 5000);
    expectObject(gc.settings, {[c.cellType]: c.settings});
  }
}

function flagsTable(): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('name', ['a', 'b', 'c', 'd']),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag', [false, true, false, true]),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'value', [1, 2, 3, 4]),
  ]);
}

/** Pixels of [g] in the device-pixel rectangle that pass [accept]. */
function countPixels(g: CanvasRenderingContext2D, r: DG.Rect,
  accept: (rgba: Uint8ClampedArray, i: number) => boolean): number {
  const data = g.getImageData(Math.round(r.x), Math.round(r.y), Math.round(r.width), Math.round(r.height)).data;
  let n = 0;
  for (let i = 0; i < data.length; i += 4)
    n += accept(data, i) ? 1 : 0;
  return n;
}

const isChip = (d: Uint8ClampedArray, i: number) => Math.abs(d[i] - DG.Color.r(CHIP_COLOR)) < 3 &&
  Math.abs(d[i + 1] - DG.Color.g(CHIP_COLOR)) < 3 && Math.abs(d[i + 2] - DG.Color.b(CHIP_COLOR)) < 3;

const isPainted = (d: Uint8ClampedArray, i: number) => d[i + 3] > 0;

/** Chip-colored pixels the grid shows in a cell. */
function chipPixels(grid: DG.Grid, column: string, row: number): number {
  const b = grid.cell(column, row).bounds;
  const dpr = window.devicePixelRatio;
  return countPixels(grid.canvas.getContext('2d')!,
    new DG.Rect(b.x * dpr, b.y * dpr, b.width * dpr, b.height * dpr), isChip);
}

category('Summary columns', () => {
  test('Layout', async () => {
    await layoutRoundTrip(summaryColumnsView(), demog(), expectSummaryColumns);
  });

  test('Project', async () => {
    await projectRoundTrip(summaryColumnsView(), expectSummaryColumns);
  });

  test('Tags marks the flagged rows and leaves the other cells alone (GROK-20888)', async () => {
    const df = flagsTable();
    const grid = grok.shell.addTableView(df).grid;
    const tags = grid.columns.add({gridColumnName: 'Tags', cellType: 'Tags'});
    await waitUntil(() => tags.renderer?.cellType === 'Tags', 'the Tags renderer did not load', 5000);
    grid.invalidate();
    await waitUntil(() => chipPixels(grid, 'Tags', 1) > 0, 'the flagged row shows no chip', 5000);
    const bounds = grid.cell('Tags', 1).bounds;
    const area = bounds.width * bounds.height * window.devicePixelRatio ** 2;
    for (let row = 0; row < df.rowCount; row++) {
      for (const column of ['name', 'flag', 'value', 'Tags']) {
        const n = chipPixels(grid, column, row);
        if (column === 'Tags' && df.get('flag', row))
          expect(n > 0 && n < area / 2, true, `the chip of row ${row} covers ${n} of ${area} pixels`);
        else
          expect(n, 0, `chip-colored pixels in ${column} of row ${row}`);
      }
    }
  });

  test('A removed source column brought back rejoins its Tags column (#4378)', async () => {
    const df = flagsTable();
    const grid = grok.shell.addTableView(df).grid;
    const tags = grid.columns.add({gridColumnName: 'Tags', cellType: 'Tags'});
    tags.settings = {Tags: {columnNames: ['flag']}};
    await waitUntil(() => tags.renderer?.cellType === 'Tags', 'the Tags renderer did not load', 5000);
    const flag = df.col('flag')!.clone();
    df.columns.remove('flag');
    expectArray(tags.settings['Tags'].columnNames, []);
    df.columns.insert(flag, 1);
    expectArray(tags.settings['Tags'].columnNames, ['flag']);
    grid.invalidate();
    await waitUntil(() => chipPixels(grid, 'Tags', 1) > 0, 'the flagged row shows no chip after the column came back', 5000);
  });

  test('Tags and Bar fill only their own shapes (GROK-20888)', async () => {
    const grid = grok.shell.addTableView(flagsTable()).grid;
    grid.columns.add({gridColumnName: 'Tags', cellType: 'Tags'});
    const cases: [DG.GridCellRenderer, string, number, boolean][] = [
      [new TagsCellRenderer(), 'Tags', 1, true],
      [new TagsCellRenderer(), 'Tags', 0, false],
      [new BarCellRenderer(), 'value', 1, true],
    ];
    for (const [renderer, column, row, paints] of cases) {
      const g = ui.canvas(200, 60).getContext('2d')!;
      // a path an earlier drawing left open: filling it would paint the whole canvas
      g.beginPath();
      g.rect(0, 0, 200, 60);
      const cell = grid.cell(column, row);
      renderer.render(g, 60, 16, 80, 28, cell, cell.style);
      const inCell = countPixels(g, new DG.Rect(60, 16, 80, 28), isPainted);
      const outside = countPixels(g, new DG.Rect(0, 0, 200, 60), isPainted) - inCell;
      const what = `${renderer.cellType} on ${column} of row ${row}`;
      expect(outside, 0, `${what}: pixels painted outside the cell`);
      expect(inCell > 0 && inCell < 80 * 28 / 2, paints, `${what}: ${inCell} pixels painted in the cell`);
    }
  });

  test('Tags column is as wide as its widest cell', async () => {
    const values = ['span < 2,N ≤ 3,%extrap > 20%', 'span < 2', ''];
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('violated', values)]);
    const col = df.col('violated')!;
    col.setTag(DG.Tags.MultiValueSeparator, ',');
    col.setTag(DG.TAGS.CELL_RENDERER, 'Tags');
    const grid = grok.shell.addTableView(df).grid;
    const gc = grid.col('violated')!;
    await waitUntil(() => gc.renderer?.cellType === 'Tags', 'the Tags renderer did not load', 5000);
    await waitUntil(() => !grid.isRenderPending, 'the grid did not settle', 5000);
    const g = grid.canvas.getContext('2d')!;
    g.font = gc.contentCellStyle.font ?? grid.props.defaultCellFont;
    const layout = new TagsCellRenderer().calculateLayout(g, gc.width, grid.props.rowHeight, values[0].split(','));
    expectArray(layout.map((item) => item.mode), ['full', 'full', 'full']);
    expect(gc.width < 250, true, `width ${gc.width}`);
  });

  test('Color fills each legend item with its own category (#4233)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('color', ['red', 'lime', 'blue', 'red']),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'value', [1, 2, 3, 4]),
    ]);
    df.col('color')!.setTag(DG.TAGS.CELL_RENDERER, 'Color');
    const view = grok.shell.addTableView(df);
    await waitUntil(() => view.grid.col('color')?.renderer?.cellType === 'Color',
      'the Color renderer did not load', 5000);
    const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT, {x: 'value', y: 'value', colorColumnName: 'color'});
    const centers = () => Array.from(sp.root.querySelectorAll<HTMLCanvasElement>('canvas.d4-legend-value'), (c) =>
      Array.from(c.getContext('2d')!.getImageData(c.width >> 1, c.height >> 1, 1, 1).data.slice(0, 3)).join(','));
    await waitUntil(() => centers().length === 3, 'the legend shows no Color items', 5000);
    expectArray(centers().sort(), ['0,0,255', '0,255,0', '255,0,0']);
  });
});
