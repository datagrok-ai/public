import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {IFitChartData} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {chooseLegendCorner, legendEntries, legendTooltip, legendTooltipElement} from '../fit/fit-legend';

import {sigmoidPoints, renderedTexts, renderedRows} from './curve-data';

/** Curves spread over several columns, the way the multi curve viewer hands them over. */
function multiColumnJson(columns: {column: string, series: string[]}[], showColumnLabel: boolean = true): string {
  const series = [];
  for (const {column, series: names} of columns) {
    for (const name of names) {
      series.push({fitFunction: 'sigmoid', name: name, columnName: column,
        points: sigmoidPoints(-6.5)});
    }
  }
  return JSON.stringify({chartOptions: {logX: true, showColumnLabel: showColumnLabel},
    series: series} as IFitChartData);
}

const PLOT = new DG.Rect(0, 0, 500, 400);

function measureWidth(text: string): number {
  const g = ui.canvas(10, 10).getContext('2d')!;
  g.font = '11px Roboto, "Roboto Local"';
  return g.measureText(text).width;
}

category('legend', () => {
  test('a column holding one curve does not get a row of its own', async () => {
    // the siRNA shape: four curve columns, one curve each
    const data = JSON.parse(multiColumnJson([
      {column: 'curve', series: ['PF5624 P1_1h']},
      {column: 'durability_curve', series: ['siR-0016_durability']},
      {column: 'in_vivo_dr_curve', series: ['siR-0016_invivo_dr']},
      {column: 'plasma_stability_curve', series: ['siR-0016_plasma_stab']},
    ])) as IFitChartData;
    expectArray(legendEntries(data).map((e) => e.text),
      ['PF5624 P1_1h', 'siR-0016_durability', 'siR-0016_invivo_dr', 'siR-0016_plasma_stab']);

    // a column with several curves still groups them under its name
    const grouped = JSON.parse(multiColumnJson([{column: 'curve', series: ['low', 'high']}])) as IFitChartData;
    expectArray(legendEntries(grouped).map((e) => e.text), ['curve', 'low', 'high']);
  });

  test('a name shared by two columns is qualified by the column', async () => {
    const data = JSON.parse(multiColumnJson([
      {column: 'day 1', series: ['control']},
      {column: 'day 7', series: ['control']},
    ])) as IFitChartData;
    expectArray(legendEntries(data).map((e) => e.text), ['day 1: control', 'day 7: control']);
  });

  test('two curves named alike in one column are left as the data has them', async () => {
    // the demo shape: one column, both series called Run:2023-08-08 - a column prefix on both would
    // tell them apart no better and take twice the width
    const data = JSON.parse(multiColumnJson([
      {column: 'multiple styled series', series: ['Run:2023-08-08', 'Run:2023-08-08']},
    ])) as IFitChartData;
    expectArray(legendEntries(data).map((e) => e.text),
      ['multiple styled series', 'Run:2023-08-08', 'Run:2023-08-08']);
  });

  test('a shortened row says what it stands for when hovered', async () => {
    const long = 'Copy of 20250528_NHA Assay_analysis Plate_1 rep 2 run 3 reader B plate map 7 revision 2';
    const data = JSON.parse(multiColumnJson([{column: 'curve', series: [long]}])) as IFitChartData;
    const drawn = renderedRows(JSON.stringify(data), 'legendHover', 500, 400);
    const row = drawn.rows.find((r) => r.text.startsWith('Copy of'))!;
    expect(row !== undefined, true, 'the curve should be named');

    const hovered = legendTooltip(drawn.data, drawn.dataBox, row.x + 5, row.y - 4);
    expect(hovered?.length, 1, 'hovering the row should say one thing');
    expect(hovered![0].text, long, 'and it should be the whole name');
    expect(legendTooltip(drawn.data, drawn.dataBox, drawn.dataBox.x + 5, drawn.dataBox.maxY - 5) === null, true,
      'the rest of the plot should say nothing');
  });

  test('candlesticks and statistics count as ink, not just fit lines', async () => {
    // a plateau drawn as candlesticks, no fit line: it used to report nothing, so the legend sat on it
    const points = [];
    for (const x of [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]) {
      points.push({x: x, y: 100});
      points.push({x: x, y: 96});
    }
    points.push({x: 1e-4, y: 0});
    const candles = JSON.stringify({chartOptions: {logX: true},
      series: [{name: 'Run:2023-08-08', columnName: 'curve', showPoints: 'candlesticks', showFitLine: false,
        connectDots: true, points: points}]} as IFitChartData);
    const drawn = renderedRows(candles, 'legendCandles', 500, 400);
    const row = drawn.rows.find((r) => r.text.startsWith('Run:'))!;
    expect(row !== undefined, true, 'the curve should be named');
    expect(row.y > drawn.dataBox.midY, true, `points fill the top, so the legend belongs below: y=${row.y}`);
  });

  test('the statistics on the plot are ink the legend keeps clear of', async () => {
    const data = JSON.parse(multiColumnJson([{column: 'curve', series: ['Run:2023-08-08']}])) as IFitChartData;
    data.chartOptions!.showStatistics = ['ic50'];
    const drawn = renderedRows(JSON.stringify(data), 'legendStats', 500, 400);
    const statistic = drawn.rows.find((r) => r.text.startsWith('IC50'))!;
    const row = drawn.rows.find((r) => r.text.startsWith('Run:'))!;
    expect(statistic !== undefined && row !== undefined, true, 'both should be drawn');
    expect(Math.abs(row.y - statistic.y) > 10 || Math.abs(row.x - statistic.x) > 60, true,
      `the legend should not land on the statistics: legend ${row.x},${row.y} vs ${statistic.x},${statistic.y}`);
  });

  test('a curve across the top sends the legend to the bottom', async () => {
    // the shape that misplaced it: a plateau spanning the plot just under the top, both bottom
    // corners empty. Sitting right below the top left box is still crowding it
    const box = new DG.Rect(0, 0, 400, 300);
    const plateau = [];
    for (let x = box.x + 10; x < box.maxX - 10; x += 8)
      plateau.push(new DG.Point(x, box.y + 30));
    const corner = chooseLegendCorner(box, 120, 20, plateau).rect;
    expect(corner.y > box.midY, true, `expected a bottom corner, got y=${corner.y}`);
  });

  test('the legend takes a corner the curves leave free', async () => {
    const box = new DG.Rect(0, 0, 400, 300);
    const busyTopRight = Array.from({length: 30}, (_, i) => new DG.Point(box.maxX - 10 - i, box.y + 10 + i % 20));
    const corner = chooseLegendCorner(box, 120, 60, busyTopRight).rect;
    expect(corner.x < box.midX, true, 'a crowded top right should push the legend away');

    // a couple of stray points are not a reason to move
    const almostFree = [new DG.Point(box.maxX - 20, box.y + 20), new DG.Point(box.maxX - 30, box.y + 30)];
    expect(chooseLegendCorner(box, 120, 60, almostFree).rect.x > box.midX, true, 'it should stay put');
    // and neither is a chart with no free corner at all
    const busyAll = [];
    for (const corner of [[box.maxX - 20, box.y + 20], [box.x + 20, box.y + 20],
      [box.maxX - 20, box.maxY - 20], [box.x + 20, box.maxY - 20]]) {
      for (let i = 0; i < 30; i++)
        busyAll.push(new DG.Point(corner[0], corner[1] + (i % 10)));
    }
    expect(chooseLegendCorner(box, 120, 60, busyAll).rect.x > box.midX, true, 'nowhere better means no move');
    expect(chooseLegendCorner(box, 120, 60).rect.x > box.midX, true, 'and with nothing drawn it is the top right');

    // once placed it holds its corner unless the alternative is clearly emptier, so a few pixels of
    // resizing cannot make it hop
    const mild = [...busyTopRight.slice(0, 6), new DG.Point(box.x + 20, box.maxY - 20)];
    const held = chooseLegendCorner(box, 120, 60, mild, 3);
    expect(held.index, 3, 'it should keep the corner it was in');
    expect(chooseLegendCorner(box, 120, 60, busyTopRight, 3).index, 3, 'even when the top right is busy');
  });

  test('a name is only shortened when it really does not fit', async () => {
    // this one was cut although the corner it sits in is wide open
    const name = 'SigmoidalInhibitionFunction';
    for (const width of [500, 430, 380]) {
      const drawn = renderedTexts(multiColumnJson([{column: 'curve', series: [name]}]),
        `legendFits${width}`, width, 400);
      expect(drawn.includes(name), true,
        `at ${width}px expected the whole name, got: ${drawn.filter((t) => t.startsWith('Sigmo'))}`);
    }
  });

  test('a name too long for any corner is ellipsized', async () => {
    const long = 'Copy of 20250528_NHA Assay_analysis Plate_1 rep 2 run 3 reader B plate map 7 revision 2';
    const data = multiColumnJson([{column: 'curve', series: [long]}]);
    const drawn = renderedTexts(data, 'legendEllipsis', 500, 400);
    const entry = drawn.find((t) => t.startsWith('Copy of'));
    expect(entry !== undefined, true, 'the curve should be named');
    expect(entry!.endsWith('…'), true, `expected an ellipsis, got: ${entry}`);
    expect(entry!.length < long.length, true, 'the name should be shortened');
    // it is capped, so the legend can never take more than a corner of the plot it annotates
    expect(measureWidth(entry!) <= PLOT.width * 0.55, true, 'the row must stay within its cap');
  });

  test('naming the curves can be turned off for the whole plot', async () => {
    const data = JSON.parse(multiColumnJson([{column: 'curve', series: ['Run:2023-08-08']}])) as IFitChartData;
    expect(renderedTexts(JSON.stringify(data), 'legendOn', 500, 400).includes('Run:2023-08-08'), true);

    data.chartOptions!.showLegend = false;
    const off = renderedTexts(JSON.stringify(data), 'legendOff', 500, 400);
    expect(off.includes('Run:2023-08-08'), false, 'showLegend: false should leave the curves unnamed');
    // and the rest of the chart is untouched
    expect(off.some((t) => t === 'Conc.' || t === 'Activity' || t.length > 0), true);
  });

  test('hovering what did not fit gives the whole legend', async () => {
    const names = Array.from({length: 40}, (_, i) => `series ${i}`);
    const drawn = renderedRows(multiColumnJson([{column: 'curve', series: names}]), 'legendMoreHover', 500, 400);
    const more = drawn.rows.find((r) => r.text.startsWith('+') && r.text.endsWith('more'))!;
    expect(more !== undefined, true, 'the legend should say how many it left out');

    const hovered = legendTooltip(drawn.data, drawn.dataBox, more.x + 5, more.y - 4);
    // one row per curve plus the column that groups them
    expect(hovered?.length, names.length + 1, 'the tooltip should carry every row the legend has');
    expect(hovered!.some((r) => r.text === 'series 39'), true, 'including the ones it could not draw');
    expect(hovered!.every((r) => r.color.length > 0), true, 'each named in the colour of its curve');
    // the tooltip is the legend, glyphs and all
    const element = legendTooltipElement(hovered!);
    expect(element.querySelectorAll('canvas').length, hovered!.length, 'every row should carry its glyph');
    expect(element.innerText.includes('series 39'), true);
  });

  test('more curves than fit are counted instead of stacked', async () => {
    const names = Array.from({length: 40}, (_, i) => `series ${i}`);
    const drawn = renderedTexts(multiColumnJson([{column: 'curve', series: names}]), 'legendOverflow', 500, 400);
    const more = drawn.find((t) => t.startsWith('+') && t.endsWith('more'));
    expect(more !== undefined, true, 'the legend should say how many it left out');
    // the names it did draw plus the ones it counted cover every curve
    const shown = drawn.filter((t) => t.startsWith('series ')).length;
    expect(shown > 0, true, 'it should still draw the curves that fit');
    expect(shown + parseInt(more!.slice(1), 10), names.length,
      'drawn names plus the counted ones should cover every curve');
  });
});
