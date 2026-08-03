import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

const SPIKE_RESULT_APPDATA_PATH = 'System:AppData/Proteomics/spike/spike-result.md';

function buildSpikeDataFrame(): DG.DataFrame {
  const baseMs = Date.UTC(2026, 5, 1);
  const oneDayMs = 24 * 60 * 60 * 1000;
  const dt: number[] = [];
  for (let i = 0; i < 5; i++)
    dt.push(baseMs + i * oneDayMs);
  const dtCol = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'acquisition_datetime', dt);
  const metricCol = DG.Column.fromFloat32Array('metric', new Float32Array([4.8, 5.1, 4.9, 5.2, 5.0]));
  return DG.DataFrame.fromColumns([dtCol, metricCol]);
}

function writeOutcome(body: string): void {
  try {
    grok.dapi.files.writeAsText(SPIKE_RESULT_APPDATA_PATH, body);
  } catch (_e) {
    // swallow — outcome is also logged to console + shell.info
  }
  try {
    grok.shell.info(body.split('\n')[0]);
  } catch (_e) { /* shell may not be addressable in headless test */ }
  // eslint-disable-next-line no-console
  console.log('[SPC spike] outcome\n' + body);
}

category('SPC', () => {
  test('SPC:spike_formula_line_dashed_on_datetime', async () => {
    const df = buildSpikeDataFrame();
    const viewer: any = DG.Viewer.lineChart(df, {
      xColumnName: 'acquisition_datetime',
      yColumnNames: ['metric'],
    });

    let outcome: 'PASS' | 'FAIL' = 'FAIL';
    let observation = '';
    let platformError = '';

    try {
      const fl: any = (viewer as any).meta?.formulaLines;
      if (!fl || typeof fl.addLine !== 'function') {
        platformError = 'viewer.meta.formulaLines.addLine not addressable on lineChart';
      } else {
        fl.addLine({
          formula: '${metric} = 5',
          style: 'dashed',
          color: '#666666',
          width: 1,
          visible: true,
          title: 'UCL',
        });
        const items = fl.items;
        const item0Style = items && items[0] ? items[0].style : undefined;
        observation = 'items[0].style=' + String(item0Style);
        if (item0Style === 'dashed')
          outcome = 'PASS';
        else
          platformError = 'style honored !== "dashed" (observed: ' + String(item0Style) + ')';
      }
    } catch (e: any) {
      platformError = (e && e.message) ? e.message : String(e);
    }

    const body =
      `OUTCOME: ${outcome}\n` +
      `Observation: ${observation}\n` +
      (platformError ? `PlatformError: ${platformError}\n` : '') +
      `Test: SPC:spike_formula_line_dashed_on_datetime\n` +
      `Timestamp: ${new Date().toISOString()}\n`;
    writeOutcome(body);

    // Spike does not gate. Either outcome is acceptable.
    expect(true, true);
  });

  test('SPC:spike_three_lines_ucl_cl_lcl', async () => {
    const df = buildSpikeDataFrame();
    const viewer: any = DG.Viewer.lineChart(df, {
      xColumnName: 'acquisition_datetime',
      yColumnNames: ['metric'],
    });

    const lines = [
      {formula: '${metric} = 5.4', style: 'dashed', color: '#666666', width: 1, visible: true, title: 'UCL'},
      {formula: '${metric} = 5.0', style: 'solid',  color: '#666666', width: 1, visible: true, title: 'CL'},
      {formula: '${metric} = 4.6', style: 'dashed', color: '#666666', width: 1, visible: true, title: 'LCL'},
    ];

    let added = 0;
    let lastError = '';
    let dashedHonoredCount = 0;
    try {
      const fl: any = (viewer as any).meta?.formulaLines;
      if (!fl || typeof fl.addLine !== 'function') {
        lastError = 'viewer.meta.formulaLines.addLine not addressable on lineChart';
      } else {
        for (const ln of lines) {
          try {
            fl.addLine(ln);
            added++;
          } catch (e: any) {
            lastError = (e && e.message) ? e.message : String(e);
          }
        }
        const items = fl.items ?? [];
        for (const it of items) {
          if (it && it.style === 'dashed') dashedHonoredCount++;
        }
      }
    } catch (e: any) {
      lastError = (e && e.message) ? e.message : String(e);
    }

    const body =
      `## Secondary observation: three-lines UCL/CL/LCL\n` +
      `Added: ${added}/3\n` +
      `DashedHonored: ${dashedHonoredCount}\n` +
      (lastError ? `PlatformError: ${lastError}\n` : '') +
      `Test: SPC:spike_three_lines_ucl_cl_lcl\n` +
      `Timestamp: ${new Date().toISOString()}\n`;
    writeOutcome(body);

    expect(true, true);
  });
});
