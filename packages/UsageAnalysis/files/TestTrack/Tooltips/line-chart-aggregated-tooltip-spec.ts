import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;

async function readTooltipAtChartPoint(page: import('@playwright/test').Page,
                                       fracX: number, fracY: number)
                                       : Promise<{display: string; text: string}> {
  const lc = page.locator('[name="viewer-Line-chart"] canvas').first();
  const box = await lc.boundingBox();
  if (!box) throw new Error('Line chart canvas not visible');
  await page.evaluate(() => (window as any).ui?.tooltip?.hide?.());
  await page.mouse.move(0, 0);
  await page.waitForTimeout(200);
  // Move off the chart, then onto the target point.
  await page.mouse.move(box.x - 50, box.y + box.height / 2);
  await page.waitForTimeout(150);
  const targetX = box.x + box.width * fracX;
  const targetY = box.y + box.height * fracY;
  await page.mouse.move(targetX, targetY, {steps: 10});
  await page.waitForTimeout(1000);
  return page.evaluate(() => {
    const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
    if (!tt) return {display: 'missing', text: ''};
    return {display: getComputedStyle(tt).display, text: tt.textContent ?? ''};
  });
}

test('Line chart: aggregated tooltip with split column', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Setup + open SPGI dataset (scenario's "Open SPGI Dataset" step)
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_: any, i: number) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Add Line Chart and set X=Chemist 521, Y=CAST Idea ID', async () => {
    await page.evaluate(() => {
      (document.querySelector('[name="icon-line-chart"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 10000});
    await page.waitForTimeout(800);
    await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      lc.props.xColumnName = 'Chemist 521';
      lc.props.yColumnNames = ['CAST Idea ID'];
    });
    await page.waitForTimeout(400);
    const cfg = await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      return {x: lc.props.xColumnName, y: lc.props.yColumnNames};
    });
    expect(cfg.x).toBe('Chemist 521');
    expect(cfg.y).toEqual(['CAST Idea ID']);
  });

  await softStep('Open tooltip settings and add aggregated tooltip entries', async () => {
    // Open Edit Aggregated Tooltip dialog via context menu Tooltip > Edit...
    await page.evaluate(async () => {
      const lc = document.querySelector('[name="viewer-Line-chart"]')!;
      const canvas = lc.querySelector('canvas')!;
      const r = canvas.getBoundingClientRect();
      ['mousedown', 'mouseup', 'contextmenu'].forEach((t) => {
        canvas.dispatchEvent(new MouseEvent(t, {
          bubbles: true, cancelable: true,
          clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, button: 2,
        }));
      });
      await new Promise((r2) => setTimeout(r2, 400));
      const edit = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item'))
        .find((it: any) => it.querySelector('.d4-menu-item-label')?.textContent.trim() === 'Edit...');
      (edit as HTMLElement)?.click();
    });
    await page.waitForTimeout(600);
    // Cancel the dialog first — the column-picker inside it isn't reachable through synthetic
    // events, and leaving the dialog open would commit an empty state on OK. Set the property
    // directly, then re-open the dialog to verify the value parsed into rows correctly.
    await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d: any) => d.offsetWidth > 0 && d.querySelector('.d4-dialog-header')?.textContent.includes('Edit Aggregated Tooltip'));
      (dialog?.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click();
    });
    await page.waitForTimeout(300);
    await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      lc.setOptions({aggTooltipColumns: 'unique(Stereo Category)\nmin(Average Mass)'});
    });
    await page.waitForTimeout(300);
    // Re-open the dialog and read back the parsed rows.
    await page.evaluate(async () => {
      const lc = document.querySelector('[name="viewer-Line-chart"]')!;
      const canvas = lc.querySelector('canvas')!;
      const r = canvas.getBoundingClientRect();
      ['mousedown', 'mouseup', 'contextmenu'].forEach((t) => {
        canvas.dispatchEvent(new MouseEvent(t, {
          bubbles: true, cancelable: true,
          clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, button: 2,
        }));
      });
      await new Promise((r2) => setTimeout(r2, 400));
      const edit = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item'))
        .find((it: any) => it.querySelector('.d4-menu-item-label')?.textContent.trim() === 'Edit...');
      (edit as HTMLElement)?.click();
    });
    await page.waitForTimeout(600);
    const rows = await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d: any) => d.offsetWidth > 0 && d.querySelector('.d4-dialog-header')?.textContent.includes('Edit Aggregated Tooltip'));
      if (!dialog) return [];
      return Array.from(dialog.querySelectorAll('tbody tr'))
        .filter((tr: any) => tr.querySelector('.d4-column-selector'))
        .map((tr: any) => ({
          column: tr.querySelector('.d4-column-selector-column')?.textContent.trim(),
          aggr: tr.querySelector('select.ui-input-editor')?.value,
        }));
    });
    expect(rows).toEqual([
      {column: 'Stereo Category', aggr: 'unique'},
      {column: 'Average Mass', aggr: 'min'},
    ]);
    // Close the dialog with CANCEL to avoid OK overwriting state.
    await page.evaluate(() => {
      const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
        .find((d: any) => d.offsetWidth > 0 && d.querySelector('.d4-dialog-header')?.textContent.includes('Edit Aggregated Tooltip'));
      (dialog?.querySelector('[name="button-CANCEL"]') as HTMLElement | null)?.click();
    });
    await page.waitForTimeout(300);
  });

  await softStep('Split chart by Stereo Category', async () => {
    await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      lc.props.splitColumnNames = ['Stereo Category'];
    });
    await page.waitForTimeout(600);
    const split = await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      return lc.props.splitColumnNames;
    });
    expect(split).toEqual(['Stereo Category']);
  });

  await softStep('Hover over Line Chart points', async () => {
    // Probe many positions across the chart — different X bins / Y series surface different
    // tooltip content. Aggregated tooltips with the configured columns only appear when the
    // hover lands on a specific marker hit-test region.
    let anyVisible = false;
    let aggregatedText = '';
    let probesSeen: string[] = [];
    for (let fx = 0.1; fx <= 0.9; fx += 0.05) {
      for (const fy of [0.25, 0.4, 0.55, 0.7]) {
        const tt = await readTooltipAtChartPoint(page, fx, fy);
        if (tt.display === 'block') anyVisible = true;
        probesSeen.push(`${fx.toFixed(2)}x${fy} ${tt.display}: ${tt.text.slice(0, 80)}`);
        if (/unique\(Stereo Category\)/.test(tt.text) && /min\(Average Mass\)/.test(tt.text)) {
          aggregatedText = tt.text;
          break;
        }
      }
      if (aggregatedText) break;
    }
    (test.info() as any).aggregatedText = aggregatedText;
    (test.info() as any).probesSeen = probesSeen;
    expect(anyVisible).toBe(true);
  });

  await softStep('Verify tooltip displays configured aggregated info, no errors', async () => {
    const aggregatedText: string = (test.info() as any).aggregatedText ?? '';
    if (!aggregatedText) {
      const probes: string[] = (test.info() as any).probesSeen ?? [];
      console.log('Probe trace (no aggregated tooltip surfaced):\n' + probes.join('\n'));
    }
    expect(aggregatedText).toContain('unique(Stereo Category)');
    expect(aggregatedText).toContain('min(Average Mass)');
    // No chart-related errors raised via grok.shell.
    const warnings = await page.evaluate(() => grok.shell.warnings?.length ?? 0);
    expect(warnings).toBe(0);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
