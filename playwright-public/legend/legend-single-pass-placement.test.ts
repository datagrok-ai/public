/* ---
sub_features_covered: [legend.placement.single-pass, legend.placement.no-settle-flash]
--- */
// A resize decides the legend placement at layout time, in the same render pass — the
// viewer lays out against the final legend inset once, instead of committing an
// intermediate state (the old slot re-estimated at the new size) and moving the legend on a
// delayed settle pass, which made hosts like the trellis re-render and visibly jump twice.
// The occupancy marks keep their chart-box fraction across a resize, so the corner decision
// extrapolates from last frame's grid instead of waiting for a repaint.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '@datagrok-libraries/test/src/playwright/spec-login';
import * as v from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

async function recordPlacements(page: any) {
  await page.evaluate(() => {
    (window as any).__legendMoves = [];
    if ((window as any).__legendSub) return;
    (window as any).__legendSub = (window as any).grok.events
      .onEvent('d4-legend-placement-changed')
      .subscribe((a: any) => (window as any).__legendMoves.push(
        `${a.args.width}x${a.args.height} ${a.args.fromMode}/${a.args.fromSlot}->${a.args.toMode}/${a.args.toSlot}`));
  });
}

const moves = (page: any): Promise<string[]> =>
  page.evaluate(() => (window as any).__legendMoves);

/** Samples the legend root's slot and box every animation frame; keeps distinct states. */
async function startSampler(page: any, viewerType: string) {
  await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    const samples: any[] = (window as any).__legendSamples = [];
    (window as any).__samplerStop = false;
    const tick = () => {
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      if (root) {
        const r = root.getBoundingClientRect();
        const s = {slot: root.getAttribute('data-legend-slot'), mode: root.getAttribute('data-legend-mode'),
          passes: Number(root.getAttribute('data-legend-passes')),
          w: Math.round(r.width), h: Math.round(r.height)};
        const last = samples[samples.length - 1];
        if (!last || last.slot !== s.slot || last.mode !== s.mode || last.passes !== s.passes
            || last.w !== s.w || last.h !== s.h)
          samples.push(s);
      }
      if (!(window as any).__samplerStop) requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  }, viewerType);
}

const stopSampler = (page: any): Promise<any[]> =>
  page.evaluate(() => { (window as any).__samplerStop = true; return (window as any).__legendSamples; });

async function legendState(page: any, viewerType: string) {
  return await page.evaluate((vt: string) => {
    const tv = (window as any).grok.shell.tv;
    const x = tv.viewers.find((q: any) => q.type === vt);
    const root = x?.root.querySelector('[name="legend"]') as HTMLElement;
    if (!root) return {noLegend: true};
    const r = root.getBoundingClientRect();
    return {mode: root.getAttribute('data-legend-mode'), slot: root.getAttribute('data-legend-slot'),
      passes: Number(root.getAttribute('data-legend-passes')),
      w: Math.round(r.width), h: Math.round(r.height)};
  }, viewerType);
}

/** The single-pass invariant for any viewer: a resize produces at most one placement move,
 * never revisits a slot, and never re-commits the old slot at the new geometry before
 * moving (the settle-pass artifact). */
async function expectSinglePassResize(page: any, viewerType: string,
  from: [number, number], to: [number, number], label: string) {
  await v.resizeViewer(page, viewerType, from[0], from[1]);
  const initial = await legendState(page, viewerType);
  await recordPlacements(page);
  await page.evaluate(() => { (window as any).__legendMoves = []; });
  await startSampler(page, viewerType);
  await v.resizeViewer(page, viewerType, to[0], to[1]);
  const samples = await stopSampler(page);
  const m = await moves(page);

  expect(m.length, `${label}: expected at most one placement move:\n${m.join('\n')}`)
    .toBeLessThanOrEqual(1);

  const slotSeq = samples.map((s: any) => `${s.mode}/${s.slot}`)
    .filter((s: string, i: number, a: string[]) => i === 0 || s !== a[i - 1]);
  expect(new Set(slotSeq).size, `${label}: a slot was revisited: ${slotSeq.join(' -> ')}`)
    .toBe(slotSeq.length);

  // the flex layout stretches the old docked legend the instant the viewer resizes — that is
  // pure CSS, not a decision; a *commit* of the old slot after the resize (passes bumped)
  // is the settle-pass artifact this spec pins
  const final = samples[samples.length - 1];
  if (final && final.slot !== initial.slot) {
    const recommitted = samples.filter((s: any) => s.slot === initial.slot && s.passes > initial.passes);
    expect(recommitted.length,
      `${label}: the old slot was re-committed after the resize before the move:\n${JSON.stringify(samples)}`)
      .toBe(0);
  }
}

test('Trellis resize decides the legend slot in one pass (SPGI)', async ({page}) => {
  test.setTimeout(900_000);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(String(e).slice(0, 300)));

  await loginToDatagrok(page);
  await v.openTable(page);

  await page.evaluate(() => {
    (window as any).grok.shell.tv.addViewer('Trellis plot');
  });
  await v.waitForLegendIdle(page, 'Trellis plot');

  await softStep('A tall trellis docks its legend on top', async () => {
    await v.resizeViewer(page, 'Trellis plot', 460, 640);
    const s = await legendState(page, 'Trellis plot');
    expect(s.mode, JSON.stringify(s)).toBe('docked');
    expect(s.slot, JSON.stringify(s)).toBe('top');
  });

  await softStep('Widening commits the legend straight to the right — no full-width top state in between', async () => {
    await recordPlacements(page);
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const x = tv.viewers.find((q: any) => q.type === 'Trellis plot');
      const root = x.root.querySelector('[name="legend"]') as HTMLElement;
      (window as any).__p0 = Number(root.getAttribute('data-legend-passes'));
    });
    await startSampler(page, 'Trellis plot');
    await v.resizeViewer(page, 'Trellis plot', 1000, 640);
    const samples = await stopSampler(page);
    const m = await moves(page);

    const s = await legendState(page, 'Trellis plot');
    expect(s.slot, `final state: ${JSON.stringify(s)}`).toBe('right');

    expect(m.length, `expected a single placement move:\n${m.join('\n')}`).toBe(1);
    expect(m[0], m[0]).toContain('docked/top->docked/right');

    // the old settle pass first re-committed the top legend at the new width, then moved it;
    // the instant CSS stretch of the old legend (same passes count) is not a commit
    const p0 = await page.evaluate(() => (window as any).__p0);
    const wideTop = samples.filter((x: any) => x.slot === 'top' && x.passes > p0);
    expect(wideTop.length,
      `legend was re-committed docked-top after the resize before moving right:\n${JSON.stringify(samples)} p0=${p0}`).toBe(0);

    // once it lands on the right, it must not revisit the top
    const rightAt = samples.findIndex((x: any) => x.slot === 'right');
    const topAfter = samples.slice(rightAt + 1).filter((x: any) => x.slot === 'top');
    expect(topAfter.length, JSON.stringify(samples)).toBe(0);
  });

  await softStep('Other viewers resize in one pass too', async () => {
    await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      tv.viewers.find((q: any) => q.type === 'Trellis plot')?.close();
    });

    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Line chart',
        {yColumnNames: ['Average Mass'], splitColumnNames: ['Series']});
    });
    await v.waitForLegendIdle(page, 'Line chart');
    await expectSinglePassResize(page, 'Line chart', [460, 640], [1000, 460], 'line chart');
    await page.evaluate(() => {
      (window as any).grok.shell.tv.viewers.find((q: any) => q.type === 'Line chart')?.close();
    });

    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Scatter plot', {colorColumnName: 'Series'});
    });
    await v.waitForLegendIdle(page, 'Scatter plot');
    await expectSinglePassResize(page, 'Scatter plot', [600, 640], [1100, 500], 'scatter plot');
    await page.evaluate(() => {
      (window as any).grok.shell.tv.viewers.find((q: any) => q.type === 'Scatter plot')?.close();
    });

    await page.evaluate(() => {
      (window as any).grok.shell.tv.addViewer('Bar chart', {stackColumnName: 'Series'});
    });
    await v.waitForLegendIdle(page, 'Bar chart');
    await expectSinglePassResize(page, 'Bar chart', [600, 640], [1100, 500], 'bar chart');
  });

  expect(errors, `page errors:\n${errors.join('\n')}`).toEqual([]);
  await v.cleanupShell(page);
  v.finishSpec('Single-pass placement failures');
});
