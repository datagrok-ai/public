import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

/**
 * Section-local twin of `helpers/viewers.addLegendViewers`: same viewers, same legend column, same
 * `legendVisibility = 'Always'` on every non-Grid viewer, and the same column-prop map. The shared
 * helper pays a flat 300ms per viewer plus a `settleMs` sleep at the end — 3.6s on the seven-viewer
 * specs. Here the add waits on the viewer appearing and the settle waits on the legend item counts
 * holding still, capped at the sleep it replaces.
 */
export async function addLegendViewers(
  page: Page, options: {column: string; viewers: string[]; capMs?: number},
): Promise<void> {
  await v.installEventWaits(page);
  await page.evaluate(async ({types, col, cap, map}) => {
    const w = window as any;
    const tv = w.grok.shell.tv;
    for (const t of types) {
      const before = tv.viewers.filter((x: any) => x.type === t).length;
      tv.addViewer(t);
      await w.__poll(() => tv.viewers.filter((x: any) => x.type === t).length,
        (c: number) => c > before, 1500, 25);
    }
    for (const view of tv.viewers) {
      if (view.type === 'Grid') continue;
      try {
        const entry = types.includes(view.type) ? (map as any)[view.type] : null;
        if (entry) view.props[entry.prop] = entry.array ? [col] : col;
        try { view.props.legendVisibility = 'Always'; } catch (_) {}
      } catch (_) {}
    }
    const stamp = () => tv.viewers.filter((x: any) => x.type !== 'Grid')
      .map((x: any) => x.root.querySelectorAll('[name="legend"] .d4-legend-item').length).join(',');
    await w.__poll(stamp, (s: string) => s.split(',').some((n) => Number(n) > 0), cap, 25);
    await w.__settledFor(stamp, 150, cap, 25);
  }, {types: options.viewers, col: options.column, cap: options.capMs ?? 1500, map: v.LEGEND_COLUMN_PROP});
}

/**
 * Holds until the row filter count has been unchanged for a quiet gap, capped at `capMs`. For the
 * negative baselines: a gesture that must NOT filter cannot be proven by a poll that returns the
 * pre-gesture value immediately, so the hold is real — but it ends on the quiet gap, not the cap.
 */
export async function holdFilterCount(page: Page, capMs: number): Promise<number> {
  return page.evaluate((cap) => {
    const w = window as any;
    const df = w.grok.shell.tv.dataFrame;
    return w.__settledFor(() => df.filter.trueCount, 200, cap, 25);
  }, capMs);
}
