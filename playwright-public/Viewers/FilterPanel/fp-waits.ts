import {Page} from '@playwright/test';

// __moved spends its whole cap whenever a gesture legitimately changes nothing, which is most of a
// menu ladder ("Select all" over an already-full card, a mode switch that must not move the rows, a
// probe click that misses the bar). This adds the one signal that says nothing more is coming: the
// frame's own ddt-rows-filtered event, armed before the act. A pass that does land is still awaited
// in full and then settled — FiltersCore debounces its criteria at 50ms (filters_core.dart:278) —
// and a gesture that raises no pass at all stops at `graceMs` instead of at the cap.
//
// __fpQuick runs the act itself, for a gesture dispatched inside the page. __fpArm / __fpWait are
// the same wait split in two, for a gesture that Playwright drives from Node between them.
export async function installFilterWaits(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__fpWait) return;

    const settle = async (read: () => any, from: any, t0: number, capMs: number, graceMs: number,
      firedAt: () => number) => {
      for (;;) {
        if (read() !== from)
          return w.__settledFor(read, 150, Math.max(0, t0 + capMs - Date.now()), 25);
        const now = Date.now();
        if (now - t0 >= capMs || now - (firedAt() || t0) > graceMs) return read();
        await new Promise((r) => setTimeout(r, 25));
      }
    };

    w.__fpArm = () => {
      w.__fpSub?.unsubscribe();
      w.__fpFired = 0;
      w.__fpSub = w.grok.shell.tv.dataFrame.onRowsFiltered.subscribe(() => { w.__fpFired = Date.now(); });
    };

    w.__fpWait = async (read: () => any, from: any, capMs = 1200, graceMs = 300) => {
      const out = await settle(read, from, Date.now(), capMs, graceMs, () => w.__fpFired ?? 0);
      w.__fpSub?.unsubscribe();
      w.__fpSub = null;
      return out;
    };

    w.__fpQuick = async (read: () => any, act: () => any, capMs = 1200, graceMs = 300) => {
      const df = w.grok.shell.tv.dataFrame;
      const from = read();
      let firedAt = 0;
      const sub = df.onRowsFiltered.subscribe(() => { firedAt = Date.now(); });
      const t0 = Date.now();
      try {
        await act();
        return await settle(read, from, t0, capMs, graceMs, () => firedAt);
      }
      finally { sub.unsubscribe(); }
    };
  });
}
