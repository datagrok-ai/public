import { test, expect } from './helpers/diff-studio';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import {
  openDiffStudio, openModelFromLibrary, toggleRibbonSwitch,
  ribbonSwitchOn, selectChoice, inputEditor, inputHost,
} from './helpers/diff-studio';

test('DiffStudio Sensitivity Analysis — Bioreactor: open SA view, Process mode cascade, Run produces 4 viewers',
  async ({ page }) => {
    test.setTimeout(300_000);
    const { softStep, assertAllPassed } = createSoftStepCollector();
    const monitor = attachErrorMonitor(page);

    await softStep('Step 1: Open Diff Studio + Bioreactor; Edit OFF', async () => {
      await openDiffStudio(page);
      await openModelFromLibrary(page, 'Bioreactor');
      if (await ribbonSwitchOn(page, 'Edit')) await toggleRibbonSwitch(page, 'Edit');
      await expect(page.locator(inputHost('Process-mode'))).toBeVisible({ timeout: 10_000 });
    });

    await softStep('Step 2: Click Sensitivity icon — SA view opens', async () => {

      await page.locator('.diff-studio-ribbon-sa-icon').first().click();

      await page.waitForTimeout(2000);
      await expect(page.locator(inputHost('Process-mode'))).toBeVisible({ timeout: 15_000 });
    });

    await softStep('Step 3: Modify Process mode; FFox/KKox/others cascade-update', async () => {
      const readAll = async () => ({
        ffox: await page.locator(inputEditor('FFox')).inputValue().catch(() => ''),
        kkox: await page.locator(inputEditor('KKox')).inputValue().catch(() => ''),
        mea: await page.locator(inputEditor('MEAthiol')).inputValue().catch(() => ''),
        gas: await page.locator(inputEditor('Gas')).inputValue().catch(() => ''),
      });
      const before = await readAll();
      await selectChoice(page, 'Process-mode', 'Mode 1');
      await page.waitForTimeout(2500);
      const after = await readAll();
      const changed = Object.keys(before).filter(k => (before as any)[k] !== (after as any)[k]);
      expect(changed.length).toBeGreaterThanOrEqual(3);
    });

    await softStep('Step 4: Enable switchers and Run — four viewers open with valid data', async () => {

      const enableSwitcher = async (safeName: string): Promise<boolean> => {
        return await page.evaluate((name) => {
          const host = document.querySelector(`[name="input-host-${name}"]`) as HTMLElement | null;
          if (!host) return false;
          const sw = host.querySelector('.sa-switch-input .ui-input-switch') as HTMLElement | null;
          if (!sw) return false;
          if (sw.classList.contains('ui-input-switch-on')) return true;
          sw.scrollIntoView({ block: 'center' });
          sw.click();
          return true;
        }, safeName);
      };

      await enableSwitcher('FFox');
      await enableSwitcher('FKox');
      await enableSwitcher('FFred');
      await page.waitForTimeout(800);

      await page.locator('.d4-ribbon-item i.fa-play').first().click();

      await page.waitForTimeout(15_000);

      const viewerCount = await page.locator('.d4-viewer').count();
      expect(viewerCount).toBeGreaterThanOrEqual(4);

      const viewerInfo = await page.evaluate(() => {
        const win = window as any;
        const view = win.grok?.shell?.v;
        if (!view) return null;
        const viewers = Array.from(view.viewers ?? []) as any[];
        if (viewers.length === 0) return null;
        return viewers.map((v) => {
          const df = v.dataFrame ?? v.table;
          let nonFiniteCount = 0;
          let totalCells = 0;
          if (df?.columns?.toList) {
            for (const col of df.columns.toList()) {
              const type: string = col.type;
              if (type !== 'double' && type !== 'float' && type !== 'int' && type !== 'qnum')
                continue;
              const len: number = col.length ?? 0;
              for (let i = 0; i < len; i++) {
                totalCells++;
                const x = col.get(i);
                if (typeof x === 'number' && !isFinite(x)) nonFiniteCount++;
              }
            }
          }
          return {
            type: v.type ?? '?',
            name: v.name ?? v.title ?? '',
            rowCount: df?.rowCount ?? 0,
            columnCount: df?.columns?.length ?? 0,
            nonFiniteCount,
            totalCells,
          };
        });
      });

      expect(viewerInfo).not.toBeNull();
      expect(viewerInfo!.length).toBeGreaterThanOrEqual(4);

      const viewersWithData = viewerInfo!.filter(v => v.rowCount > 0 && v.columnCount > 0);
      expect(viewersWithData.length).toBeGreaterThanOrEqual(1);
      for (const v of viewerInfo!) {
        expect(v.nonFiniteCount).toBe(0);

        expect(v.name).not.toMatch(/^untitled$/i);
      }
    });

    assertAllPassed();
    monitor.assertNone();
  });
