import { test, expect } from '@playwright/test';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import {
  BASE, openDiffStudio, openModelFromLibrary, setInputValue, inputEditor, inputHost,
} from './helpers/diff-studio';

/**
 * Test Track scenario: DiffStudio/catalog.md
 * 1. Open a model from the library (Open model icon → Library → PK-PD).
 * 2. Save the model to the Model Hub catalog (Save to Model Hub icon).
 * 3. Access the Model Hub (Apps → Model Hub).
 * 4. Click the Refresh icon — the model list refreshes.
 * 5. Run the PK-PD model from the Model Hub catalog.
 * 6. Modify input parameters and verify the results update.
 */
test('DiffStudio Catalog — PK-PD: load → Save to Model Hub → refresh → run from catalog → modify dose',
  async ({ page }) => {
    test.setTimeout(300_000);
    const { softStep, assertAllPassed } = createSoftStepCollector();
    const monitor = attachErrorMonitor(page);

    await softStep('Step 1: Open PK-PD from Diff Studio Library', async () => {
      await openDiffStudio(page);
      await openModelFromLibrary(page, 'PK-PD');
      await expect(page.locator(inputHost('dose'))).toBeVisible({ timeout: 15_000 });
    });

    await softStep('Step 2: Click "Save to Model Hub" icon', async () => {
      await page.locator('.diff-studio-ribbon-save-to-model-catalog-icon').first().click();
      // The platform shows a balloon notification on success
      const balloon = page.locator('.d4-balloon, .grok-notification').filter({
        hasText: /Saved to Library|PK-PD\.ivp/,
      });
      await expect(balloon.first()).toBeVisible({ timeout: 15_000 });
    });

    await softStep('Step 3: Open the Model Hub (Apps → Compute → Model Hub) — PK-PD is listed', async () => {
      // Walk the Browse tree: Apps → Compute → Model Hub. URL `/apps/Compute2/ModelCatalog`
      // returns "Application not found"; the tree navigation is the canonical UI path.
      await page.goto(BASE);
      await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
      await page.waitForTimeout(1500);

      const clickTreeLabel = async (label: string): Promise<boolean> => {
        return await page.evaluate((text) => {
          const candidates = Array.from(document.querySelectorAll(
            '.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item-label')) as HTMLElement[];
          const el = candidates.find(e => e.textContent?.trim() === text);
          if (!el) return false;
          el.scrollIntoView({ behavior: 'instant', block: 'center' });
          el.click();
          return true;
        }, label);
      };

      expect(await clickTreeLabel('Apps')).toBe(true);
      await page.waitForTimeout(800);
      expect(await clickTreeLabel('Compute')).toBe(true);
      await page.waitForTimeout(800);
      expect(await clickTreeLabel('Model Hub')).toBe(true);
      await page.waitForTimeout(3000);

      // Model Hub renders saved models as cards labelled with the model name.
      await expect(page.getByText('PK-PD', { exact: true }).first()).toBeVisible({ timeout: 20_000 });
    });

    await softStep('Step 4: Refresh actually reloads — count changes when files appear/disappear', async () => {
      const refresh = page.locator(
        '.d4-ribbon i.fa-sync, .d4-ribbon i.fa-sync-alt, .d4-ribbon i.fa-refresh, .d4-ribbon i.fa-redo',
      ).first();

      const cardCount = async (): Promise<number> =>
        await page.getByText('PK-PD', { exact: true }).count();

      const waitForCatalogReady = async (): Promise<void> => {
        await page.waitForFunction(() => {
          return document.body.innerText.includes('PK-PD');
        }, null, { timeout: 15_000 });
        await page.waitForTimeout(1500);
      };

      // Baseline: refresh once, wait for catalog to repopulate, count PK-PD cards.
      if (await refresh.count() > 0) await refresh.click({ force: true });
      await waitForCatalogReady();
      const before = await cardCount();
      expect(before).toBeGreaterThan(0);

      // Delete the freshly-saved user copy from Step 2. The platform stores user catalog entries
      // as scripts (tag `model`); deletion via the file API removes the underlying .ivp written
      // by `Save to Model Hub`. UI alternative: right-click → Delete on the card. The API is
      // faster and deterministic for asserting "Refresh observably reloads".
      const deleted = await page.evaluate(async () => {
        try {
          const grok = (window as any).grok;
          // Find user-saved PK-PD models. Names vary (e.g. PK-PD copies); match by name.
          const scripts = await grok.dapi.scripts.filter('name = "PK-PD"').list();
          if (!scripts || scripts.length === 0) return 0;
          // Delete just the latest one — the test only needs ONE delta to prove refresh works.
          const last = scripts[scripts.length - 1];
          await grok.dapi.scripts.delete(last);
          return 1;
        } catch {
          return 0;
        }
      });

      if (deleted > 0) {
        await refresh.click({ force: true });
        await waitForCatalogReady();
        const after = await cardCount();
        expect(after).toBe(before - deleted);
      } else {
        // Couldn't delete (no user-saved PK-PD; e.g. the file save in Step 2 has not yet
        // materialised). Fall back to the weaker check: PK-PD still listed after refresh.
        await expect(page.getByText('PK-PD', { exact: true }).first())
          .toBeVisible({ timeout: 15_000 });
      }
    });

    await softStep('Step 5: Run the PK-PD model from the Model Hub catalog', async () => {
      // Pick the last "PK-PD" entry — the freshly saved user copy
      const items = page.getByText('PK-PD', { exact: true });
      const count = await items.count();
      expect(count).toBeGreaterThan(0);
      await items.nth(count - 1).dblclick();
      await page.locator(inputHost('dose')).waitFor({ timeout: 30_000 });
      await expect(page.locator(inputHost('count'))).toBeVisible({ timeout: 10_000 });
    });

    await softStep('Step 6: Modify dose input; value updates live', async () => {
      // Models opened via Model Hub run inside a Compute2 RichFunctionView. Its URL does NOT
      // mirror inputs and its chart is wrapped in Vue components that hide the canvas from
      // every selector strategy tried (.d4-viewer, .echarts, .rfv-chart, broad screenshots).
      // Chart redraw verification therefore stays manual — see ui-only.md M-1.7.
      const ed = page.locator(inputEditor('dose'));
      const before = await ed.inputValue();
      await setInputValue(page, 'dose', '5000');
      const after = await ed.inputValue();
      expect(after).toBe('5000');
      expect(after).not.toBe(before);
    });

    assertAllPassed();
    monitor.assertNone();
  });
