import {test, expect, Page} from '@playwright/test';

const BASE_URL = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const SPGI_PATH = 'System:AppData/Chem/tests/spgi-100.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

async function setup(page: Page) {
  await page.goto(BASE_URL);
  await page.waitForSelector('[name="Browse"], [name="Toolbox"]', {timeout: 120000});
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showHelp = false;
    grok.shell.closeAll();
  });
  await page.waitForTimeout(1000);
}

async function getFilteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

async function getTotalCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
}

test.describe('Filter Panel — Cloned Views', () => {
  test.setTimeout(600_000);

  test('Cloned view preserves filter state and layout save/restore works', async ({page}) => {
    stepErrors.length = 0;

    await setup(page);

    let totalRows = 0;

    // Step 1: Open SPGI dataset
    await softStep('1. Open SPGI dataset', async () => {
      await page.evaluate(async (path) => {
        const df = await grok.data.files.openTable(path);
        grok.shell.addTableView(df);
      }, SPGI_PATH);
      await page.waitForFunction(() => {
        const grid = grok.shell.tv?.grid;
        if (!grid) return false;
        const gc = grid.columns.byName('Structure');
        return gc && gc.cellType !== 'html' && gc.cellType !== 'string';
      }, {timeout: 60000});
      await page.waitForTimeout(4000);
      totalRows = await getTotalCount(page);
      expect(totalRows).toBeGreaterThan(0);
    });

    // Step 2: Open the Filter Panel
    await softStep('2. Open the Filter Panel', async () => {
      await page.evaluate(() => grok.shell.tv.getFiltersGroup());
      await page.waitForTimeout(3000);
      const hasFilters = await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'));
      expect(hasFilters).toBe(true);
    });

    // Step 3: Navigate to Competition assay filter
    await softStep('3. Navigate to Competition assay filter', async () => {
      const found = await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers) {
          if (h.textContent!.trim().startsWith('Competition assay-')) {
            if (!h.textContent!.trim().includes('Date')) {
              h.scrollIntoView({behavior: 'instant', block: 'center'});
              return true;
            }
          }
        }
        return false;
      });
      expect(found).toBe(true);
    });

    // Step 4: Filter out missing values on Competition assay
    let afterMissing = 0;
    await softStep('4. Filter out missing values for Competition assay', async () => {
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'Competition assay', filterOutMissingValues: true});
      });
      await page.waitForTimeout(2000);

      afterMissing = await getFilteredCount(page);
      expect(afterMissing).toBeLessThan(totalRows);
    });

    // Step 5: Set Stereo Category filter to S_ACHIR
    let afterStereo = 0;
    await softStep('5. Set Stereo Category filter to S_ACHIR', async () => {
      await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stereo Category', selected: ['S_ACHIR']});
      });
      await page.waitForTimeout(1000);
      afterStereo = await getFilteredCount(page);
      expect(afterStereo).toBeLessThan(afterMissing);
    });

    // Step 6: Sketch c1ccccc1 (benzene) in Structure filter
    let afterBenzene = 0;
    await softStep('6. Sketch c1ccccc1 structure', async () => {
      // Click Sketch link on Structure filter
      await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers) {
          if (h.textContent!.trim() === 'Structure') {
            const card = h.closest('.d4-filter')!;
            const sketchLink = card.querySelector('.sketch-link') as HTMLElement;
            if (sketchLink) sketchLink.click();
            break;
          }
        }
      });
      await page.waitForTimeout(2000);

      // Type SMILES and confirm
      const smilesInput = page.getByPlaceholder('SMILES, MOLBLOCK, Inchi, ChEMBL id, etc');
      await smilesInput.fill('c1ccccc1');
      await page.keyboard.press('Enter');
      await page.waitForTimeout(1000);
      await page.locator('button:has-text("OK")').click();
      await page.waitForTimeout(3000);

      afterBenzene = await getFilteredCount(page);
      expect(afterBenzene).toBeLessThan(afterStereo);
    });

    // Step 7: Clone View via View > Layout > Clone View
    await softStep('7. Clone View via top menu', async () => {
      await page.locator('text="View"').first().click();
      await page.waitForTimeout(500);
      await page.locator('text="Layout"').first().hover();
      await page.waitForTimeout(500);
      await page.locator('text="Clone View"').first().click();
      await page.waitForTimeout(3000);

      const viewCount = await page.evaluate(() => {
        let count = 0;
        for (const v of grok.shell.views)
          if (v.type === 'TableView') count++;
        return count;
      });
      expect(viewCount).toBe(2);
    });

    // Step 8: Check Filter Panel on clone matches original
    await softStep('8. Check cloned view filter state matches', async () => {
      const cloneFiltered = await getFilteredCount(page);
      expect(cloneFiltered).toBe(afterBenzene);

      const hasFilterPanel = await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'));
      expect(hasFilterPanel).toBe(true);
    });

    // Step 9: Turn all filters off
    await softStep('9. Turn all filters off (global checkbox)', async () => {
      await page.evaluate(() => {
        const filterViewer = document.querySelector('[name="viewer-Filters"]')!;
        const cb = filterViewer.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
        cb.click();
      });
      await page.waitForTimeout(1000);

      const allVisible = await getFilteredCount(page);
      expect(allVisible).toBe(totalRows);
    });

    // Step 10: Turn filters on again
    await softStep('10. Turn filters on again', async () => {
      await page.evaluate(() => {
        const filterViewer = document.querySelector('[name="viewer-Filters"]')!;
        const cb = filterViewer.querySelector('.d4-filter-group-header input[type="checkbox"]') as HTMLInputElement;
        cb.click();
      });
      await page.waitForTimeout(1000);

      const refiltered = await getFilteredCount(page);
      expect(refiltered).toBe(afterBenzene);
    });

    // Step 11: Clear Structure filter, set to C1CCCCC1
    let afterCyclohexane = 0;
    await softStep('11. Clear Structure filter, set to C1CCCCC1', async () => {
      // Click Clear on Structure filter
      await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers) {
          if (h.textContent!.trim() === 'Structure') {
            const card = h.closest('.d4-filter')!;
            const clearSpan = Array.from(card.querySelectorAll('span')).find(s => s.textContent!.trim() === 'Clear');
            if (clearSpan) clearSpan.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true}));
            break;
          }
        }
      });
      await page.waitForTimeout(2000);

      const afterClear = await getFilteredCount(page);
      expect(afterClear).toBe(afterStereo);

      // Now set C1CCCCC1
      await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers) {
          if (h.textContent!.trim() === 'Structure') {
            const card = h.closest('.d4-filter')!;
            const sketchLink = card.querySelector('.sketch-link') as HTMLElement;
            if (sketchLink) sketchLink.click();
            break;
          }
        }
      });
      await page.waitForTimeout(2000);

      const smilesInput = page.getByPlaceholder('SMILES, MOLBLOCK, Inchi, ChEMBL id, etc');
      await smilesInput.fill('C1CCCCC1');
      await page.keyboard.press('Enter');
      await page.waitForTimeout(1000);
      await page.locator('button:has-text("OK")').click();
      await page.waitForTimeout(3000);

      afterCyclohexane = await getFilteredCount(page);
      expect(afterCyclohexane).not.toBe(afterStereo);
    });

    // Step 12: Remove Structure filter from cloned view
    await softStep('12. Remove Structure filter (click X)', async () => {
      await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers) {
          if (h.textContent!.trim().startsWith('Structure')) {
            const card = h.closest('.d4-filter')!;
            const closeIcon = card.querySelector('[name="icon-times"]') as HTMLElement;
            if (closeIcon) {
              closeIcon.scrollIntoView({behavior: 'instant', block: 'center'});
              closeIcon.click();
            }
            break;
          }
        }
      });
      await page.waitForTimeout(1000);

      // Verify structure filter removed from clone
      const hasStructure = await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers)
          if (h.textContent!.trim() === 'Structure') return true;
        return false;
      });
      expect(hasStructure).toBe(false);

      // Filtered state should not change
      const filtered = await getFilteredCount(page);
      expect(filtered).toBe(afterCyclohexane);
    });

    // Step 13: Save the layout
    let savedLayout: string;
    await softStep('13. Save the layout', async () => {
      savedLayout = await page.evaluate(() => {
        const layout = grok.shell.tv.saveLayout();
        return JSON.stringify(layout.toJson());
      });
      expect(savedLayout).toBeTruthy();
    });

    // Step 14: Close the Filter Panel
    await softStep('14. Close the Filter Panel', async () => {
      await page.evaluate(() => {
        for (const v of grok.shell.tv.viewers)
          if (v.type === 'Filters') { v.close(); break; }
      });
      await page.waitForTimeout(1000);

      const filterPanel = await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'));
      expect(filterPanel).toBe(false);
    });

    // Step 15: Apply the saved layout
    await softStep('15. Apply saved layout — Filter Panel opens without Structure', async () => {
      await page.evaluate((json) => {
        const layout = DG.ViewLayout.fromJson(JSON.parse(json));
        grok.shell.tv.loadLayout(layout);
      }, savedLayout!);
      await page.waitForTimeout(5000);

      // Filter panel should be open
      const filterPanelOpen = await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'));
      expect(filterPanelOpen).toBe(true);

      // Structure filter should NOT be present
      const hasStructure = await page.evaluate(() => {
        const headers = document.querySelectorAll('.d4-filter-header');
        for (const h of headers)
          if (h.textContent!.trim() === 'Structure') return true;
        return false;
      });
      expect(hasStructure).toBe(false);

      // Filter state should be preserved
      const filtered = await getFilteredCount(page);
      expect(filtered).toBe(afterCyclohexane);
    });

    // Cleanup
    await page.evaluate(() => grok.shell.closeAll());

    if (stepErrors.length > 0) {
      const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
});
